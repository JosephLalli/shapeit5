/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <phaser/phaser_header.h>

#include <io/genotype_reader/genotype_reader_header.h>
#include <io/haplotype_writer.h>
#include <io/gmap_reader.h>
#include <io/pedigree_reader.h>
#include <io/haploid_reader.h>
#include <modules/genotype_builder.h>
#include <objects/super_site_builder.h>
#include <objects/supersite_debug.h>

using namespace std;

void phaser::read_files_and_initialise() {
	//step0: Initialize seed and multi-threading
	rng.setSeed(options["seed"].as < int > ());
	if (options["thread"].as < int > () > 1) {
		i_workers = 0; i_jobs = 0;
		id_workers = vector < pthread_t > (options["thread"].as < int > ());
		pthread_mutex_init(&mutex_workers, NULL);
	}

	//step1: Set up the genotype reader
	vrb.title("Reading genotype data:");
	genotype_reader readerG(H, G, V);
	readerG.setThreads(options["thread"].as < int > ());
	readerG.setRegion(options["region"].as < string > ());
	readerG.setMainFilename(options["input"].as < string > ());
	if (options.count("reference")) readerG.addReferenceFilename(options["reference"].as < string > ());
	if (options.count("scaffold")) readerG.addScaffoldFilename(options["scaffold"].as < string > ());
	if (options.count("filter-snp")) readerG.setFilterSNP();
	if (!options["filter-maf"].defaulted()) readerG.setFilterMAF(options["filter-maf"].as < double > ());

	//step2: Read the genotype data
	readerG.scanGenotypes();
	readerG.allocateGenotypes();
	readerG.readGenotypes();

	//step3: Read haploid samples
	if (options.count("haploids")) {
		haploid_reader readerH;
		readerH.readHaploidFile(options["haploids"].as < string > ());
		G.resetHaploidHeterozgotes(readerH.samples);
	}

	//step4: Read pedigrees and scaffold diploid kids
	if (options.count("pedigree")) {
		pedigree_reader readerP;
		readerP.readPedigreeFile(options["pedigree"].as < string > ());
		G.scaffoldUsingPedigrees(readerP);
	}

	//step5: Read and initialise genetic map
	vrb.title("Setting up genetic map:");
	if (options.count("map")) {
		gmap_reader readerGM;
		readerGM.readGeneticMapFile(options["map"].as < string > ());
		V.setGeneticMap(readerGM);
	} else V.setGeneticMap();
	M.initialise(V, options["hmm-ne"].as < int > (), (readerG.n_main_samples+readerG.n_ref_samples)*2);

	//step6: Initialize haplotype set
	vrb.title("Initializing data structures:");
	G.imputeMonomorphic(V);
	H.updateHaplotypes(G, true);
	H.transposeHaplotypes_H2V(true);
	// Seed per-sample RNGs deterministically so multithreaded runs are reproducible
	G.seedRngs(rng.getSeed());

	//step7: Initialize PBWT for selecting states
	if (pbwt_auto) {
		unsigned int cumulative_sample_size = readerG.n_main_samples + readerG.n_ref_samples;
		pbwt_depth = max(min((int)round(9-log10(cumulative_sample_size)), 8), 2);
		pbwt_modulo = max(min((log(cumulative_sample_size) - log(50) + 1) * 0.01, 0.15), 0.005);
		vrb.bullet("PBWT parameters auto setting : [modulo = " + stb.str(pbwt_modulo, 3) + " / depth = " + stb.str(pbwt_depth, 3) + "]");
	} else {
		pbwt_depth = options["pbwt-depth"].as < int > ();
		pbwt_modulo = options["pbwt-modulo"].as < double > ();
	}

	const int pbwt_mac = options["pbwt-mac"].as < int > ();
	supersite_mac_threshold = std::max(0, pbwt_mac);

	H.initialize(V,	pbwt_modulo,
					options["pbwt-window"].as < double > (),
					options["pbwt-mdr"].as < double > (),
					pbwt_depth,
					pbwt_mac,
					options["thread"].as < int > ());

	if (!options.count("pbwt-disable-init")) H.solve(&G);

	//step8: Build super-sites for multiallelic positions (one-time, before genotype building)
	// This collapses split biallelic records at identical (chr,bp) positions into
	// super-sites and packs per-haplotype codes (2 codes per byte).
	// Uses global panel order (H.n_hap), independent of per-window conditioning sets.
	// MUST be done before genotype_builder so build() can use locus_to_super_idx for segment boundaries.
		if (enable_supersites) {
			vrb.title("Building super-sites");
			rebuildSupersiteMetadata("initialise");
			
			vrb.bullet("Built " + stb.str(super_sites.size()) + " super-sites covering " + 
			           stb.str(std::count(is_super_site.begin(), is_super_site.end(), true)) + " variant positions");
			
			// Set supersite context for all genotypes BEFORE building segments
			const auto ss_cfg = supersite_invariants::SupersiteDebugConfig::from_env();
			for (unsigned int i = 0; i < G.n_ind; i++) {
				G.vecG[i]->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
			}
			
			// Update anchor variant encoding to reflect supersite genotype status
			// This must run after context is set but before segment building
			updateSuperSiteAnchorEncoding(G, super_sites, super_site_var_index);

			for (unsigned int i = 0; i < G.n_ind; i++) {
				G.vecG[i]->snapshotSupersiteClasses(super_sites, super_site_var_index);
				G.vecG[i]->snapshotSupersiteBaseClasses(super_sites, super_site_var_index);

				// Lightweight invariant check: ensure initial c0/c1 snapshots are compatible
				if (ss_cfg.guards_enabled) {
					supersite_invariants::SupersiteInvariantViolation viol;
					if (!supersite_invariants::check_supersite_consistency_for_sample(
							*G.vecG[i],
							super_sites,
							super_site_var_index,
							ss_cfg,
							&viol)) {
						if (ss_cfg.verbose) {
							std::fprintf(stderr,
								"[supersite-invariant] init sample=%s ss_idx=%u bp=%u: %s\n",
								G.vecG[i]->name.c_str(),
								static_cast<unsigned>(viol.ss_idx),
								static_cast<unsigned>(super_sites[viol.ss_idx].bp),
								viol.message.c_str());
						}
					}
				}
			}
		}

	//step9: Initialize genotype structures
	genotype_builder(G, options["thread"].as < int > ()).build();

	//step10: Allocate data structures for computations
	unsigned int max_number_transitions = G.largestNumberOfTransitions();
	unsigned int max_number_missing = G.largestNumberOfMissings();
	threadData = vector < compute_job >(options["thread"].as < int > (), 
	                                     compute_job(V, G, H, max_number_transitions, max_number_missing,
	                                                enable_supersites ? &super_sites : nullptr,
	                                                enable_supersites ? &locus_to_super_idx : nullptr,
	                                                enable_supersites ? &super_site_var_index : nullptr));
}

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
#include <objects/super_site_builder.h>

using namespace std;

void * phaseWindow_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_ind) vrb.progress("  * HMM computations", id_job*1.0/S->G.n_ind);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_ind) S->phaseWindow(id_worker, id_job);
		else pthread_exit(NULL);
	}
}

void phaser::phaseWindow(int id_worker, int id_job) {
	threadData[id_worker].make(id_job, options["hmm-window"].as < double > ());

	//HMM compute in windows
	for (int w = 0 ; w < threadData[id_worker].size() ; w ++) {
		const size_t kstates_size = threadData[id_worker].Kstates[w].size();
		if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
		statH.push(kstates_size * 1.0);
		statS.push(threadData[id_worker].Windows.W[w].lengthBP(V) * 1.0e-6);
		if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);

		int outcome = 0;

		if (G.vecG[id_job]->double_precision) {
			//Run using double precision as underflow happened previously
			haplotype_segment_double HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M,
				enable_supersites ? &super_sites : nullptr,
				enable_supersites ? &is_super_site : nullptr,
				enable_supersites ? &locus_to_super_idx : nullptr,
				enable_supersites ? packed_allele_codes.data() : nullptr,
				enable_supersites ? packed_allele_codes.size() : 0,
				enable_supersites ? &super_site_var_index : nullptr,
				nullptr);
			HS.forward();
			outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M,
			                      enable_supersites ? &threadData[id_worker].SC : nullptr,
			                      enable_supersites ? &threadData[id_worker].anchor_has_missing : nullptr,
			                      enable_supersites ? &threadData[id_worker].supersite_sc_offset : nullptr);
		} else {
			//Try single precision as this is faster
			haplotype_segment_single HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M,
				enable_supersites ? &super_sites : nullptr,
				enable_supersites ? &is_super_site : nullptr,
				enable_supersites ? &locus_to_super_idx : nullptr,
				enable_supersites ? packed_allele_codes.data() : nullptr,
				enable_supersites ? packed_allele_codes.size() : 0,
				enable_supersites ? &super_site_var_index : nullptr);
			HS.forward();
			outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M,
			                      enable_supersites ? &threadData[id_worker].SC : nullptr,
			                      enable_supersites ? &threadData[id_worker].anchor_has_missing : nullptr,
			                      enable_supersites ? &threadData[id_worker].supersite_sc_offset : nullptr);
			const auto* shared_ss_panel_matrix = enable_supersites ? HS.get_ss_panel_matrix() : nullptr;

			//Underflow happening with single precision, rerun using double precision
			if (outcome != 0) {
				haplotype_segment_double HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M,
					enable_supersites ? &super_sites : nullptr,
					enable_supersites ? &is_super_site : nullptr,
					enable_supersites ? &locus_to_super_idx : nullptr,
					enable_supersites ? packed_allele_codes.data() : nullptr,
					enable_supersites ? packed_allele_codes.size() : 0,
					enable_supersites ? &super_site_var_index : nullptr,
					shared_ss_panel_matrix);
				HS.forward();
				outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M,
				                      enable_supersites ? &threadData[id_worker].SC : nullptr,
					              enable_supersites ? &threadData[id_worker].anchor_has_missing : nullptr,
					              enable_supersites ? &threadData[id_worker].supersite_sc_offset : nullptr);
				G.vecG[id_job]->double_precision = true;
				n_underflow_recovered_precision++;
			}
		}

		//
		switch (outcome) {
		case -2: vrb.error("Diploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		case -1: vrb.error("Haploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		}
		n_underflow_recovered_summing += outcome;
	}

	//Copy over new IBD2 constraints into H
	if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
	H.Kbanned.pushIBD2(id_job, threadData[id_worker].Kbanned);
	if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);

	// Phase 3: Set supersite context for multivariant imputation
	if (enable_supersites) {
		G.vecG[id_job]->setSuperSiteContext(
			&super_sites,
			&locus_to_super_idx,
			&super_site_var_index,
			&threadData[id_worker].SC,
			&threadData[id_worker].anchor_has_missing,
			&threadData[id_worker].supersite_sc_offset
		);
	}

	//Sampling / Merging / Storing
	vector < bool > flagMerges;
	switch (iteration_types[iteration_stage]) {
	case STAGE_BURN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						break;
	case STAGE_PRUN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						G.vecG[id_job]->mapMerges(threadData[id_worker].T, options["mcmc-prune"].as < double > (), flagMerges);
						G.vecG[id_job]->performMerges(threadData[id_worker].T, flagMerges);
						break;
	case STAGE_MAIN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						G.vecG[id_job]->store(threadData[id_worker].T, threadData[id_worker].M);
						break;
	}
}

void phaser::phaseWindow() {
	tac.clock();
	int n_thread = options["thread"].as < int > ();
	n_underflow_recovered_summing = 0;
	n_underflow_recovered_precision = 0;
	i_workers = 0; i_jobs = 0;
	statH.clear(); statS.clear();
	storedKsizes.clear();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, phaseWindow_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_ind ; i ++) {
		phaseWindow(0, i);
		vrb.progress("  * HMM computations", (i+1)*1.0/G.n_ind);
	}
	vrb.bullet("HMM computations [K=" + stb.str(statH.mean(), 1) + "+/-" + stb.str(statH.sd(), 1) + " / W=" + stb.str(statS.mean(), 2) + "Mb / US=" + stb.str(n_underflow_recovered_summing) + " / UP=" + stb.str(n_underflow_recovered_precision) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

void phaser::phase() {
	unsigned long n_old_segments = G.numberOfSegments(), n_new_segments = 0;
	for (iteration_stage = 0 ; iteration_stage < iteration_counts.size() ; iteration_stage ++) {
		for (int iter = 0 ; iter < iteration_counts[iteration_stage] ; iter ++) {
			//VERBOSE
			std::string stage_label = "unknown";
			switch (iteration_types[iteration_stage]) {
			case STAGE_BURN:	stage_label = "burn"; vrb.title("Burn-in iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_PRUN:	stage_label = "prune"; vrb.title("Pruning iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_MAIN:	stage_label = "main"; vrb.title("Main iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			}
			const std::string iteration_label = stage_label + "-iter" + stb.str(iter+1) + "of" + stb.str(iteration_counts[iteration_stage]);
            //SELECT NEW STATES WITH PBWT
            H.select();
			
			//PHASE DATA
			phaseWindow();
			//MERGE IBD2 PAIRS
			H.Kbanned.collapse();
            //UPDATE H with new sampled haplotypes
            H.updateHaplotypes(G);
            //TRANSPOSE H from Hfirst to Vfirst (for next PBWT compute)
            H.transposeHaplotypes_H2V(false);
            // Rebuild supersite metadata AFTER transpose so H.H_opt_var reflects new panel state
            if (enable_supersites) {
                rebuildSupersiteMetadata(iteration_label + "/post-update");
            }
            //if (options.count("pedigree")) H.checkScaffoldPedigrees(G, options["pedigree"].as < string > ());
			//UPDATE PS after prunning
			if (iteration_types[iteration_stage] == STAGE_PRUN) {
				n_new_segments = G.numberOfSegments();
				vrb.bullet("Trimming [pc=" + stb.str((1-n_new_segments*1.0/n_old_segments)*100, 2) + "%]");
			}
		}
	}
}

void phaser::rebuildSupersiteMetadata(const std::string& context) {
	if (!enable_supersites) {
		H.clearSupersitePBWTContext();
		return;
	}

	super_sites.clear();
	is_super_site.clear();
	packed_allele_codes.clear();
	locus_to_super_idx.clear();
	super_site_var_index.clear();

	buildSuperSites(V, H,
	                super_sites,
	                is_super_site,
	                packed_allele_codes,
	                locus_to_super_idx,
	                super_site_var_index);
	M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

	supersite_build_last_context = context;
	++supersite_build_counter;
	vrb.bullet("Supersite metadata build #" + stb.str(supersite_build_counter, 0) +
	           " [" + context + "] (n_sites=" + stb.str(super_sites.size(), 0) +
	           ", packed_bytes=" + stb.str(packed_allele_codes.size(), 0) + ")");

	applySupersiteAnchorGuards();
	if (supersite_pbwt_enabled) {
		H.setSupersitePBWTContext(&super_sites, &locus_to_super_idx, &packed_allele_codes);
	}

	for (unsigned int i = 0; i < G.n_ind; ++i) {
		G.vecG[i]->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
		G.vecG[i]->snapshotSupersitePhasedGts(super_sites, super_site_var_index);
	}
}

void phaser::applySupersiteAnchorGuards() {
	if (!enable_supersites || super_sites.empty()) return;
	H.applySupersiteAnchorMask(super_sites, super_site_var_index);
	supersite_anchor_redirect = buildSupersiteAnchorMap(super_sites, super_site_var_index, V.size());
	H.setSupersiteAnchorRedirect(supersite_anchor_redirect);
}

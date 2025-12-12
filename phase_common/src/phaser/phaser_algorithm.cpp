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
#include <algorithm>
#include <sstream>
#include <iostream>
#include <unordered_set>

using namespace std;

namespace {

bool supersite_metadata_trace_enabled() {
	static int flag = -1;
	if (flag < 0) {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
		flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
	}
	return flag == 1;
}

// Parse comma-separated bp targets from SHAPEIT5_TRACE_BP
const std::unordered_set<int>& bp_trace_targets() {
	static std::unordered_set<int> targets;
	static bool initialized = false;
	if (!initialized) {
		initialized = true;
		const char* env = std::getenv("SHAPEIT5_TRACE_BP");
		if (env && env[0]) {
			std::stringstream ss(env);
			std::string token;
			while (std::getline(ss, token, ',')) {
				try {
					int bp = std::stoi(token);
					targets.insert(bp);
				} catch (const std::exception&) {
					// ignore parse errors
				}
			}
		}
	}
	return targets;
}

bool should_trace_bp(int bp) {
	const auto& targets = bp_trace_targets();
	return !targets.empty() && targets.find(bp) != targets.end();
}

const SuperSite* find_superdebug_supersite(const std::vector<SuperSite>& super_sites) {
	if (debug::SUPERDEBUG_BP <= 0) return nullptr;
	for (const auto& ss : super_sites) {
		if (static_cast<int>(ss.global_site_id) == debug::SUPERDEBUG_BP) {
			return &ss;
		}
	}
	return nullptr;
}

int find_superdebug_sample_index(const genotype_set& G) {
	if (debug::SUPERDEBUG_SAMPLENAME.empty()) return -1;
	for (int i = 0; i < G.n_ind; ++i) {
		if (G.vecG[i] && G.vecG[i]->name == debug::SUPERDEBUG_SAMPLENAME) {
			return i;
		}
	}
	return -1;
}

void log_superdebug_genotype_state(const std::string& context,
                                   const genotype_set& G,
                                   const std::vector<SuperSite>& super_sites,
                                   const std::vector<int>& super_site_var_index) {
	if (debug::SUPERDEBUG_SAMPLENAME.empty()) return;
	const SuperSite* target = find_superdebug_supersite(super_sites);
	if (!target) return;
	for (int i = 0; i < G.n_ind; ++i) {
		const genotype* sample = G.vecG[i];
		if (sample && sample->name == debug::SUPERDEBUG_SAMPLENAME) {
			debug::print_supersite_state(sample, *target, super_site_var_index, context);
			break;
		}
	}
}

void log_superdebug_panel_state(const std::string& context,
                                const haplotype_set& H,
                                const genotype_set& G,
                                const std::vector<SuperSite>& super_sites,
                                const std::vector<int>& super_site_var_index) {
	if (debug::SUPERDEBUG_SAMPLENAME.empty()) return;
	const SuperSite* target = find_superdebug_supersite(super_sites);
	if (!target) return;
	const int sample_index = find_superdebug_sample_index(G);
	if (sample_index < 0) return;
	const int hap0 = 2 * sample_index;
	const int hap1 = hap0 + 1;
	if (hap0 >= static_cast<int>(H.n_hap) || hap1 >= static_cast<int>(H.n_hap)) return;

	std::ostringstream hap0_line;
	std::ostringstream hap1_line;
	hap0_line << "[SUPERDEBUG] Sample=" << debug::SUPERDEBUG_SAMPLENAME
	          << " Pos=" << target->global_site_id
	          << " Context='" << context << "' Panel HAP0:";
	hap1_line << "  Panel HAP1:";

	for (uint16_t ai = 0; ai < target->var_count; ++ai) {
		const size_t offset = target->var_start + ai;
		if (offset >= super_site_var_index.size()) break;
		const int locus = super_site_var_index[offset];
		int bit0 = -1;
		int bit1 = -1;
		if (locus >= 0 && locus < static_cast<int>(H.H_opt_var.n_rows)) {
			if (hap0 >= 0 && hap0 < static_cast<int>(H.H_opt_var.n_cols)) {
				bit0 = H.H_opt_var.get(static_cast<unsigned int>(locus), static_cast<unsigned int>(hap0));
			}
			if (hap1 >= 0 && hap1 < static_cast<int>(H.H_opt_var.n_cols)) {
				bit1 = H.H_opt_var.get(static_cast<unsigned int>(locus), static_cast<unsigned int>(hap1));
			}
		}
		hap0_line << " " << bit0;
		hap1_line << " " << bit1;
	}

	std::cout << hap0_line.str() << std::endl;
	std::cout << hap1_line.str() << std::endl;
}

void log_superdebug_packed_codes(const std::string& context,
                                 const genotype_set& G,
                                 const std::vector<SuperSite>& super_sites,
                                 const std::vector<uint8_t>& packed_codes) {
	if (debug::SUPERDEBUG_SAMPLENAME.empty()) return;
	if (packed_codes.empty()) return;
	const SuperSite* target = find_superdebug_supersite(super_sites);
	if (!target) return;
	const int sample_index = find_superdebug_sample_index(G);
	if (sample_index < 0) return;
	const int hap0 = 2 * sample_index;
	const int hap1 = hap0 + 1;
	uint8_t code0 = unpackSuperSiteCode(packed_codes.data(), target->panel_offset, static_cast<uint32_t>(hap0));
	uint8_t code1 = unpackSuperSiteCode(packed_codes.data(), target->panel_offset, static_cast<uint32_t>(hap1));

	std::cout << "[SUPERDEBUG] Sample=" << debug::SUPERDEBUG_SAMPLENAME
	          << " Pos=" << target->global_site_id
	          << " Context='" << context << "' PackedCodes h0="
	          << static_cast<unsigned>(code0)
	          << " h1=" << static_cast<unsigned>(code1)
	          << std::endl;
}

} // namespace

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
		if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
		statH.push(threadData[id_worker].Kstates[w].size()*1.0);
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
				enable_supersites ? &super_site_var_index : nullptr);
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

			//Underflow happening with single precision, rerun using double precision
			if (outcome != 0) {
				haplotype_segment_double HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M,
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
				G.vecG[id_job]->double_precision = true;
				n_underflow_recovered_precision++;
			}
		}

		//
		switch (outcome) {
		case -2: vrb.error("Diploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		case -1: vrb.error("Haploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		}
		if (enable_supersites && threadData[id_worker].sc_buffer_active() && !threadData[id_worker].verify_sc_guards()) {
			vrb.error("Supersite SC guard corrupted for [" + G.vecG[id_job]->name + "]");
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
	std::vector<uint8_t> supersite_codes_before_update;
	static bool bp_trace_emitted = false;
	if (!bp_trace_emitted && !bp_trace_targets().empty()) {
		bp_trace_emitted = true;
		for (int locus = 0; locus < V.size(); ++locus) {
			const variant* vp = V.vec_pos[locus];
			if (!vp) continue;
			if (should_trace_bp(vp->bp)) {
				std::ostringstream oss;
				oss << "[BP_INDEX_TRACE] locus=" << locus
				    << " chr=" << vp->chr
				    << " bp=" << vp->bp
				    << " id=" << vp->id
				    << " ref=" << vp->ref
				    << " alt=" << vp->alt;
				vrb.bullet(oss.str());
			}
		}
	}
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
			supersite_codes_before_update.clear();
            //SELECT NEW STATES WITH PBWT
            H.select();
            // Snapshot current packed codes before haplotype refresh
            if (enable_supersites) {
                supersite_codes_before_update = packed_allele_codes;
            }
			
			//PHASE DATA
			phaseWindow();
			//MERGE IBD2 PAIRS
			H.Kbanned.collapse();
			if (enable_supersites && !debug::SUPERDEBUG_SAMPLENAME.empty()) {
				log_superdebug_genotype_state(iteration_label + "/pre-hap-update sample", G, super_sites, super_site_var_index);
				log_superdebug_panel_state(iteration_label + "/pre-hap-update panel", H, G, super_sites, super_site_var_index);
			}
            //UPDATE H with new sampled haplotypes
            H.updateHaplotypes(G);
            if (enable_supersites && !debug::SUPERDEBUG_SAMPLENAME.empty()) {
				log_superdebug_panel_state(iteration_label + "/post-hap-update panel", H, G, super_sites, super_site_var_index);
			}
            //TRANSPOSE H from Hfirst to Vfirst (for next PBWT compute)
            H.transposeHaplotypes_H2V(false);
            // Rebuild supersite metadata AFTER transpose so H.H_opt_var reflects new panel state
            if (enable_supersites) {
                rebuildSupersiteMetadata(iteration_label + "/post-update", &supersite_codes_before_update);
                if (!debug::SUPERDEBUG_SAMPLENAME.empty()) {
					log_superdebug_panel_state(iteration_label + "/post-rebuild panel", H, G, super_sites, super_site_var_index);
					log_superdebug_packed_codes(iteration_label + "/post-rebuild packed", G, super_sites, packed_allele_codes);
				}
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

void phaser::logPackedCodeDiff(const std::string& context, const std::vector<uint8_t>& before, const std::vector<uint8_t>& after) const {
	if (!enable_supersites) return;
	const size_t min_size = std::min(before.size(), after.size());
	size_t changed = 0;
	for (size_t i = 0; i < min_size; ++i) {
		if (before[i] != after[i]) ++changed;
	}
	const size_t added = (after.size() > before.size()) ? (after.size() - before.size()) : 0;
	const size_t removed = (before.size() > after.size()) ? (before.size() - after.size()) : 0;
	vrb.bullet("Supersite packed-code diff [" + context + "] changed=" + stb.str(changed, 0) +
	           " added=" + stb.str(added, 0) + " removed=" + stb.str(removed, 0));
	if (supersite_metadata_trace_enabled() && changed > 0) {
		std::ostringstream oss;
		oss << "  first-changed-bytes";
		size_t reported = 0;
		for (size_t i = 0; i < min_size && reported < 8; ++i) {
			if (before[i] != after[i]) {
				oss << " [" << i << ":" << static_cast<int>(before[i]) << "->" << static_cast<int>(after[i]) << "]";
				++reported;
			}
		}
		vrb.bullet(oss.str());
	}
}

void phaser::traceSupersiteAnchors(const std::string& context, const std::vector<uint8_t>& codes_snapshot, size_t max_sites, size_t max_haps) const {
	if (!enable_supersites || !supersite_metadata_trace_enabled()) return;
	if (super_sites.empty() || codes_snapshot.empty()) return;
	const size_t n_sites = std::min(max_sites, super_sites.size());
	const size_t n_haps = std::min(max_haps, static_cast<size_t>(H.n_hap));
	for (size_t ss_idx = 0; ss_idx < n_sites; ++ss_idx) {
		const SuperSite& ss = super_sites[ss_idx];
		std::ostringstream oss;
		oss << "[SupersiteTrace] " << context << " ss_idx=" << ss_idx << " locus=" << ss.global_site_id;
		for (size_t hap = 0; hap < n_haps; ++hap) {
			uint8_t code = unpackSuperSiteCode(codes_snapshot.data(), ss.panel_offset, hap);
			int panel_bit = -1;
			if (ss.global_site_id < H.H_opt_var.n_rows && hap < H.H_opt_var.n_cols) {
				panel_bit = H.H_opt_var.get(ss.global_site_id, hap);
			}
			oss << " h" << hap << "(code=" << static_cast<int>(code) << ",panel=" << panel_bit << ")";
		}
		vrb.bullet(oss.str());
	}
}

void phaser::rebuildSupersiteMetadata(const std::string& context, const std::vector<uint8_t>* diff_against) {
	if (!enable_supersites) return;

	super_sites.clear();
	is_super_site.clear();
	packed_allele_codes.clear();
	locus_to_super_idx.clear();
	super_site_var_index.clear();
	sample_supersite_genotypes.clear();

	buildSuperSites(V, H,
	                super_sites,
	                is_super_site,
	                packed_allele_codes,
	                locus_to_super_idx,
	                super_site_var_index,
	                sample_supersite_genotypes,
	                supersite_mac_threshold);
	M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

	supersite_build_last_context = context;
	++supersite_build_counter;
	vrb.bullet("Supersite metadata build #" + stb.str(supersite_build_counter, 0) +
	           " [" + context + "] (n_sites=" + stb.str(super_sites.size(), 0) +
	           ", packed_bytes=" + stb.str(packed_allele_codes.size(), 0) + ")");

	if (diff_against) {
		logPackedCodeDiff(context, *diff_against, packed_allele_codes);
	}

	traceSupersiteAnchors(context, packed_allele_codes);

	applySupersiteAnchorGuards();

	for (unsigned int i = 0; i < G.n_ind; ++i) {
		G.vecG[i]->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
		G.vecG[i]->snapshotSupersiteClasses(super_sites, super_site_var_index);
	}
}

void phaser::applySupersiteAnchorGuards() {
	if (!enable_supersites || super_sites.empty()) return;
	H.applySupersiteAnchorMask(super_sites, super_site_var_index);
	supersite_anchor_redirect = buildSupersiteAnchorMap(super_sites, super_site_var_index, V.size());
	H.setSupersiteAnchorRedirect(supersite_anchor_redirect);
}

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

#include <objects/compute_job.h>
#include <models/supersite_trace_utils.h>

#include <cmath>
#include <cstdio>

using namespace std;

#define MAX_OVERLAP_HETS 0.75f
#define N_RANDOM_HAPS 100

compute_job::compute_job(variant_map & _V, genotype_set & _G, conditioning_set & _H, unsigned int n_max_transitions, unsigned int n_max_missing,
												 const std::vector<SuperSite>* ss,
												 const std::vector<int>* locus_ss_idx,
												 const std::vector<int>* ss_var_idx) 
		: V(_V), G(_G), H(_H), super_sites(ss), locus_to_super_idx(locus_ss_idx), super_site_var_index(ss_var_idx),
			sc_guard_value(std::numeric_limits<float>::quiet_NaN()), sc_guard_active(false) {
	T = vector < double > (n_max_transitions, 0.0);
	M = vector < float > (n_max_missing , 0.0);
	Ordering = vector < unsigned int > (H.n_hap);
	iota(Ordering.begin(), Ordering.end(), 0);
	Oiterator = 0;
}

compute_job::~compute_job() {
	free();
}

void compute_job::free () {
	vector < double > ().swap(T);
	vector < float > ().swap(M);
	vector < vector < unsigned int > > ().swap(Kstates);
	Kbanned.clear();
	Windows.clear();
	SC.clear();
	sc_guard_active = false;
}

void compute_job::make(unsigned int ind, double min_window_size) {
	//1. Mapping coordinates of each segment
	int n_windows = Windows.build (V, G.vecG[ind], min_window_size);

	//2. Update conditional haps
	unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
	Kstates = vector < vector < unsigned int > > (n_windows, vector < unsigned int >());
	unsigned long curr_hap0 = 2*ind+0, curr_hap1 = 2*ind+1;
	for (int w = 0 ; w < n_windows ; w++) {
		vector < int > phap = vector < int > (2 * H.depth, -1);
		for (int l = Windows.W[w].start_locus ; l <= Windows.W[w].stop_locus ; l++) {
			if (H.sites_pbwt_selection[l]) {
				for (int s = 0 ; s < H.depth ; s ++) {
						// Diagnostic guard: compute the flat indices and check bounds before accessing
						size_t vec_size = H.indexes_pbwt_neighbour.size();
						unsigned long idx0 = static_cast<unsigned long>(s) * addr_offset + curr_hap0 * static_cast<unsigned long>(H.sites_pbwt_ngroups) + static_cast<unsigned long>(H.sites_pbwt_grouping[l]);
						unsigned long idx1 = static_cast<unsigned long>(s) * addr_offset + curr_hap1 * static_cast<unsigned long>(H.sites_pbwt_ngroups) + static_cast<unsigned long>(H.sites_pbwt_grouping[l]);
						if (idx0 >= vec_size || idx1 >= vec_size) {
							std::fprintf(stderr, "PBWT index OOB: idx0=%lu idx1=%lu vec_size=%zu addr_offset=%lu s=%d curr_hap0=%lu curr_hap1=%lu grouping=%d l=%d H.sites_pbwt_ngroups=%d H.n_ind=%d H.n_hap=%d H.depth=%d n_windows=%d\n",
								idx0, idx1, vec_size, addr_offset, s, curr_hap0, curr_hap1, H.sites_pbwt_grouping[l], l, H.sites_pbwt_ngroups, H.n_ind, H.n_hap, H.depth, n_windows);
							std::abort();
						}
						int cond_hap0 = H.indexes_pbwt_neighbour[idx0];
						int cond_hap1 = H.indexes_pbwt_neighbour[idx1];
					if ((cond_hap0 >= 0) && (cond_hap0 != phap[2*s+0])) { Kstates[w].push_back(cond_hap0); phap[2*s+0] = cond_hap0; };
					if ((cond_hap1 >= 0) && (cond_hap1 != phap[2*s+1])) { Kstates[w].push_back(cond_hap1); phap[2*s+1] = cond_hap1; };
				}
			}
		}
		sort(Kstates[w].begin(), Kstates[w].end());
		Kstates[w].erase(unique(Kstates[w].begin(), Kstates[w].end()), Kstates[w].end());
	}

	//3. Protect for IBD2
	Kbanned.clear();
	for (int w = 0 ; w < n_windows; w++) {
		vector < int > toBeRemoved;

		//3.1. Identify potential IBD2 pairs
		for (int k = 1; k < Kstates[w].size() ; k++) {
			unsigned int ind0 = Kstates[w][k-1]/2;
			unsigned int ind1 = Kstates[w][k]/2;
			if (ind0 == ind1 && ind0 < G.n_ind && !G.vecG[ind0]->haploid) {
				float het_overlap = H.H_opt_hap.getMatchHets(ind, ind0, Windows.W[w].start_locus, Windows.W[w].stop_locus);
				if (het_overlap > 0.75) {
					toBeRemoved.push_back(k-1);
					toBeRemoved.push_back(k);
					Kbanned.emplace_back(ind0, Windows.W[w].start_locus, Windows.W[w].stop_locus);
					//cout << "IBD2 : " << G.vecG[ind]->name << " vs " << G.vecG[ind0]->name << " / P = " << stb.str(het_overlap*100.0, 2) << endl;
				}
			}
		}

		//3.2. Remove potential IBD2 states from conditioning set
		if (toBeRemoved.size() > 0) {
			vector < unsigned int > Ktmp;
			Ktmp.reserve(Kstates[w].size() - toBeRemoved.size());
			for (int k = 0, p = 0; k < Kstates[w].size() ; k++) {
				if (p < toBeRemoved.size() && toBeRemoved[p] == k) p++;
				else Ktmp.push_back(Kstates[w][k]);
			}

			sort(Ktmp.begin(), Ktmp.end());
			Ktmp.erase(unique(Ktmp.begin(), Ktmp.end()), Ktmp.end());
			Kstates[w] = Ktmp;
		}
	}

	//4. Protect for #states = 0
	for (int w = 0 ; w < n_windows; w++) {
		if (Kstates[w].size() < 2) {
			for (int i = 0 ; i < N_RANDOM_HAPS ; i++) {
				int random_state = Ordering[Oiterator];
				if (random_state/2 != ind) Kstates[w].push_back(random_state);
				Oiterator=((Oiterator+1)==H.n_hap)?0:(Oiterator+1);
			}
			sort(Kstates[w].begin(), Kstates[w].end());
			Kstates[w].erase(unique(Kstates[w].begin(), Kstates[w].end()), Kstates[w].end());
			vrb.warning("No PBWT states found [" + G.vecG[ind]->name  + " / w=" + stb.str(w) + "] / Using " + stb.str(Kstates[w].size()) + " random states");
		}
	}
	
	//5. Phase 3: Populate anchor_has_missing and allocate SC with thread-local offsets
	if (super_sites && locus_to_super_idx && super_site_var_index) {
		anchor_has_missing.assign(super_sites->size(), false);
		supersite_sc_offset.assign(super_sites->size(), 0);
		
		// Check each supersite: if ALL members are missing for this sample, set flag
		for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
			const SuperSite& ss = (*super_sites)[ss_idx];
			bool all_missing = true;
			
			for (uint32_t i = 0; i < ss.var_count; ++i) {
				int v_idx = (*super_site_var_index)[ss.var_start + i];
				unsigned char v = G.vecG[ind]->Variants[DIV2(v_idx)];
				if (!VAR_GET_MIS(MOD2(v_idx), v)) {
					all_missing = false;
					break;
				}
			}
			
			anchor_has_missing[ss_idx] = all_missing;
		}
		if (supersite_trace_enabled()) {
			size_t flagged = std::count(anchor_has_missing.begin(), anchor_has_missing.end(), true);
			supersite_trace_log("[SupersiteSC] sample=%s anchors_missing=%zu/%zu\n",
			                    G.vecG[ind]->name.c_str(),
			                    flagged,
			                    anchor_has_missing.size());
		}
		
		// Allocate SC: compute total size and set thread-local offsets for each supersite
		uint32_t total_size = 0;
		for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
			const SuperSite& ss = (*super_sites)[ss_idx];
			// Note: ss.n_classes is now set immutably during buildSuperSites()
			
			if (anchor_has_missing[ss_idx]) {
				supersite_sc_offset[ss_idx] = total_size;
				total_size += HAP_NUMBER * ss.n_classes;  // 8 lanes × C classes
			} else {
				supersite_sc_offset[ss_idx] = 0;  // Not used, but set for consistency
			}
		}
		
		if (total_size > 0) {
			SC.assign(total_size + 2, 0.0f);
			sc_guard_active = true;
			const float guard = sc_guard_value;
			SC.front() = guard;
			SC.back() = guard;
			for (size_t ss_idx = 0; ss_idx < supersite_sc_offset.size(); ++ss_idx) {
				if (anchor_has_missing[ss_idx]) supersite_sc_offset[ss_idx] += 1;
			}
			if (supersite_trace_enabled()) {
				supersite_trace_log("[SupersiteSC] sample=%s allocated=%zu floats\n",
				                    G.vecG[ind]->name.c_str(),
				                    SC.size());
				for (size_t ss_idx = 0, reported = 0; ss_idx < supersite_sc_offset.size() && reported < 4; ++ss_idx) {
					if (anchor_has_missing[ss_idx]) {
						supersite_trace_log("  offset ss_idx=%zu -> %u (classes=%u)\n",
						                    ss_idx,
						                    supersite_sc_offset[ss_idx],
						                    static_cast<unsigned>((*super_sites)[ss_idx].n_classes));
						++reported;
					}
				}
			}
		} else {
			SC.clear();
			sc_guard_active = false;
			if (supersite_trace_enabled()) {
				supersite_trace_log("[SupersiteSC] sample=%s no missing supersite anchors (SC cleared)\n",
				                    G.vecG[ind]->name.c_str());
			}
		}
	} else {
		SC.clear();
		sc_guard_active = false;
	}
}

bool compute_job::verify_sc_guards() const {
	if (!sc_guard_active || SC.size() < 2) return true;
	return std::isnan(SC.front()) && std::isnan(SC.back());
}

bool compute_job::sc_buffer_active() const {
	return sc_guard_active && SC.size() >= 2;
}

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

#include <containers/conditioning_set/conditioning_set_header.h>
#include <models/super_site_accessor.h>
#include <algorithm>

using namespace std;

bool conditioning_set::isSupersiteAnchor(int locus) const {
	if (locus < 0) return false;
	if (supersite_pbwt_super_sites == nullptr || supersite_pbwt_locus_to_super_idx == nullptr) return false;
	if (static_cast<size_t>(locus) >= supersite_pbwt_locus_to_super_idx->size()) return false;
	const int ss_idx = (*supersite_pbwt_locus_to_super_idx)[locus];
	if (ss_idx < 0 || static_cast<size_t>(ss_idx) >= supersite_pbwt_super_sites->size()) return false;
	return (*supersite_pbwt_super_sites)[ss_idx].global_site_id == static_cast<uint32_t>(locus);
}

uint8_t conditioning_set::class_code(int locus, int hap) const {
	if (hap < 0 || hap >= n_hap) return 0;
	if (!isSupersiteAnchor(locus)) {
		return static_cast<uint8_t>(H_opt_var.get(locus, hap));
	}

	// Supersite anchor path
	if (supersite_pbwt_super_sites == nullptr || supersite_pbwt_packed_codes == nullptr) return 0;
	const int ss_idx = (*supersite_pbwt_locus_to_super_idx)[locus];
	if (ss_idx < 0 || static_cast<size_t>(ss_idx) >= supersite_pbwt_super_sites->size()) return 0;
	const SuperSite& ss = (*supersite_pbwt_super_sites)[ss_idx];
	const uint8_t code = unpackSuperSiteCode(supersite_pbwt_packed_codes, ss.panel_offset, static_cast<uint32_t>(hap));
	if (code >= ss.n_classes) return 0;
	return code;
}

void * selecter_callback(void * ptr) {
	conditioning_set * S = static_cast< conditioning_set * >( ptr );

	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_worker ++;
	pthread_mutex_unlock(&S->mutex_workers);

	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_job ++;
		pthread_mutex_unlock(&S->mutex_workers);

		if (id_job <= S->sites_pbwt_mthreading.back()) {
			S->select(id_job);
			pthread_mutex_lock(&S->mutex_workers);
			vrb.progress("  * PBWT selection", (++S->d_job)*1.0/(S->sites_pbwt_mthreading.back()+1));
			pthread_mutex_unlock(&S->mutex_workers);
		}
		else pthread_exit(NULL);
	}
}

void conditioning_set::transposePBWTneighbours() {
	int block = 32;
	unsigned long addr_tar, addr_src;
	unsigned long addr_offset = sites_pbwt_ngroups * n_ind * 2UL;
	for (int d = 0; d < depth ; d ++) {
		for (int s = 0; s < sites_pbwt_ngroups ; s += block) {
			for(int h = 0; h < n_ind * 2; ++h) {
				for(int b = 0; b < block && s + b < sites_pbwt_ngroups ; ++b) {
					addr_tar = depth * addr_offset + h*sites_pbwt_ngroups + s + b;
					addr_src = d * addr_offset + (s + b)*n_ind*2UL + h;
					indexes_pbwt_neighbour[addr_tar] = indexes_pbwt_neighbour[addr_src];
				}
			}
		}
		std::copy(indexes_pbwt_neighbour.begin() + depth * addr_offset , indexes_pbwt_neighbour.end(), indexes_pbwt_neighbour.begin() + d * addr_offset );
	}
}


void conditioning_set::select(int chunk) {
	vector < int > A = vector < int > (n_hap, 0);
	vector < int > B = vector < int > (n_hap, 0);
	vector < int > C = vector < int > (n_hap, 0);
	vector < int > D = vector < int > (n_hap, 0);
	vector < uint8_t > bucket_id;
	vector < int > bucket_p;
	iota(A.begin(), A.end(), 0);
	fill(C.begin(), C.end(), 0);

	for (int l = 0 ; l < n_site ; l ++) {
		bool eval = sites_pbwt_evaluation[l];
		bool selc = sites_pbwt_selection[l];
		bool chnk = (sites_pbwt_mthreading[l] == chunk);
		bool buff = (sites_pbwt_mthreading[l] < chunk) && (l >= starts_pbwt_mthreading[chunk]);

		if (eval && (chnk || buff)) {
			if (supersite_pbwt_active() && isSupersiteAnchor(l)) {
				const int ss_idx = (*supersite_pbwt_locus_to_super_idx)[l];
				const SuperSite& ss = (*supersite_pbwt_super_sites)[ss_idx];
				const int t = static_cast<int>(ss.n_classes);
				if (static_cast<int>(bucket_id.size()) < n_hap) bucket_id.resize(n_hap);
				if (static_cast<int>(bucket_p.size()) < t) bucket_p.resize(t);
				// Start each bucket divergence at current locus (matches binary p/q init).
				std::fill(bucket_p.begin(), bucket_p.begin() + t, l);

				// Generalized Durbin update: update all p[c] per scanned hap, emit p[class], reset that bucket.
				int u = 0;
				for (int h = 0; h < n_hap; h++) {
					const uint8_t cls = class_code(l, A[h]);
					const int dlookup = C[h];
					for (int c = 0; c < t; ++c) {
						if (dlookup > bucket_p[c]) bucket_p[c] = dlookup;
					}
					bucket_id[u] = cls;
					B[u] = A[h];
					D[u] = bucket_p[cls];
					bucket_p[cls] = 0;
					++u;
				}

				// Stable bucket emit in class-id order 0..t-1
				int out = 0;
				for (int cls = 0; cls < t; ++cls) {
					for (int h = 0; h < n_hap; ++h) {
						if (bucket_id[h] == cls) {
							A[out] = B[h];
							C[out] = D[h];
							++out;
						}
					}
				}
			} else {
				int u = 0, v = 0, p = l, q = l;
				for (int h = 0 ; h < n_hap ; h ++) {
					int alookup = A[h], dlookup = C[h];
					if (dlookup > p) p = dlookup;
					if (dlookup > q) q = dlookup;
					if (!H_opt_var.get(l, alookup)) {
						A[u] = alookup;
						C[u] = p;
						p = 0;
						u++;
					} else {
						B[v] = alookup;
						D[v] = q;
						q = 0;
						v++;
					}
				}
				std::copy(B.begin(), B.begin()+v, A.begin()+u);
				std::copy(D.begin(), D.begin()+v, C.begin()+u);
			}
			if (selc && chnk) store(l, A, C);
		}
	}
}

void conditioning_set::store(int l, vector < int > & A, vector < int > & C) {
	unsigned long addr_offset = sites_pbwt_ngroups * n_ind * 2UL;
	for (int h = 0 ; h < n_hap ; h ++) {
		int chap = A[h];
		int cind = chap / 2;
		if (cind < n_ind) {
			int add_guess0 = 0, add_guess1 = 0, offset0 = 1, offset1 = 1, hap_guess0 = -1, hap_guess1 = -1, div_guess0 = -1, div_guess1 = -1;
			unsigned long tar_idx = sites_pbwt_grouping[l] * 2UL * n_ind + chap;
			for (int n_added = 0 ; n_added < depth ; ) {
				if ((h-offset0)>=0) {
					hap_guess0 = A[h-offset0];
					div_guess0 = max(C[h-offset0+1], div_guess0);
					add_guess0 = Kbanned.noIBD2(chap, hap_guess0, l);
				} else { add_guess0 = 0; div_guess0 = l+1; }
				if ((h+offset1)<n_hap) {
					hap_guess1 = A[h+offset1];
					div_guess1 = max(C[h+offset1], div_guess1);
					add_guess1 = Kbanned.noIBD2(chap, hap_guess1, l);
				} else { add_guess1 = 0; div_guess1 = l+1; }
				if (add_guess0 && add_guess1) {
					if (div_guess0 < div_guess1) {
						indexes_pbwt_neighbour[n_added*addr_offset+tar_idx] = hap_guess0;
						offset0++; n_added++;
					} else {
						indexes_pbwt_neighbour[n_added*addr_offset+tar_idx] = hap_guess1;
						offset1++; n_added++;
					}
				} else if (add_guess0) {
					indexes_pbwt_neighbour[n_added*addr_offset+tar_idx] = hap_guess0;
					offset0++; n_added++;
				} else if (add_guess1) {
					indexes_pbwt_neighbour[n_added*addr_offset+tar_idx] = hap_guess1;
					offset1++; n_added++;
				} else {
					offset0++;
					offset1++;
				}
			}
		}
	}
}

void conditioning_set::select() {
	tac.clock();
	i_worker = 0; i_job = 0, d_job = 0;

	//Select new sites at which to trigger storage
	vector < vector < int > > candidates = vector < vector < int > > (sites_pbwt_grouping.back() + 1);
	for (int l = 0 ; l < n_site ; l++) {
		if (!sites_pbwt_evaluation[l]) continue;
		int locus = l;
		if (supersite_anchor_redirect_enabled &&
		    l < static_cast<int>(supersite_anchor_redirect.size()) &&
		    supersite_anchor_redirect[l] >= 0) {
			locus = supersite_anchor_redirect[l];
		}
		candidates[sites_pbwt_grouping[l]].push_back(locus);
	}
	if (supersite_anchor_redirect_enabled) {
		for (auto& vec : candidates) {
			if (vec.size() > 1) {
				std::sort(vec.begin(), vec.end());
				vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
			}
		}
	}
	sites_pbwt_selection = vector < bool > (n_site , false);
	int n_selected_loci = 0;
	for (int g = 0 ; g < candidates.size() ; g++) {
		if (candidates[g].size() > 0) {
			sites_pbwt_selection[candidates[g][rng.getInt(candidates[g].size())]] = true;
			++n_selected_loci;
		}
	}

	//Clean up previous selected states
	fill(indexes_pbwt_neighbour.begin(), indexes_pbwt_neighbour.end() , -1);

	//Perform multi-threaded selection
	vrb.progress("  * PBWT selection", 0.0f);
	if (nthread > 1) {
		for (int t = 0 ; t < nthread ; t++) pthread_create( &id_workers[t] , NULL, selecter_callback, static_cast<void *>(this));
		for (int t = 0 ; t < nthread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int c = 0 ; c  <= sites_pbwt_mthreading.back() ; c ++) {
		select(c);
		vrb.progress("  * PBWT selection", c*1.0/(sites_pbwt_mthreading.back()+1));
	}

	//Transpose matrix with selected states
	transposePBWTneighbours();

	vrb.bullet("PBWT selection (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

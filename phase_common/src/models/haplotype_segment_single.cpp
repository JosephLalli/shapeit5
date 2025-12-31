// (Sibling helper implementations moved below after includes so required
//  types like SiteView and anonymous helpers are available.)
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

#include <models/haplotype_segment_single.h>
#include <models/site_emission_adapter.h>
#include <cstdio>
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

namespace {
// Map lane -> desired class for supersites; for biallelic, lane_class is 0/1
inline uint8_t lane_expected_class(const SiteView& sv, int lane) {
	switch (sv.kind) {
		case SiteKind::SuperAnchor:
		case SiteKind::SuperSibling:
			return sv.lane_class[lane];
		case SiteKind::Biallelic:
		default:
			return sv.lane_class[lane];
	}
}

}

// Sibling bookkeeping: update trace, maintain cursor/index sanity, no DP math
void haplotype_segment_single::handle_sibling_bookkeeping(const SiteView& site_view) {
	(void)site_view;
}

void haplotype_segment_single::INIT_SIB(const SiteView& site_view) {
	handle_sibling_bookkeeping(site_view);
	// No DP math, no prev_abs_locus advance, no AlphaMissing/AlphaSumMissing store
}

void haplotype_segment_single::RUN_SIB(const SiteView& site_view) {
	handle_sibling_bookkeeping(site_view);
	// No DP math, no prev_abs_locus advance, no AlphaMissing/AlphaSumMissing store
}

void haplotype_segment_single::COLLAPSE_SIB(const SiteView& site_view) {
	handle_sibling_bookkeeping(site_view);
	// No DP math, no prev_abs_locus advance, no AlphaMissing/AlphaSumMissing store
}



haplotype_segment_single::haplotype_segment_single(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M) :
	G(_G), M(_M),
	super_sites(_G ? _G->super_sites : nullptr),
	locus_to_super_idx(_G ? _G->locus_to_super_idx : nullptr),
	panel_codes(_G ? _G->supersite_panel_codes : nullptr),
	super_site_var_index(_G ? _G->super_site_var_index : nullptr),
	cond_idx(&idxH),
    supersite_sc_offset(nullptr),
    supersites_enabled_flag(_G && _G->super_sites && _G->locus_to_super_idx && _G->super_site_var_index && _G->supersite_panel_codes) {
	if (!supersites_enabled_flag) {
		super_sites = nullptr;
		locus_to_super_idx = nullptr;
		super_site_var_index = nullptr;
	} else if (G && G->ss_observed_gts.empty()) {
		G->snapshotSupersiteObservedGts(*super_sites, *super_site_var_index);
	}
	segment_first = W.start_segment;
	segment_last = W.stop_segment;
	locus_first = W.start_locus;
	locus_last = W.stop_locus;
	ambiguous_first = W.start_ambiguous;
	ambiguous_last = W.stop_ambiguous;
	missing_first = W.start_missing;
	missing_last = W.stop_missing;
	transition_first = W.start_transition;
	transition_last = W.stop_transition;
	n_cond_haps = idxH.size();
	n_missing = missing_last - missing_first + 1;

	probSumT = 0.0f;
	prob = aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0f);
	probSumH = aligned_vector32 < float > (HAP_NUMBER, 0.0f);
    probSumK = aligned_vector32 < float > (n_cond_haps, 0.0f);
	Alpha = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0f));
	AlphaLocus = vector < int > (segment_last - segment_first + 1, 0);
	AlphaSum = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
	AlphaSumSum = aligned_vector32 < float > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0f));
		AlphaSumMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
	}
	// Map: locus -> rel_missing index for anchors
	missing_index_by_locus.assign(locus_last - locus_first + 1, -1);
	//Cache efficient data transfer for conditioning haplotypes
	curr_rel_locus_offset = Hhap.subset(H, idxH, locus_first, locus_last);
	Hvar.allocateFast(Hhap.n_cols, Hhap.n_rows);
    Hhap.transpose(Hvar);

	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_locus = locus_first;
	prev_abs_locus = locus_first;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;

    if (super_sites && locus_to_super_idx && panel_codes && super_site_var_index) {
        ss_cond_codes = aligned_vector32<uint8_t>(n_cond_haps, 0);
        ss_emissions = aligned_vector32<float>(n_cond_haps, 1.0f);
        ss_emissions_h1 = aligned_vector32<float>(n_cond_haps, 1.0f);

        // Populate static panel code matrix (like Hvar for biallelic)
        // Take a snapshot of current cond_idx to avoid dynamic dependency
        if (cond_idx) {
            ss_panel_matrix.resize(super_sites->size());

            for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
                const SuperSite& ss = (*super_sites)[ss_idx];
                ss_panel_matrix[ss_idx].resize(n_cond_haps);

                // Unpack panel codes for all conditioning haplotypes at this supersite
                for (unsigned int k = 0; k < n_cond_haps; ++k) {
                    unsigned int gh = (*cond_idx)[k];  // Snapshot: read cond_idx ONCE here
                    ss_panel_matrix[ss_idx][k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
                }
            }
        }
    }

}

haplotype_segment_single::~haplotype_segment_single() {
	G = NULL;
	segment_first = 0;
	segment_last = 0;
	locus_first = 0;
	locus_last = 0;
	ambiguous_first = 0;
	ambiguous_last = 0;
	transition_first = 0;
	n_cond_haps = 0;
	n_missing = 0;
	curr_segment_index = 0;
	curr_segment_locus = 0;
	curr_abs_locus = 0;
	curr_rel_locus = 0;
	curr_abs_ambiguous = 0;
	curr_abs_transition = 0;
	probSumT = 0.0;
	prob.clear();
	probSumK.clear();
	probSumH.clear();
	Alpha.clear();
	AlphaSum.clear();
	AlphaSumSum.clear();
	AlphaMissing.clear();
	AlphaSumMissing.clear();
}

void haplotype_segment_single::forward() {
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;
	prev_abs_locus = locus_first;

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index);

    for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
        // Keep segment cursor in range before any Lengths_bio access (protects trailing siblings)
        if (curr_segment_index < segment_first) curr_segment_index = segment_first;
        if (curr_segment_index > segment_last) curr_segment_index = segment_last;
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool update_prev_locus = true;
		char rare_allele = M.rare_allele[curr_abs_locus];

		SiteView site_view{};
		bool has_supersite = supersites_enabled && supersite_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		if (!has_supersite) {
			bial_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		}
			const bool is_anchor = (site_view.kind == SiteKind::SuperAnchor);
			const bool is_sibling = (site_view.kind == SiteKind::SuperSibling);
			const EmitKind emit = site_view.emit_kind;
			const bool hmm_mis = (emit == EmitKind::Mis);
			const bool hmm_amb = (emit == EmitKind::Amb);
			const bool hmm_hom = (emit == EmitKind::Hom);
	// Bookkeeping flags (excluding siblings for DP paths); sibling bookkeeping
	// is handled explicitly via RUN_SIB/INIT_SIB/COLLAPSE_SIB to avoid altering
	// transition probabilities.
	const bool data_mis = hmm_mis && !is_sibling;
	const bool data_amb = hmm_amb && !is_sibling;
		yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;
        if (curr_rel_locus == 0) {
            if (is_anchor) {
                switch (emit) {
                    case EmitKind::Hom:
                        INIT_HOM();
                        break;
                    case EmitKind::Amb:
                        INIT_AMB();
                        break;
                    case EmitKind::Mis:
                        INIT_MIS();
                        break;
                }
            } else if (is_sibling) {
                // Sibling at window start: neutral init but no prev_locus advance
                INIT_MIS();
                update_prev_locus = false;
            } else {
				if (hmm_hom) INIT_HOM();
				else if (hmm_amb) INIT_AMB();
				else INIT_MIS();
			}
        } else if (curr_segment_locus != 0) {
            if (is_anchor) {
                switch (emit) {
                    case EmitKind::Hom:
                        update_prev_locus = RUN_HOM(rare_allele);
                        break;
                    case EmitKind::Amb:
                        RUN_AMB();
                        break;
                    case EmitKind::Mis:
                        RUN_MIS();
                        break;
                }
			} else if (is_sibling) {
				// Sibling within window: no-op DP but perform bookkeeping
				RUN_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) update_prev_locus = RUN_HOM(rare_allele);
				else if (hmm_amb) RUN_AMB();
				else RUN_MIS();
			}
        } else {
            if (is_anchor) {
                switch (emit) {
                    case EmitKind::Hom:
                        COLLAPSE_HOM();
                        break;
                    case EmitKind::Amb:
                        COLLAPSE_AMB();
                        break;
                    case EmitKind::Mis:
                        COLLAPSE_MIS();
                        break;
                }
			} else if (is_sibling) {
				// Sibling at segment boundary: no-op DP with bookkeeping
				COLLAPSE_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) COLLAPSE_HOM();
				else if (hmm_amb) COLLAPSE_AMB();
				else  COLLAPSE_MIS();
			}
        }
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;
		const int prev_after = prev_abs_locus;
		if (curr_segment_locus == (G->Lengths_bio[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths_bio[curr_segment_index] - 1) {
			const int rel_seg = curr_segment_index - segment_first;
			// Safety assert to capture the exact failure during testing
			assert(rel_seg >= 0 && rel_seg < static_cast<int>(Alpha.size()) && "rel_seg out of range when writing Alpha");
			Alpha[rel_seg] = prob;
			AlphaSum[rel_seg] = probSumH;
			AlphaSumSum[rel_seg] = probSumT;
			AlphaLocus[rel_seg] = prev_abs_locus;
		}
					if (data_mis) {
						// Supersite siblings are bookkeeping-only; do not cache AlphaMissing/AlphaSumMissing
						// or consume rel-missing slots for them, otherwise anchors may read sibling
						// caches during backward SC.
						if (is_sibling && supersites_enabled_flag) {
							// Do NOT advance curr_abs_missing / curr_rel_missing
						} else {
							AlphaMissing[curr_rel_missing] = prob;
							AlphaSumMissing[curr_rel_missing] = probSumH;
							// Record the rel-missing index for backward SC (supersites) or IMPUTE (biallelic)
							int idx = curr_abs_locus - locus_first;
							if (idx >= 0 && idx < (int)missing_index_by_locus.size()) missing_index_by_locus[idx] = curr_rel_missing;
							curr_abs_missing ++;
						}
					}
		// Only increment segment locus for non-sibling loci (siblings don't count toward segment length)
		if (!is_sibling) {
			curr_segment_locus ++;
		}
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
		const bool can_advance_amb = data_amb && has_amb_range && (curr_abs_ambiguous <= ambiguous_last);
		if (can_advance_amb) {
			curr_abs_ambiguous++;
		}
		if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last + 1)) {
			assert(false && "forward ambiguous cursor moved out of window bounds");
		}
			if (curr_segment_locus >= G->Lengths_bio[curr_segment_index]) {
				// Avoid advancing past the last segment when trailing siblings are present
				if (curr_segment_index < segment_last) {
					curr_segment_index++;
					curr_segment_locus = 0;
				} else {
					// Pin to last valid slot to keep Lengths_bio accesses in range
					curr_segment_locus = G->Lengths_bio[curr_segment_index] - 1;
				}
			}

	}
		if (curr_segment_index > segment_last) {
			curr_segment_index = segment_last;
			curr_segment_locus = (segment_last >= segment_first) ? G->Lengths_bio[segment_last] - 1 : 0;
		}
		curr_abs_locus = locus_last;
}

int haplotype_segment_single::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities, 
                                       vector<float>* SC, const vector<bool>* anchor_has_missing, const vector<uint32_t>* supersite_sc_offset) {
	int n_underflow_recovered = 0;
	// Set thread-local offset storage for IMPUTE_SUPERSITE_MULTIVARIATE calls
	this->supersite_sc_offset = supersite_sc_offset;
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths_bio[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;
	prev_abs_locus = locus_last;
	const bool has_transition_buffer = !transition_probabilities.empty();
	const size_t expected_transitions = static_cast<size_t>(G->n_transitions);
	if (!has_transition_buffer) {
		assert(expected_transitions == 0);
	} else {
		assert(transition_probabilities.size() >= expected_transitions);
	}

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index);

	// Flag: backward pass always starts with initialization; if locus_last is a sibling, defer until first non-sibling
	bool need_init = true;
	bool pending_collapse = false;
	yt = 0.0f;
	nt = 1.0f;

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		if (curr_segment_locus == G->Lengths_bio[curr_segment_index] - 1) pending_collapse = true;
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		char rare_allele = M.rare_allele[curr_abs_locus];
		bool update_prev_locus = true;

		SiteView site_view{};
		bool has_supersite = supersites_enabled && supersite_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		if (!has_supersite) {
			bial_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		}
		const bool is_anchor = (site_view.kind == SiteKind::SuperAnchor);
		const bool is_sibling = (site_view.kind == SiteKind::SuperSibling);

		// Deferred initialization: if we started on a sibling, keep doing INIT_MIS until we hit a non-sibling
		if (need_init && is_sibling) {
			// Minimal interaction: just advance the counters as if this locus was processed, but do no HMM math.
			// This ensures that when the first anchor is reached, prev_abs_locus is equal to it, making yt=0,
			// which is correct for the start of a backward pass.
			INIT_SIB(site_view);
			// Siblings don't count toward segment length, so don't decrement curr_segment_locus
			// NOTE: No SET_OTHER_TRANS here - deferred init doesn't write transitions
			continue;
		}

		const EmitKind emit = site_view.emit_kind;
		const bool hmm_mis = (emit == EmitKind::Mis);
		const bool hmm_amb = (emit == EmitKind::Amb);
		const bool hmm_hom = (emit == EmitKind::Hom);
	// Do not treat supersite siblings as missing/ambiguous here; the anchor handles
	// imputation for the whole supersite. Sibling positions are bookkeeping-only.
	const bool data_mis = hmm_mis && !is_sibling;
	const bool data_amb = hmm_amb && !is_sibling;
		// Guarded transition probability calculation
		if (!is_sibling) {
			yt = (curr_abs_locus == locus_last || curr_abs_locus == prev_abs_locus) ? 0.0f : M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
			nt = 1.0f - yt;
		} else {
		}

		// Deferred initialization: first non-sibling after starting on a sibling - use INIT
		if (need_init && is_anchor) {
			need_init = false;
			pending_collapse = false;
			// FIX: Set prev_abs_locus to curr_abs_locus to ensure yt=0 for the first anchor
			prev_abs_locus = curr_abs_locus;
			switch (emit) {
				case EmitKind::Hom:
					INIT_HOM();
					break;
				case EmitKind::Amb:
					INIT_AMB();
					break;
				case EmitKind::Mis:
					INIT_MIS();
					break;
			}
		} else if (need_init && !is_sibling) {
			// Deferred initialization: first biallelic after starting on a sibling - use INIT
			need_init = false;
			pending_collapse = false;
			if (hmm_hom) INIT_HOM();
			else if (hmm_amb) INIT_AMB();
			else INIT_MIS();
		} else if (!pending_collapse) {
			if (is_anchor) {
				switch (emit) {
					case EmitKind::Hom:
						update_prev_locus = RUN_HOM(rare_allele);
						break;
					case EmitKind::Amb:
						RUN_AMB();
						break;
					case EmitKind::Mis:
						RUN_MIS();
						break;
				}
			} else if (is_sibling) {
				// Sibling within window (backward): apply transition only
				RUN_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) update_prev_locus = RUN_HOM(rare_allele);
				else if (hmm_amb) RUN_AMB();
				else RUN_MIS();
			}
		} else {
			if (is_anchor) {
				switch (emit) {
					case EmitKind::Hom:
						COLLAPSE_HOM();
						break;
					case EmitKind::Amb:
						COLLAPSE_AMB();
						break;
					case EmitKind::Mis:
						COLLAPSE_MIS();
						break;
				}
				pending_collapse = false;
			} else if (is_sibling) {
				// Sibling at segment boundary (backward): bookkeeping only
				COLLAPSE_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) COLLAPSE_HOM();
				else if (hmm_amb) COLLAPSE_AMB();
				else COLLAPSE_MIS();
				pending_collapse = false;
			}
		}
		if (curr_segment_locus == 0) SUMK();
		const int prev_after = update_prev_locus ? curr_abs_locus : prev_abs_locus;
		prev_abs_locus = prev_after;

			// Only emit transition probabilities if caller provided storage
			if (!transition_probabilities.empty()) {
				if (curr_abs_locus == 0) {
					SET_FIRST_TRANS(transition_probabilities);
				}
				if (!is_sibling && curr_segment_locus == 0 && curr_abs_locus != locus_first && curr_segment_index > segment_first) {
					int ret = SET_OTHER_TRANS(transition_probabilities);
					if (ret < 0) return ret;
					else n_underflow_recovered += ret;
				}
			}

		if (data_mis) {
			bool handled = false;
			// First, special-case supersite anchors that have SC posteriors prepared
			if (is_anchor && anchor_has_missing && SC && site_view.supersite_index >= 0 && (*anchor_has_missing)[site_view.supersite_index]) {
				int idx = -1;
				int map_i = curr_abs_locus - locus_first;
				if (map_i >= 0 && map_i < (int)missing_index_by_locus.size()) idx = missing_index_by_locus[map_i];
				if (idx >= 0) {
					// Ensure supersite codes are loaded before any debug numerators
					ss_load_cond_codes(*site_view.supersite, site_view.supersite_index);
					IMPUTE_SUPERSITE_MULTIVARIATE(*SC, *site_view.supersite, site_view.supersite_index, idx, supersite_sc_offset);
					handled = true;
				}
			}

			// For non-supersite or regular missing sites, only call IMPUTE if
			// forward pass actually recorded a missing index for this locus.
			if (!handled) {
				int map_i = curr_abs_locus - locus_first;
				int rel_missing_idx = -1;
				if (map_i >= 0 && map_i < (int)missing_index_by_locus.size()) rel_missing_idx = missing_index_by_locus[map_i];
				if (rel_missing_idx >= 0) {
					// Forward recorded a missing slot for this locus -> perform imputation
					IMPUTE(missing_probabilities);
				} else {
					// Nothing recorded in forward (e.g., sibling or no-missing window) -> skip IMPUTE
				}
			}
			// Only decrement missing index for non-siblings in supersite mode
			// Forward pass skips incrementing for siblings (lines 724-730), so backward must match
			if (!(is_sibling && supersites_enabled_flag)) {
				curr_abs_missing--;
			}
		}


		// Only decrement segment locus for non-sibling loci (siblings don't count toward segment length)
		if (!is_sibling) {
			curr_segment_locus--;
		}
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
		const bool can_retreat_amb = data_amb && has_amb_range && (curr_abs_ambiguous >= ambiguous_first);
		if (can_retreat_amb) curr_abs_ambiguous--;
		const int lower_bound = ambiguous_first - 1;
		if (has_amb_range && (curr_abs_ambiguous < lower_bound || curr_abs_ambiguous > ambiguous_last)) {
			assert(false && "backward ambiguous cursor moved out of window bounds");
		}
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths_bio[curr_segment_index] - 1;
		}
	}
		curr_segment_index = segment_first;
		curr_segment_locus = 0;
		curr_abs_locus = locus_first;
	return n_underflow_recovered;
}

haplotype_segment_single::LaneWeights haplotype_segment_single::compute_lane_match_weights(const SiteView& site_view) const {
	LaneWeights lw{};
	for (int h = 0; h < HAP_NUMBER; ++h) lw.match[h] = 0.0f;
	lw.neither = 0.0f;

	for (unsigned int k = 0; k < n_cond_haps; ++k) {
		const uint8_t donor_code = (site_view.kind == SiteKind::SuperAnchor || site_view.kind == SiteKind::SuperSibling)
			? ss_cond_codes[k]
			: Hvar.get(locus_first, k);
		for (int h = 0; h < HAP_NUMBER; ++h) {
			const uint8_t want = lane_expected_class(site_view, h);
			if (donor_code == want) {
				lw.match[h] += 1.0f;
			} else {
				lw.neither += 1.0f;
			}
		}
	}
	return lw;
}

void haplotype_segment_single::build_lane_priors_first(const SiteView& site_view, double lane_probs[HAP_NUMBER]) const {
	LaneWeights lw = compute_lane_match_weights(site_view);

	double total = 0.0;
	const double neither_share = (HAP_NUMBER > 0) ? (0.5 * static_cast<double>(lw.neither) / static_cast<double>(HAP_NUMBER)) : 0.0;
	for (int h = 0; h < HAP_NUMBER; ++h) {
		double w = static_cast<double>(lw.match[h]) + neither_share;
		lane_probs[h] = w;
		total += w;
	}
	if (total <= std::numeric_limits<double>::min()) {
		for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] = 1.0 / static_cast<double>(HAP_NUMBER);
	} else {
		const double inv = 1.0 / total;
		for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] *= inv;
	}
}

void haplotype_segment_single::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	const unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	std::vector<double> cprobs(n_transitions, 0.0);

	// Loudly enforce buffer sizing instead of skipping writes.
	assert(n_transitions > 0);
	assert(transition_probabilities.size() >= n_transitions);

	double lane_probs[HAP_NUMBER];

	// Build SiteView for the first locus to derive lane_class
	SiteView sv{};
	bool has_supersite = false;
	if (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx) {
		SupersiteEmissionAdapter sup_adapt(G, super_sites, locus_to_super_idx, super_site_var_index);
		has_supersite = sup_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	if (!has_supersite) {
		BiallelicEmissionAdapter bial_adapt(G, &Hvar);
		bial_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	build_lane_priors_first(sv, lane_probs);

	double scaleDip = 0.0;
	for (unsigned int d = 0, t = 0; d < 64; ++d) {
		if (DIP_GET(G->Diplotypes[0], d)) {
			const int hap0 = DIP_HAP0(d);
			const int hap1 = DIP_HAP1(d);
			const double val = lane_probs[hap0] * lane_probs[hap1];
			cprobs[t++] = val;
			scaleDip += val;
		}
	}

	if (scaleDip <= std::numeric_limits<double>::min()) {
		const double uniform = 1.0 / static_cast<double>(n_transitions);
		for (unsigned int t = 0; t < n_transitions; ++t) transition_probabilities[t] = uniform;
		return;
	}

	const double inv_scale = 1.0 / scaleDip;
	for (unsigned int t = 0; t < n_transitions; ++t) transition_probabilities[t] = cprobs[t] * inv_scale;
}

int haplotype_segment_single::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
    int underflow_recovered = 0;
    if (TRANS_HAP()) {
        return -1;
    }

    if (TRANS_DIP_MULT()) {
        if (TRANS_DIP_ADD()) return -2;
        else underflow_recovered = 1;
    }

	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;
	// Loudly enforce buffer sizing instead of silently skipping writes.
	assert(n_transitions > 0);
	const int start = curr_abs_transition - static_cast<int>(n_transitions - 1);
	assert(start >= 0);
	assert(static_cast<size_t>(start + n_transitions) <= transition_probabilities.size());

	curr_abs_transition = start;
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}

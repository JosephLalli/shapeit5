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

#include <models/haplotype_segment_double.h>
#include <models/site_emission_adapter.h>
#include <cstdio>
#include <string>
#include <limits>
#include <algorithm>

using namespace std;

haplotype_segment_double::haplotype_segment_double(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M) :
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

	probSumT = 0.0;
	prob = aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0);
	probSumH = aligned_vector32 < double > (HAP_NUMBER, 0.0);
	probSumK = aligned_vector32 < double > (n_cond_haps, 0.0);
	Alpha = vector < aligned_vector32 < double > > (segment_last - segment_first + 1, aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0));
	AlphaLocus = vector < int > (segment_last - segment_first + 1, 0);
	AlphaSum = vector < aligned_vector32 < double > > (segment_last - segment_first + 1, aligned_vector32 < double > (HAP_NUMBER, 0.0));
	AlphaSumSum = aligned_vector32 < double > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0));
		AlphaSumMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER, 0.0));
	}
	//Cache efficient data transfer for conditioning haplotypes
	curr_rel_locus_offset = Hhap.subset(H, idxH, locus_first, locus_last);
	Hvar.allocateFast(Hhap.n_cols, Hhap.n_rows);
	Hhap.transpose(Hvar);

	// allocate scratch for supersites if enabled
	if (supersites_enabled_flag) {
		ss_cond_codes = aligned_vector32<uint8_t>(n_cond_haps, 0);
		ss_emissions = aligned_vector32<double>(n_cond_haps, 1.0);
		ss_emissions_h1 = aligned_vector32<double>(n_cond_haps, 1.0);

		if (panel_codes) {
			ss_panel_matrix.resize(super_sites->size());
			for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
				const SuperSite& ss = (*super_sites)[ss_idx];
				ss_panel_matrix[ss_idx].resize(n_cond_haps);
				for (unsigned int k = 0; k < n_cond_haps; ++k) {
					unsigned int gh = (*cond_idx)[k];
					ss_panel_matrix[ss_idx][k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
				}
			}
		}
	}

}

haplotype_segment_double::~haplotype_segment_double() {
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

void haplotype_segment_double::forward() {
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;
	prev_abs_locus = locus_first;
	// Use biological lengths (anchors only) when available so sibling loci do not skew
	// segment bookkeeping; mirrors the single-precision path.
	const std::vector<unsigned short>& lengths_bio = G->Lengths_bio.empty() ? G->Lengths : G->Lengths_bio;
	assert(!lengths_bio.empty());

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
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb && !is_sibling;

		// Guard against non-monotone prev_abs_locus (can happen with sibling bookkeeping
		// when tests skip filling Lengths_bio): fall back to yt=0 in that case.
		if (curr_abs_locus <= prev_abs_locus) {
			yt = 0.0;
		} else {
			yt = M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		}
		nt = 1.0 - yt;

            if (curr_rel_locus == 0) {
                if (is_anchor) {
					switch (emit) {
						case EmitKind::Hom: INIT_HOM(); break;
						case EmitKind::Amb: INIT_AMB(); break;
						case EmitKind::Mis: INIT_MIS(); break;
					}
			} else if (is_sibling) {
				// Sibling at window start: initialize neutrally but do not advance prev_abs_locus
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
				// Sibling within window: no-op propagation (avoid renormalization)
				update_prev_locus = false;
			} else {
				if (hmm_hom) update_prev_locus = RUN_HOM(rare_allele);
				else if (hmm_amb) RUN_AMB();
				else RUN_MIS();
			}
            } else {
                if (is_anchor) {
					switch (emit) {
						case EmitKind::Hom: COLLAPSE_HOM(); break;
						case EmitKind::Amb: COLLAPSE_AMB(); break;
						case EmitKind::Mis: COLLAPSE_MIS(); break;
					}
			} else if (is_sibling) {
				// Sibling at segment boundary: no-op and do not advance prev_abs_locus
				update_prev_locus = false;
			} else {
				if (hmm_hom) COLLAPSE_HOM();
				else if (hmm_amb) COLLAPSE_AMB();
				else  COLLAPSE_MIS();
			}
		}
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

        if (curr_segment_locus == (lengths_bio[curr_segment_index] - 1)) SUMK();

	if (curr_segment_locus == lengths_bio[curr_segment_index] - 1) {
			const int rel_seg = curr_segment_index - segment_first;
			Alpha[rel_seg] = prob;
			AlphaSum[rel_seg] = probSumH;
			AlphaSumSum[rel_seg] = probSumT;
			AlphaLocus[rel_seg] = prev_abs_locus;
	}
        if (data_mis) {
            if (curr_rel_missing < 0 || curr_rel_missing >= (int)AlphaMissing.size()) {
                assert(false && "AlphaMissing index out of bounds");
            }
            AlphaMissing[curr_rel_missing] = prob;
            AlphaSumMissing[curr_rel_missing] = probSumH;
            if (is_anchor) {
                int map_i = curr_abs_locus - locus_first;
                if (map_i >= 0 && map_i < (int)missing_index_by_locus.size()) missing_index_by_locus[map_i] = curr_rel_missing;
            }
            curr_abs_missing ++;
        }

	// Only increment segment locus for non-sibling loci to mirror single-precision
	if (!is_sibling) curr_segment_locus ++;
	const bool has_amb_range = (ambiguous_first <= ambiguous_last);
	const bool can_advance_amb = data_amb && has_amb_range && (curr_abs_ambiguous <= ambiguous_last);
	if (can_advance_amb) {
		curr_abs_ambiguous++;
	}
	if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last + 1)) {
		assert(false && "forward ambiguous cursor moved out of window bounds");
	}
		if (curr_segment_locus >= lengths_bio[curr_segment_index]) {
			// Avoid advancing past the last segment when trailing siblings are present
			if (curr_segment_index < segment_last) {
				curr_segment_index++;
				curr_segment_locus = 0;
			} else {
				// Pin to last valid slot to keep Lengths_bio accesses in range
				curr_segment_locus = lengths_bio[curr_segment_index] - 1;
			}
		}

	}
}

int haplotype_segment_double::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities,
                                       vector < float > * SC, const vector < bool > * anchor_has_missing, const vector<uint32_t>* supersite_sc_offset) {
	int n_underflow_recovered = 0;
	// Use biological lengths (anchors only) when available so sibling loci do not skew
	// segment bookkeeping; mirrors the single-precision path.
	const std::vector<unsigned short>& lengths_bio = G->Lengths_bio.empty() ? G->Lengths : G->Lengths_bio;
	assert(!lengths_bio.empty());
	// Set thread-local offset storage for IMPUTE_SUPERSITE_MULTIVARIATE calls
	this->supersite_sc_offset = supersite_sc_offset;
	curr_segment_index = segment_last;
	curr_segment_locus = lengths_bio[segment_last] - 1;
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

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		if (curr_segment_locus == lengths_bio[curr_segment_index] - 1) pending_collapse = true;
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
		// Deferred initialization: if we started on a sibling, keep doing INIT_SIB until we hit a non-sibling
		if (need_init && is_sibling) {
			INIT_SIB(site_view);
			update_prev_locus = false;
			prev_abs_locus = update_prev_locus ? curr_abs_locus : prev_abs_locus;
			continue;
		}
		const EmitKind emit = site_view.emit_kind;
		const bool hmm_mis = (emit == EmitKind::Mis);
		const bool hmm_amb = (emit == EmitKind::Amb);
		const bool hmm_hom = (emit == EmitKind::Hom);
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb && !is_sibling;

		if (curr_abs_locus >= prev_abs_locus) {
			yt = 0.0;
		} else {
			yt = M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		}
		nt = 1.0f - yt;

		// Deferred initialization: first non-sibling after starting on a sibling - use INIT
		if (need_init && is_anchor) {
			need_init = false;
			pending_collapse = false;
			switch (emit) {
				case EmitKind::Hom: INIT_HOM(); break;
				case EmitKind::Amb: INIT_AMB(); break;
				case EmitKind::Mis: INIT_MIS(); break;
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
				// Sibling within window (backward): no-op propagation
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
					case EmitKind::Hom: COLLAPSE_HOM(); break;
					case EmitKind::Amb: COLLAPSE_AMB(); break;
					case EmitKind::Mis: COLLAPSE_MIS(); break;
				}
				pending_collapse = false;
			} else if (is_sibling) {
				// Sibling at segment boundary (backward): no-op
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
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

		if (curr_abs_locus == 0) {
			SET_FIRST_TRANS(transition_probabilities);
		}
		if (!is_sibling && curr_segment_locus == 0 && curr_abs_locus != locus_first && curr_segment_index > segment_first) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (data_mis) {
			bool supersite_handled = false;
			if (is_anchor && anchor_has_missing && SC && site_view.supersite_index >= 0 && (*anchor_has_missing)[site_view.supersite_index]) {
				int map_i = curr_abs_locus - locus_first;
				int rel_idx = -1;
				if (map_i >= 0 && map_i < (int)missing_index_by_locus.size()) rel_idx = missing_index_by_locus[map_i];
				if (rel_idx >= 0) {
					IMPUTE_SUPERSITE_MULTIVARIATE(*SC, *site_view.supersite, site_view.supersite_index, rel_idx, supersite_sc_offset);
					supersite_handled = true;
				}
			}
			if (!supersite_handled) {
				IMPUTE(missing_probabilities);
			}
			curr_abs_missing--;
		}

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
			curr_segment_locus = lengths_bio[curr_segment_index] - 1;
		}
	}

	// DIAGNOSTIC: Check if transitions_stored matches expected for THIS WINDOW only.
	return n_underflow_recovered;
}

void haplotype_segment_double::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	const unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	std::vector<double> cprobs(n_transitions, 0.0);

	// Loudly enforce buffer sizing instead of skipping writes.
	assert(n_transitions > 0);
	assert(transition_probabilities.size() >= n_transitions);

	double lane_probs[HAP_NUMBER];
	const double lane_total = probSumT;
	const double inv_total = (lane_total > std::numeric_limits<double>::min()) ? (1.0 / lane_total) : 0.0;
	for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] = probSumH[h] * inv_total;

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

int haplotype_segment_double::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
    int underflow_recovered = 0;
    if (TRANS_HAP()) {
        return -1;
    }
    if (TRANS_DIP_MULT()) {
        bool under = TRANS_DIP_ADD();
        if (under) {
            return -2;
        } else underflow_recovered = 1;
    }
	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;

	// Loudly enforce buffer sizing instead of skipping writes.
	const long long start_idx = static_cast<long long>(curr_abs_transition) - static_cast<long long>(n_transitions) + 1;
	const size_t buf_size = transition_probabilities.size();
	assert(!transition_probabilities.empty());
	assert(start_idx >= 0);
	assert(static_cast<size_t>(start_idx) + n_transitions <= buf_size);

	curr_abs_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}

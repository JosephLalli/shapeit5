// Sibling bookkeeping: update trace, maintain cursor/index sanity, no DP math
void haplotype_segment_single::handle_sibling_bookkeeping(const SiteView& site_view) {
	// Only log ambiguous cursor, do not advance or store AlphaMissing/AlphaSumMissing
	// This keeps ambiguous/missing indices aligned for anchors, but does not touch DP state
	trace_ambiguous_cursor(trace_forward_active ? "fwd_post" : "bwd_post", curr_abs_locus, true, 0);
	// Optionally, add asserts here to ensure ambiguous/missing cursor is in bounds
	assert(curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last);
	assert(curr_abs_missing >= missing_first && curr_abs_missing <= missing_last);
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
#include <mutex>
#include <cstdio>
#include <string>
#include <limits>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

// Debug underflow tracing (gated by SHAPEIT5_DEBUG_UNDERFLOW)
namespace {
static std::mutex g_underflow_log_mutex;
static bool underflow_trace_enabled() {
    static int flag = -1;
    if (flag < 0) {
        const char* env = std::getenv("SHAPEIT5_DEBUG_UNDERFLOW");
        flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return flag == 1;
}
static void ensure_logs_dir() {
    struct stat st{};
    if (stat("logs", &st) != 0) {
        mkdir("logs", 0777);
    }
}

static bool supersite_trace_enabled() {
    static int flag = -1;
    if (flag < 0) {
        const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
        flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return flag == 1;
}

} // namespace

void haplotype_segment_single::trace_ambiguous_cursor(const char* stage, int locus, bool is_sibling, int expected_delta) const {
	if (!supersite_trace_enabled()) return;
	if (!stage) stage = "";
	if (ambiguous_first > ambiguous_last) return;

	const bool forward_stage = (std::strcmp(stage, "fwd_pre") == 0) || (std::strcmp(stage, "fwd_post") == 0);
	const bool backward_stage = (std::strcmp(stage, "bwd_pre") == 0) || (std::strcmp(stage, "bwd_post") == 0);

	int actual_delta = 0;
	bool delta_valid = false;

	if (forward_stage && trace_forward_active) {
		if (std::strcmp(stage, "fwd_pre") == 0) {
			trace_forward_pre_cursor = curr_abs_ambiguous;
			trace_forward_pre_locus = locus;
			trace_forward_pre_valid = true;
		} else if (std::strcmp(stage, "fwd_post") == 0) {
			if (trace_forward_pre_valid && trace_forward_pre_locus == locus) {
				actual_delta = curr_abs_ambiguous - trace_forward_pre_cursor;
				delta_valid = true;
			}
			trace_forward_pre_valid = false;
		}
	} else if (backward_stage && trace_backward_active) {
		if (std::strcmp(stage, "bwd_pre") == 0) {
			trace_backward_pre_cursor = curr_abs_ambiguous;
			trace_backward_pre_locus = locus;
			trace_backward_pre_valid = true;
		} else if (std::strcmp(stage, "bwd_post") == 0) {
			if (trace_backward_pre_valid && trace_backward_pre_locus == locus) {
				actual_delta = curr_abs_ambiguous - trace_backward_pre_cursor;
				delta_valid = true;
			}
			trace_backward_pre_valid = false;
		}
	} else {
		if (std::strcmp(stage, "bwd_pre") == 0) {
			trace_backward_pre_valid = false;
		} else if (std::strcmp(stage, "fwd_pre") == 0) {
			trace_forward_pre_valid = false;
		}
	}

	if (delta_valid && actual_delta != expected_delta) {
		std::fprintf(stdout,
				 "[ss-amb-delta][single] stage=%s locus=%d curr_abs_ambiguous=%d expected_delta=%d actual_delta=%d seg_idx=%d is_sibling=%d\n",
				 stage,
				 locus,
				 curr_abs_ambiguous,
				 expected_delta,
				 actual_delta,
				 curr_segment_index,
				 static_cast<int>(is_sibling));
		assert(actual_delta == expected_delta && "ambiguous cursor delta mismatch");
	}

	const int lower = ambiguous_first;
	const int upper = ambiguous_last;
	if (curr_abs_ambiguous < lower || curr_abs_ambiguous > upper) {
		std::fprintf(stdout,
				 "[ss-amb-drift][single] stage=%s locus=%d curr_abs_ambiguous=%d range=[%d,%d] rel=%d seg_idx=%d is_sibling=%d\n",
				 stage,
				 locus,
				 curr_abs_ambiguous,
				 lower,
				 upper,
				 curr_abs_ambiguous - lower,
				 curr_segment_index,
				 static_cast<int>(is_sibling));
	}

	std::string delta_buffer;
	const char* delta_repr = "NA";
	if (delta_valid) {
		delta_buffer = std::to_string(actual_delta);
		delta_repr = delta_buffer.c_str();
	}
	std::fprintf(stdout,
			 "[ss-amb-cursor][single] stage=%s locus=%d curr_abs_ambiguous=%d range=[%d,%d] actual_delta=%s expected_delta=%d seg_idx=%d is_sibling=%d\n",
			 stage,
			 locus,
			 curr_abs_ambiguous,
			 lower,
			 upper,
			 delta_repr,
			 expected_delta,
			 curr_segment_index,
			 static_cast<int>(is_sibling));
}

haplotype_segment_single::haplotype_segment_single(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M,
    const std::vector<SuperSite>* _super_sites,
    const std::vector<bool>* _is_super_site,
    const std::vector<int>* _locus_to_super_idx,
	const uint8_t* _panel_codes,
	size_t _panel_codes_size,
	const std::vector<int>* _super_site_var_index) :
	G(_G), M(_M), super_sites(_super_sites), is_super_site(_is_super_site),
	locus_to_super_idx(_locus_to_super_idx), panel_codes(_panel_codes), panel_codes_size(_panel_codes_size), super_site_var_index(_super_site_var_index), cond_idx(&idxH),
    supersite_sc_offset(nullptr),
    supersites_enabled_flag(_super_sites && _locus_to_super_idx && _super_site_var_index && _panel_codes) {
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
	AlphaLaneSum = vector<LaneMarginal>(segment_last - segment_first + 1, LaneMarginal{});
	AlphaSumSum = aligned_vector32 < float > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER * n_cond_haps, 0.0f));
		AlphaSumMissing = vector < aligned_vector32 < float > > (n_missing, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
	}
	// Map: locus -> rel_missing index for anchors
	missing_index_by_locus.assign(locus_last - locus_first + 1, -1);
	init_match_mask.resize(static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER);
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
        // Initialize cache tracking - sized to number of supersites, all false
        // Cache is per-segment and automatically reset when new segment created
        // (segments are created fresh for each window, so no explicit invalidation needed)
        if (super_sites) {
            ss_cached.resize(super_sites->size(), false);
        }
    }

	trace_forward_pre_cursor = 0;
	trace_forward_pre_locus = -1;
	trace_forward_pre_valid = false;
	trace_backward_pre_cursor = 0;
	trace_backward_pre_locus = -1;
	trace_backward_pre_valid = false;
	trace_forward_active = false;
	trace_backward_active = false;

    // Test-only diagnostics
    if (supersite_trace_enabled()) {
        std::fprintf(stdout,
            "HS_ctor seg=[%d,%d] loci=[%d,%d] n_cond_haps=%u n_missing=%u supersites=%s panel_codes=%s cond_idx=%s\n",
            segment_first, segment_last, locus_first, locus_last, n_cond_haps, n_missing,
            (super_sites?"Y":"N"), (panel_codes?"Y":"N"), (cond_idx?"Y":"N"));
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
	AlphaLaneSum.clear();
	AlphaSumSum.clear();
	AlphaMissing.clear();
	AlphaSumMissing.clear();
}

void haplotype_segment_single::forward() {
    if (supersite_trace_enabled()) std::fprintf(stderr, "FWD0\n");
	curr_segment_index = segment_first;
	curr_segment_locus = 0;
	curr_abs_ambiguous = ambiguous_first;
	curr_abs_missing = missing_first;
	prev_abs_locus = locus_first;
	if (supersite_trace_enabled()) std::fprintf(stderr, "FWD1 seg=%d lf=%d ll=%d\n", curr_segment_index, locus_first, locus_last);
	trace_forward_active = true;
	trace_backward_active = false;

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx, panel_codes_size);
    if (supersite_trace_enabled()) std::fprintf(stderr, "FWD2 after adapters\n");

    if (supersite_trace_enabled()) {
        std::fprintf(stdout,
            "FWD_start loci=[%d,%d] seg=[%d,%d] n_cond_haps=%u supersites_enabled=%d prob_size=%zu probSumH=%zu\n",
            locus_first, locus_last, segment_first, segment_last, n_cond_haps, (int)supersites_enabled,
            prob.size(), (size_t)HAP_NUMBER);
    }

    for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
        if (supersite_trace_enabled()) std::fprintf(stderr, "FWD3 loop enter abs=%d rel_off=%d\n", curr_abs_locus, curr_rel_locus_offset);
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool update_prev_locus = true;
		char rare_allele = M.rare_allele[curr_abs_locus];

		SiteView site_view{};
        if (supersite_trace_enabled()) std::fprintf(stderr, "FWD4 before build_view supersites_enabled=%d\n", (int)supersites_enabled);
		bool has_supersite = supersites_enabled && supersite_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
        if (supersite_trace_enabled()) std::fprintf(stderr, "FWD5 after build_view has_ss=%d kind=%d emit=%d\n", (int)has_supersite, (int)site_view.kind, (int)site_view.emit_kind);
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
		trace_ambiguous_cursor("fwd_pre", curr_abs_locus, is_sibling, 0);

		if (supersite_trace_enabled()) {
			const bool dbg_is_anchor = (site_view.kind == SiteKind::SuperAnchor);
			const bool dbg_is_sibling = (site_view.kind == SiteKind::SuperSibling);
			std::fprintf(stdout,
				"FWD locus=%d rel=%d seg_locus=%d kind=%d emit=%d is_anchor=%d is_sibling=%d amb=%d mis=%d sample_cls0=%u sample_cls1=%u amb_mask=0x%02x curr_abs_amb=%d range=[%d,%d]\n",
				curr_abs_locus, curr_rel_locus, curr_segment_locus,
				(int)site_view.kind, (int)emit,
				dbg_is_anchor ? 1 : 0,
				dbg_is_sibling ? 1 : 0,
				(int)hmm_amb, (int)hmm_mis,
				site_view.sample_class0, site_view.sample_class1,
				site_view.amb_mask,
				curr_abs_ambiguous,
				ambiguous_first,
				ambiguous_last);
		}

        if (curr_rel_locus == 0) {
            if (is_anchor) {
                if (M.ss_anchor_split_emissions) {
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                    INIT_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
                } else {
                    switch (emit) {
                        case EmitKind::Hom:
                            SS_INIT_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
                            break;
                        case EmitKind::Amb:
                            SS_INIT_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
                            break;
                        case EmitKind::Mis:
                            SS_INIT_MIS();
                            break;
                    }
                }
            } else if (is_sibling) {
                // Sibling at window start: neutral init but no prev_locus advance
                INIT_MIS();
                update_prev_locus = false;
            } else {
				if (hmm_mis) {
					INIT_MIS();
				} else {
					bial_adapter.build_match_mask(site_view, n_cond_haps, curr_rel_locus + curr_rel_locus_offset, init_match_mask);
					INIT_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
				}
			}
        } else if (curr_segment_locus != 0) {
            if (is_anchor) {
                if (M.ss_anchor_split_emissions) {
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                    RUN_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
                } else {
                    switch (emit) {
                        case EmitKind::Hom:
                            update_prev_locus = SS_RUN_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
                            break;
                        case EmitKind::Amb:
                            update_prev_locus = SS_RUN_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
                            break;
                        case EmitKind::Mis:
                            update_prev_locus = SS_RUN_MIS();
                            break;
                    }
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
                if (M.ss_anchor_split_emissions) {
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                    COLLAPSE_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
                } else {
                    switch (emit) {
                        case EmitKind::Hom:
                            SS_COLLAPSE_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
                            break;
                        case EmitKind::Amb:
                            SS_COLLAPSE_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
                            break;
                        case EmitKind::Mis:
                            SS_COLLAPSE_MIS();
                            break;
                    }
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

		if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			const int rel_seg = curr_segment_index - segment_first;
			Alpha[rel_seg] = prob;
			AlphaSum[rel_seg] = probSumH;
			__m256 lane_vec = _mm256_load_ps(&probSumH[0]);
			_mm256_store_ps(AlphaLaneSum[rel_seg].lane, lane_vec);
			AlphaSumSum[rel_seg] = probSumT;
			AlphaLocus[rel_seg] = prev_abs_locus;
		}
		if (data_mis) {
			AlphaMissing[curr_rel_missing] = prob;
			AlphaSumMissing[curr_rel_missing] = probSumH;
			// If this is a supersite anchor, record the rel-missing index for backward SC
			if (is_anchor) {
				int idx = curr_abs_locus - locus_first;
				if (idx >= 0 && idx < (int)missing_index_by_locus.size()) missing_index_by_locus[idx] = curr_rel_missing;
			}
			curr_abs_missing ++;
		}

		curr_segment_locus ++;
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
		const int cursor_before = curr_abs_ambiguous;
	const bool can_advance_amb = data_amb && has_amb_range && (curr_abs_ambiguous < ambiguous_last);
	const bool sib_advance_amb = is_sibling && hmm_amb && has_amb_range && (curr_abs_ambiguous < ambiguous_last);
	const int expected_delta = (can_advance_amb || sib_advance_amb) ? 1 : 0;
		if (supersite_trace_enabled()) {
			std::fprintf(stdout,
				"FWD.delta locus=%d cursor_before=%d expected_delta=%d data_amb=%d data_mis=%d is_sibling=%d has_range=%d\n",
				curr_abs_locus,
				cursor_before,
				expected_delta,
				data_amb ? 1 : 0,
				data_mis ? 1 : 0,
				is_sibling ? 1 : 0,
				has_amb_range ? 1 : 0);
		}
	if (can_advance_amb || sib_advance_amb) curr_abs_ambiguous++;
		trace_ambiguous_cursor("fwd_post", curr_abs_locus, is_sibling, expected_delta);
		if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last)) {
			if (supersite_trace_enabled()) {
				int seg_len = (curr_segment_index >= 0 && curr_segment_index < (int)G->Lengths.size()) ? G->Lengths[curr_segment_index] : -1;
				std::fprintf(stderr,
					"[ss-amb-oob][single] stage=fwd_post locus=%d before=%d after=%d expected_delta=%d range=[%d,%d] seg_idx=%d seg_len=%d seg_loc=%d is_sibling=%d\n",
					curr_abs_locus,
					cursor_before,
					curr_abs_ambiguous,
					expected_delta,
					ambiguous_first,
					ambiguous_last,
					curr_segment_index,
					seg_len,
					curr_segment_locus,
					static_cast<int>(is_sibling));
				std::fflush(stderr);
			}
			assert(false && "forward ambiguous cursor moved out of window bounds");
		}
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
		if (curr_segment_index > segment_last) {
			curr_segment_index = segment_last;
			curr_segment_locus = (segment_last >= segment_first) ? G->Lengths[segment_last] - 1 : 0;
		}
		curr_abs_locus = locus_last;
	trace_forward_active = false;
}

int haplotype_segment_single::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities, 
                                       vector<float>* SC, const vector<bool>* anchor_has_missing, const vector<uint32_t>* supersite_sc_offset) {
	int n_underflow_recovered = 0;
	// Set thread-local offset storage for IMPUTE_SUPERSITE_MULTIVARIATE calls
	this->supersite_sc_offset = supersite_sc_offset;
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;
	prev_abs_locus = locus_last;
	trace_forward_active = false;
	trace_backward_active = true;

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx, panel_codes_size);

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
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
		const EmitKind emit = site_view.emit_kind;
		const bool hmm_mis = (emit == EmitKind::Mis);
		const bool hmm_amb = (emit == EmitKind::Amb);
		const bool hmm_hom = (emit == EmitKind::Hom);
	// Include siblings in bookkeeping so backward cursor and missing
	// indices match the forward pass expectations.
	const bool data_mis = hmm_mis;
	const bool data_amb = hmm_amb;
		yt = (curr_abs_locus == locus_last)?0.0:M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;
		trace_ambiguous_cursor("bwd_pre", curr_abs_locus, is_sibling, 0);

		if (curr_abs_locus == locus_last) {
			if (is_anchor) {
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
					INIT_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
				} else {
					switch (emit) {
						case EmitKind::Hom:
							SS_INIT_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
							break;
						case EmitKind::Amb:
							SS_INIT_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
							break;
						case EmitKind::Mis:
							SS_INIT_MIS();
							break;
					}
				}
			} else if (is_sibling) {
				// Sibling at window end (backward): neutral init but no prev_locus advance
				INIT_MIS();
				update_prev_locus = false;
			} else {
				if (hmm_mis) {
					INIT_MIS();
				} else {
					bial_adapter.build_match_mask(site_view, n_cond_haps, curr_rel_locus + curr_rel_locus_offset, init_match_mask);
					INIT_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
				}
			}
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (is_anchor) {
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
					RUN_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
				} else {
					switch (emit) {
						case EmitKind::Hom:
							update_prev_locus = SS_RUN_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
							break;
						case EmitKind::Amb:
							update_prev_locus = SS_RUN_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
							break;
						case EmitKind::Mis:
							update_prev_locus = SS_RUN_MIS();
							break;
					}
				}
			} else if (is_sibling) {
				// Sibling within window (backward): pure bookkeeping
				RUN_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) update_prev_locus = RUN_HOM(rare_allele);
				else if (hmm_amb) RUN_AMB();
				else RUN_MIS();
			}
		} else {
			if (is_anchor) {
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
					COLLAPSE_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
				} else {
					switch (emit) {
						case EmitKind::Hom:
							SS_COLLAPSE_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
							break;
						case EmitKind::Amb:
							SS_COLLAPSE_AMB(*site_view.supersite, site_view.supersite_index, site_view.sample_class0, site_view.sample_class1);
							break;
						case EmitKind::Mis:
							SS_COLLAPSE_MIS();
							break;
					}
				}
			} else if (is_sibling) {
				// Sibling at segment boundary (backward): bookkeeping only
				COLLAPSE_SIB(site_view);
				update_prev_locus = false;
			} else {
				if (hmm_hom) COLLAPSE_HOM();
				else if (hmm_amb) COLLAPSE_AMB();
				else COLLAPSE_MIS();
			}
		}
		if (curr_segment_locus == 0) SUMK();
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

		if (curr_abs_locus == 0) SET_FIRST_TRANS(transition_probabilities);
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (data_mis) {
			bool supersite_missing_handled = false;
			if (is_anchor && anchor_has_missing && SC && site_view.supersite_index >= 0 && (*anchor_has_missing)[site_view.supersite_index]) {
				int idx = -1;
				int map_i = curr_abs_locus - locus_first;
				if (map_i >= 0 && map_i < (int)missing_index_by_locus.size()) idx = missing_index_by_locus[map_i];
				if (idx >= 0) {
					if (supersite_trace_enabled()) {
						std::fprintf(stdout, "SCTrace anchor=%d rel_missing=%d\n", curr_abs_locus, idx);
						std::fprintf(stdout, "AlphaSumMissing:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", AlphaSumMissing[idx][h]);
						std::fprintf(stdout, "\nBeta:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", prob[h + idx * HAP_NUMBER]);
						std::fprintf(stdout, "\n");
					}
					IMPUTE_SUPERSITE_MULTIVARIATE(*SC, *site_view.supersite, site_view.supersite_index, idx, supersite_sc_offset);
				} else {
					IMPUTE(missing_probabilities);
				}
				supersite_missing_handled = true;
			}
			if (!supersite_missing_handled) {
				IMPUTE(missing_probabilities);
			}
			curr_abs_missing--;
		}


		curr_segment_locus--;
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
	const bool can_retreat_amb = data_amb && has_amb_range && (curr_abs_ambiguous > ambiguous_first);
	const bool sib_retreat_amb = is_sibling && hmm_amb && has_amb_range && (curr_abs_ambiguous > ambiguous_first);
	const int expected_delta = (can_retreat_amb || sib_retreat_amb) ? -1 : 0;
		const int cursor_before_bwd = curr_abs_ambiguous;
	if (expected_delta) curr_abs_ambiguous--;
		trace_ambiguous_cursor("bwd_post", curr_abs_locus, is_sibling, expected_delta);
		if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last)) {
			if (supersite_trace_enabled()) {
				int seg_len = (curr_segment_index >= 0 && curr_segment_index < (int)G->Lengths.size()) ? G->Lengths[curr_segment_index] : -1;
				std::fprintf(stderr,
					"[ss-amb-oob][single] stage=bwd_post locus=%d before=%d after=%d expected_delta=%d range=[%d,%d] seg_idx=%d seg_len=%d seg_loc=%d is_sibling=%d\n",
					curr_abs_locus,
					cursor_before_bwd,
					curr_abs_ambiguous,
					expected_delta,
					ambiguous_first,
					ambiguous_last,
					curr_segment_index,
					seg_len,
					curr_segment_locus,
					static_cast<int>(is_sibling));
				std::fflush(stderr);
			}
			assert(false && "backward ambiguous cursor moved out of window bounds");
		}
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
		curr_segment_index = segment_first;
		curr_segment_locus = 0;
		curr_abs_locus = locus_first;
		trace_backward_active = false;
		return n_underflow_recovered;
}

void haplotype_segment_single::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	const unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	std::vector<double> cprobs(n_transitions, 0.0);

	double lane_probs[HAP_NUMBER];
	bool use_outer = supersites_enabled_flag && !AlphaSum.empty();

	if (use_outer) {
		const int rel_seg = 0; // first segment in this window
		const float* alpha_lane = AlphaSum[rel_seg].data();
		double lane_total = 0.0;
		for (int h = 0; h < HAP_NUMBER; ++h) {
			double weight = static_cast<double>(alpha_lane[h]) * static_cast<double>(probSumH[h]);
			lane_probs[h] = weight;
			lane_total += weight;
		}
		if (lane_total <= std::numeric_limits<double>::min()) {
			use_outer = false;
		} else {
			const double inv_total = 1.0 / lane_total;
			for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] *= inv_total;
		}
	}

	if (!use_outer) {
		const double lane_total = probSumT;
		const double inv_total = (lane_total > std::numeric_limits<double>::min()) ? (1.0 / lane_total) : 0.0;
		for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] = static_cast<double>(probSumH[h]) * inv_total;
	}

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
        if (underflow_trace_enabled()) {
            std::lock_guard<std::mutex> lk(g_underflow_log_mutex);
            ensure_logs_dir();
            FILE* f = std::fopen("logs/underflow.tsv", "a");
            if (f) {
                static bool header_written = false;
                if (!header_written) {
                    std::fprintf(f, "sample\tprecision\tstage\tcurr_abs_locus\tprev_abs_locus\tcm_curr\tcm_prev\tyt\tnt\tprobSumT\tsumHProbs\tAlphaSumSum_prev");
                    for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\tAlphaSum_prev[%d]", h);
                    std::fprintf(f, "\tis_anchor\tis_sibling\tamb\tmis\thom\tss_idx\tn_cond_haps\n");
                    header_written = true;
                }
                int curr = curr_abs_locus;
                int prev = prev_abs_locus;
                double cm_curr = (curr >= 0 && curr < (int)M.cm.size()) ? M.cm[curr] : NAN;
                double cm_prev = (prev >= 0 && prev < (int)M.cm.size()) ? M.cm[prev] : NAN;
                bool amb = VAR_GET_AMB(MOD2(curr), G->Variants[DIV2(curr)]);
                bool mis = VAR_GET_MIS(MOD2(curr), G->Variants[DIV2(curr)]);
                bool hom = !(amb || mis);
                genotype::SuperSiteContext log_ctx = G->getSuperSiteContext(curr);
                int ss_idx = log_ctx.ss_idx;
                bool is_anchor = log_ctx.is_anchor;
                bool is_sibling = log_ctx.is_member && !log_ctx.is_anchor;
                if (log_ctx.is_member) {
                    amb = log_ctx.has_het || log_ctx.has_sca;
                    mis = log_ctx.all_missing;
                    hom = !(amb || mis);
                }
                unsigned int rel_prev_seg = (curr_segment_index - segment_first) - 1;
                double alphaSumSum_prev = (rel_prev_seg < AlphaSumSum.size()) ? AlphaSumSum[rel_prev_seg] : NAN;
                std::fprintf(f, "%s\tfloat\t%s\t%d\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g",
                             G->name.c_str(), "SET_OTHER_TRANS", curr, prev, cm_curr, cm_prev,
                             (double)yt, (double)nt, (double)probSumT, (double)sumHProbs, alphaSumSum_prev);
                if (rel_prev_seg < AlphaSum.size()) {
                    for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\t%.8g", (double)AlphaSum[rel_prev_seg][h]);
                } else {
                    for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\tNA");
                }
                std::fprintf(f, "\t%d\t%d\t%d\t%d\t%d\t%d\t%u\n", (int)is_anchor, (int)is_sibling, (int)amb, (int)mis, (int)hom, ss_idx, n_cond_haps);
                std::fclose(f);
            }
        }
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
	curr_abs_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}

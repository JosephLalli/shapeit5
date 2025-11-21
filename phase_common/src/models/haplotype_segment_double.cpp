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
#include <mutex>
#include <cstdio>
#include <string>
#include <limits>
#include <algorithm>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

// Debug underflow tracing (gated by SHAPEIT5_DEBUG_UNDERFLOW)
namespace {
static std::mutex g_underflow_log_mutex_d;
static bool underflow_trace_enabled_d() {
    static int flag = -1;
    if (flag < 0) {
        const char* env = std::getenv("SHAPEIT5_DEBUG_UNDERFLOW");
        flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return flag == 1;
}
static void ensure_logs_dir_d() {
    struct stat st{};
    if (stat("logs", &st) != 0) {
        mkdir("logs", 0777);
    }
}

static bool supersite_trace_enabled_d() {
	static int flag = -1;
	if (flag < 0) {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
		flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
	}
	return flag == 1;
}

static bool trans_trace_enabled_d() {
	static int flag = -1;
	if (flag < 0) {
		const char* env = std::getenv("SHAPEIT5_TRANS_TRACE");
		flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
	}
	return flag == 1;
}

static const char* trans_trace_sample_d() {
	static const char* sample = std::getenv("SHAPEIT5_TRANS_TRACE_SAMPLE");
	return sample && sample[0] ? sample : nullptr;
}

} // namespace

void haplotype_segment_double::trace_ambiguous_cursor(const char* stage, int locus, bool is_sibling, int expected_delta) const {
	if (!supersite_trace_enabled_d()) return;
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
		// Stage triggered outside its matching pass: invalidate stored snapshot to avoid stale comparisons.
		if (std::strcmp(stage, "bwd_pre") == 0) {
			trace_backward_pre_valid = false;
		} else if (std::strcmp(stage, "fwd_pre") == 0) {
			trace_forward_pre_valid = false;
		}
	}

	if (delta_valid && actual_delta != expected_delta) {
		std::fprintf(stdout,
				 "[ss-amb-delta][double] stage=%s locus=%d curr_abs_ambiguous=%d expected_delta=%d actual_delta=%d seg_idx=%d is_sibling=%d\n",
				 stage,
				 locus,
				 curr_abs_ambiguous,
				 expected_delta,
				 actual_delta,
				 curr_segment_index,
				 static_cast<int>(is_sibling));
		std::fflush(stdout);
		assert(actual_delta == expected_delta && "ambiguous cursor delta mismatch");
	}

	const int lower_base = ambiguous_first;
	const int upper_base = ambiguous_last;
	const bool allow_lower_exclusive = backward_stage && trace_backward_active;
	const bool allow_upper_exclusive = forward_stage && trace_forward_active;
	const int lower_limit = allow_lower_exclusive ? (lower_base - 1) : lower_base;
	const int upper_limit = allow_upper_exclusive ? (upper_base + 1) : upper_base;
	if (curr_abs_ambiguous < lower_limit || curr_abs_ambiguous > upper_limit) {
		std::fprintf(stdout,
				 "[ss-amb-drift][double] stage=%s locus=%d curr_abs_ambiguous=%d lower_limit=%d upper_limit=%d rel=%d seg_idx=%d is_sibling=%d\n",
				 stage,
				 locus,
				 curr_abs_ambiguous,
				 lower_limit,
				 upper_limit,
				 curr_abs_ambiguous - lower_base,
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
			 "[ss-amb-cursor][double] stage=%s locus=%d curr_abs_ambiguous=%d range=[%d,%d] lower_limit=%d upper_limit=%d actual_delta=%s expected_delta=%d seg_idx=%d is_sibling=%d\n",
			 stage,
			 locus,
			 curr_abs_ambiguous,
			 lower_base,
			 upper_base,
			 lower_limit,
			 upper_limit,
			 delta_repr,
			 expected_delta,
			 curr_segment_index,
			 static_cast<int>(is_sibling));
}

haplotype_segment_double::haplotype_segment_double(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M,
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
	// Tests may skip explicit supersite setup; attach context and snapshot immutable
	// supersite base classes here to keep emissions well-defined.
	if (G && supersites_enabled_flag) {
		if (!G->super_sites || !G->locus_to_super_idx || !G->super_site_var_index) {
			G->setSuperSiteContext(super_sites, locus_to_super_idx, super_site_var_index, nullptr, nullptr, nullptr);
		}
		if (G->supersite_class_pairs_base.empty()) {
			G->snapshotSupersiteBaseClasses(*super_sites, *super_site_var_index);
		}
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
    AlphaLaneSum = vector<LaneMarginal>(segment_last - segment_first + 1, LaneMarginal{});
	AlphaSumSum = aligned_vector32 < double > (segment_last - segment_first + 1, 0.0);
	if (n_missing > 0) {
		AlphaMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER * n_cond_haps, 0.0));
		AlphaSumMissing = vector < aligned_vector32 < double > > (n_missing, aligned_vector32 < double > (HAP_NUMBER, 0.0));
	}
	init_match_mask.resize(static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER);
	//Cache efficient data transfer for conditioning haplotypes
	curr_rel_locus_offset = Hhap.subset(H, idxH, locus_first, locus_last);
	Hvar.allocateFast(Hhap.n_cols, Hhap.n_rows);
	Hhap.transpose(Hvar);

	// allocate scratch for supersites if enabled
	if (super_sites && locus_to_super_idx && panel_codes && super_site_var_index) {
		ss_cond_codes = aligned_vector32<uint8_t>(n_cond_haps, 0);
		ss_emissions = aligned_vector32<double>(n_cond_haps, 1.0);
		ss_emissions_h1 = aligned_vector32<double>(n_cond_haps, 1.0);
		ss_cached.resize(super_sites->size(), false);  // Phase 3: initialize cache flags
	}

	trace_forward_pre_cursor = 0;
	trace_forward_pre_locus = -1;
	trace_forward_pre_valid = false;
	trace_backward_pre_cursor = 0;
	trace_backward_pre_locus = -1;
	trace_backward_pre_valid = false;
	trace_forward_active = false;
	trace_backward_active = false;

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
	AlphaLaneSum.clear();
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
	trace_forward_active = true;
	trace_backward_active = false;
	// Diagnostics: track ambiguous-site bookkeeping counts to verify cursor correctness
	int diag_expected_amb_sites = 0; // number of data_amb sites (excluding siblings) encountered
	int diag_advanced_amb = 0;      // number of times we actually advanced the ambiguous cursor

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx, panel_codes_size);

    const char* tr_d = std::getenv("SHAPEIT5_TEST_TRACE");
    if (tr_d && tr_d[0] != '\0' && tr_d[0] != '0') {
        std::fprintf(stdout,
            "D.FWD_start loci=[%d,%d] seg=[%d,%d] n_cond_haps=%u prob_size=%zu probSumH_size=%zu ptr_prob=%p ptr_probSumH=%p\n",
            locus_first, locus_last, segment_first, segment_last, n_cond_haps,
            prob.size(), probSumH.size(), (void*)prob.data(), (void*)probSumH.data());
    }

	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool update_prev_locus = true;
		char rare_allele = M.rare_allele[curr_abs_locus];

		SiteView site_view{};
		bool has_supersite = supersites_enabled && supersite_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		if (!has_supersite) {
			bial_adapter.build_view(curr_abs_locus, curr_abs_ambiguous, site_view);
		}
		if (tr_d && tr_d[0] != '\0' && tr_d[0] != '0') {
			const bool dbg_is_anchor = (site_view.kind == SiteKind::SuperAnchor);
			const bool dbg_is_sibling = (site_view.kind == SiteKind::SuperSibling);
			std::fprintf(stdout,
				"D.site locus=%d kind=%d emit=%d is_anchor=%d is_sibling=%d ss_idx=%d anchor_gid=%d sample_cls0=%u sample_cls1=%u amb_mask=0x%02x curr_abs_amb=%d range=[%d,%d] seg_idx=%d seg_loc=%d\n",
				curr_abs_locus,
				static_cast<int>(site_view.kind),
				static_cast<int>(site_view.emit_kind),
				dbg_is_anchor ? 1 : 0,
				dbg_is_sibling ? 1 : 0,
				site_view.supersite_index,
				site_view.supersite ? static_cast<int>(site_view.supersite->global_site_id) : -1,
				site_view.sample_class0,
				site_view.sample_class1,
				site_view.amb_mask,
				curr_abs_ambiguous,
				ambiguous_first,
				ambiguous_last,
				curr_segment_index,
				curr_segment_locus);
		}
		const bool is_anchor = (site_view.kind == SiteKind::SuperAnchor);
		const bool is_sibling = (site_view.kind == SiteKind::SuperSibling);
		const EmitKind emit = site_view.emit_kind;
		const bool hmm_mis = (emit == EmitKind::Mis);
		const bool hmm_amb = (emit == EmitKind::Amb);
		const bool hmm_hom = (emit == EmitKind::Hom);
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb && !is_sibling;
		trace_ambiguous_cursor("bwd_pre", curr_abs_locus, is_sibling, 0);
		trace_ambiguous_cursor("fwd_pre", curr_abs_locus, is_sibling, 0);

		yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0 - yt;

            if (curr_rel_locus == 0) {
                if (is_anchor) {
                    if (M.ss_anchor_split_emissions) {
                        const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
                        if (tr && tr[0] != '\0' && tr[0] != '0') {
                            std::fprintf(stdout, "forward: INIT via mask at anchor locus=%d (double)\n", curr_abs_locus);
                        }
                        supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                        INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
                    } else {
					switch (emit) {
						case EmitKind::Hom: SS_INIT_HOM(); break;
						case EmitKind::Amb: SS_INIT_AMB(); break;
						case EmitKind::Mis: SS_INIT_MIS(); break;
					}
				}
			} else if (is_sibling) {
				// Sibling at window start: initialize neutrally but do not advance prev_abs_locus
				INIT_MIS();
				update_prev_locus = false;
			} else {
				if (hmm_mis) {
					INIT_MIS();
				} else {
					bial_adapter.build_match_mask(site_view, n_cond_haps, curr_rel_locus + curr_rel_locus_offset, init_match_mask);
					INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
				}
			}
            } else if (curr_segment_locus != 0) {
                if (is_anchor) {
                    if (M.ss_anchor_split_emissions) {
                        const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
                        if (tr && tr[0] != '\0' && tr[0] != '0') {
                            std::fprintf(stdout, "forward: RUN via mask at anchor locus=%d (double)\n", curr_abs_locus);
                        }
                        supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                        RUN_FROM_MASK(init_match_mask, M.ed/M.ee);
                    } else {
					switch (emit) {
						case EmitKind::Hom:
							update_prev_locus = SS_RUN_HOM();
							break;
						case EmitKind::Amb:
							SS_RUN_AMB();
							break;
						case EmitKind::Mis:
							SS_RUN_MIS();
							break;
					}
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
                    if (M.ss_anchor_split_emissions) {
                        const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
                        if (tr && tr[0] != '\0' && tr[0] != '0') {
                            std::fprintf(stdout, "forward: COLLAPSE via mask at anchor locus=%d (double)\n", curr_abs_locus);
                        }
                        supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
                        COLLAPSE_FROM_MASK(init_match_mask, M.ed/M.ee);
                    } else {
					switch (emit) {
						case EmitKind::Hom: SS_COLLAPSE_HOM(); break;
						case EmitKind::Amb: SS_COLLAPSE_AMB(); break;
						case EmitKind::Mis: SS_COLLAPSE_MIS(); break;
					}
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

        if (curr_segment_locus == (G->Lengths[curr_segment_index] - 1)) SUMK();

        // Forward-side diagnostic: detect non-finite or underflowed probSumT early
        if (underflow_trace_enabled_d()) {
            if (!std::isfinite(probSumT) || probSumT <= std::numeric_limits<double>::min()) {
                std::lock_guard<std::mutex> lk(g_underflow_log_mutex_d);
                ensure_logs_dir_d();
                FILE* f = std::fopen("logs/underflow.tsv", "a");
                if (f) {
                    static bool header_written_fwd = false;
                    if (!header_written_fwd) {
                        std::fprintf(f, "sample\tprecision\tstage\tcurr_abs_locus\tprev_abs_locus\tcm_curr\tcm_prev\tyt\tnt\tprobSumT\tAlphaSumSum_curr");
                        for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\tprobSumH[%d]", h);
                        std::fprintf(f, "\tamb\tmis\thom\tss_idx\tn_cond_haps\n");
                        header_written_fwd = true;
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
					if (log_ctx.is_member) {
						amb = log_ctx.has_het || log_ctx.has_sca;
						mis = log_ctx.all_missing;
						hom = !(amb || mis);
					}
					double alphaSumSum_curr = probSumT; // at this locus, same scalar
					std::fprintf(f, "%s\tdouble\tFORWARD_NAN\t%d\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.17g\t%.17g",
						         G->name.c_str(), curr, prev, cm_curr, cm_prev,
						         yt, nt, probSumT, alphaSumSum_curr);
                    for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\t%.17g", probSumH[h]);
                    std::fprintf(f, "\t%d\t%d\t%d\t%d\t%u\n", (int)amb, (int)mis, (int)hom, ss_idx, n_cond_haps);
                    std::fclose(f);
                }
            }
        }
	if (curr_segment_locus == G->Lengths[curr_segment_index] - 1) {
			const int rel_seg = curr_segment_index - segment_first;
			Alpha[rel_seg] = prob;
			AlphaSum[rel_seg] = probSumH;
			__m256d lane_lo = _mm256_load_pd(&probSumH[0]);
			__m256d lane_hi = _mm256_load_pd(&probSumH[4]);
			_mm256_store_pd(&AlphaLaneSum[rel_seg].lane[0], lane_lo);
			_mm256_store_pd(&AlphaLaneSum[rel_seg].lane[4], lane_hi);
			AlphaSumSum[rel_seg] = probSumT;
			AlphaLocus[rel_seg] = prev_abs_locus;
	}
        if (data_mis) {
            if (curr_rel_missing < 0 || curr_rel_missing >= (int)AlphaMissing.size()) {
                const char* tr_d = std::getenv("SHAPEIT5_TEST_TRACE");
                if (tr_d && tr_d[0] != '\0' && tr_d[0] != '0') {
                    std::fprintf(stdout, "D.FWD store AlphaMissing OOB: rel_missing=%d size=%zu locus=%d\n",
                                 curr_rel_missing, AlphaMissing.size(), curr_abs_locus);
                }
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

	curr_segment_locus ++;
	const bool has_amb_range = (ambiguous_first <= ambiguous_last);
	const int cursor_before = curr_abs_ambiguous;
	const bool can_advance_amb = data_amb && has_amb_range && (curr_abs_ambiguous <= ambiguous_last);
	const int expected_delta = can_advance_amb ? 1 : 0;

	// Diagnostic accounting: count data_amb sites and actual advances
	if (expected_delta) diag_expected_amb_sites++;

	if (tr_d && tr_d[0] != '\0' && tr_d[0] != '0') {
		std::fprintf(stdout,
			"D.delta locus=%d cursor_before=%d expected_delta=%d data_amb=%d data_mis=%d is_sibling=%d has_range=%d\n",
			curr_abs_locus,
			cursor_before,
			expected_delta,
			data_amb ? 1 : 0,
			data_mis ? 1 : 0,
			is_sibling ? 1 : 0,
			has_amb_range ? 1 : 0);
	}
	if (can_advance_amb) {
		curr_abs_ambiguous++;
		diag_advanced_amb++;
	}
	trace_ambiguous_cursor("fwd_post", curr_abs_locus, is_sibling, expected_delta);
	if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last + 1)) {
		if (supersite_trace_enabled_d()) {
			int seg_len = (curr_segment_index >= 0 && curr_segment_index < (int)G->Lengths.size()) ? G->Lengths[curr_segment_index] : -1;
			std::fprintf(stderr,
				"[ss-amb-oob][double] stage=fwd_post locus=%d before=%d after=%d expected_delta=%d range=[%d,%d] seg_idx=%d seg_len=%d seg_loc=%d is_sibling=%d\n",
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

		// End-of-window diagnostic: verify ambiguous-site bookkeeping parity
		if (supersite_trace_enabled_d()) {
			const bool has_amb_range = (ambiguous_first <= ambiguous_last);
			int slots = has_amb_range ? (ambiguous_last - ambiguous_first + 1) : 0;
			if (tr_d && tr_d[0] != '\0' && tr_d[0] != '0') {
				std::fprintf(stdout, "D.diag-check locus=%d diag_expected=%d diag_advanced=%d slots=%d\n",
					curr_abs_locus, diag_expected_amb_sites, diag_advanced_amb, slots);
			}
			if (diag_expected_amb_sites > slots) {
				std::fprintf(stderr, "[ss-amb-diag][double] expected_amb_sites=%d slots=%d\n", diag_expected_amb_sites, slots);
				assert(false && "more ambiguous data sites than available slots");
			}
			if (diag_expected_amb_sites != diag_advanced_amb) {
				std::fprintf(stderr, "[ss-amb-diag][double] expected_amb_sites=%d advanced_amb=%d slots=%d\n",
					 diag_expected_amb_sites, diag_advanced_amb, slots);
				assert(false && "ambiguous cursor advanced count mismatch");
			}
			else {
				std::fprintf(stdout, "[ss-amb-diag-PASS][double] expected_amb_sites=%d advanced_amb=%d slots=%d\n",
					 diag_expected_amb_sites, diag_advanced_amb, slots);
				std::fflush(stdout);
			}
		}
	}
	trace_forward_active = false;
}

int haplotype_segment_double::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities,
                                       vector < float > * SC, const vector < bool > * anchor_has_missing, const vector<uint32_t>* supersite_sc_offset) {
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

	// Flag: backward pass always starts with initialization; if locus_last is a sibling, defer until first non-sibling
	bool need_init = true;

	// DIAGNOSTIC: Track actual transition storage count vs expected
	unsigned int debug_transitions_written = 0;
	const bool debug_track_transitions = !transition_probabilities.empty();
	const bool trans_trace = trans_trace_enabled_d();
	const char* trans_trace_sample = trans_trace_sample_d();

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

		// Deferred initialization: if we started on a sibling, keep doing INIT_SIB until we hit a non-sibling
		if (need_init && is_sibling) {
			INIT_SIB(site_view);
			update_prev_locus = false;
			trace_ambiguous_cursor("bwd_post", curr_abs_locus, is_sibling, 0);
			prev_abs_locus = update_prev_locus ? curr_abs_locus : prev_abs_locus;
			// CRITICAL: Check for segment boundary BEFORE decrementing, since we'll skip the normal check with continue
			if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
				if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
					unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
					unsigned int prev_dipcount = (curr_segment_index > 0) ? G->countDiplotypes(G->Diplotypes[curr_segment_index-1]) : 0;
					unsigned int n_trans = curr_dipcount * prev_dipcount;
					std::fprintf(stdout,
					             "[TRANS_TRACE][double][skip-sib] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d stored_so_far=%d/%u buf_size=%zu\n",
					             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, debug_transitions_written, G->n_transitions, transition_probabilities.size());
					std::fflush(stdout);
				}
				unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
				unsigned int prev_dipcount = (curr_segment_index > 0) ? G->countDiplotypes(G->Diplotypes[curr_segment_index-1]) : 0;
				unsigned int n_trans = curr_dipcount * prev_dipcount;
				int ret = SET_OTHER_TRANS(transition_probabilities);
				if (debug_track_transitions) debug_transitions_written += n_trans;
				if (ret < 0) return ret;
				else n_underflow_recovered += ret;
			}
			// CRITICAL: Decrement curr_segment_locus even for skipped siblings
			// Otherwise segments ending on siblings will never hit curr_segment_locus==0
			// and SET_OTHER_TRANS() will never be called, causing vector size mismatch
			curr_segment_locus--;
			if (curr_segment_locus < 0 && curr_segment_index > 0) {
				curr_segment_index--;
				curr_segment_locus = G->Lengths[curr_segment_index] - 1;
			}
			continue;
		}
		const EmitKind emit = site_view.emit_kind;
		const bool hmm_mis = (emit == EmitKind::Mis);
		const bool hmm_amb = (emit == EmitKind::Amb);
		const bool hmm_hom = (emit == EmitKind::Hom);
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb && !is_sibling;
		trace_ambiguous_cursor("bwd_pre", curr_abs_locus, is_sibling, 0);

		yt = (curr_abs_locus == locus_last)?0.0:M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		// Deferred initialization: first non-sibling after starting on a sibling - use INIT
		if (need_init && is_anchor) {
			need_init = false;
			if (M.ss_anchor_split_emissions) {
				supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
				INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
			} else {
				switch (emit) {
					case EmitKind::Hom: SS_INIT_HOM(); break;
					case EmitKind::Amb: SS_INIT_AMB(); break;
					case EmitKind::Mis: SS_INIT_MIS(); break;
				}
			}
		} else if (need_init && !is_sibling) {
			// Deferred initialization: first biallelic after starting on a sibling - use INIT
			need_init = false;
			if (hmm_mis) {
				INIT_MIS();
			} else {
				bial_adapter.build_match_mask(site_view, n_cond_haps, curr_rel_locus + curr_rel_locus_offset, init_match_mask);
				INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
			}
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (is_anchor) {
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
					RUN_FROM_MASK(init_match_mask, M.ed/M.ee);
				} else {
					switch (emit) {
						case EmitKind::Hom:
							update_prev_locus = SS_RUN_HOM();
							break;
						case EmitKind::Amb:
							SS_RUN_AMB();
							break;
						case EmitKind::Mis:
							SS_RUN_MIS();
							break;
					}
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
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/true, init_match_mask);
					COLLAPSE_FROM_MASK(init_match_mask, M.ed/M.ee);
				} else {
					switch (emit) {
						case EmitKind::Hom: SS_COLLAPSE_HOM(); break;
						case EmitKind::Amb: SS_COLLAPSE_AMB(); break;
						case EmitKind::Mis: SS_COLLAPSE_MIS(); break;
					}
				}
			} else if (is_sibling) {
				// Sibling at segment boundary (backward): no-op
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

		if (curr_abs_locus == 0) {
			unsigned int n_trans = G->countDiplotypes(G->Diplotypes[0]);
			if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
				std::fprintf(stdout,
				             "[TRANS_TRACE][double][first] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d stored_so_far=%d/%u buf_size=%zu\n",
				             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, debug_transitions_written, G->n_transitions, transition_probabilities.size());
				std::fflush(stdout);
			}
			SET_FIRST_TRANS(transition_probabilities);
			if (debug_track_transitions) debug_transitions_written += n_trans;
		}
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
				unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
				unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
				unsigned int n_trans = curr_dipcount * prev_dipcount;
				std::fprintf(stdout,
				             "[TRANS_TRACE][double][other] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d stored_so_far=%d/%u buf_size=%zu\n",
				             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, debug_transitions_written, G->n_transitions, transition_probabilities.size());
				std::fflush(stdout);
			}
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (debug_track_transitions) {
				unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
				unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
				debug_transitions_written += curr_dipcount * prev_dipcount;
			}
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
					if (supersite_trace_enabled_d()) {
						std::fprintf(stdout, "SCTrace(double) anchor=%d rel_missing=%d\n", curr_abs_locus, rel_idx);
						std::fprintf(stdout, "AlphaSumMissing:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.9f", AlphaSumMissing[rel_idx][h]);
						std::fprintf(stdout, "\nBeta:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.9f", prob[h + rel_idx * HAP_NUMBER]);
						std::fprintf(stdout, "\n");
					}
					IMPUTE_SUPERSITE_MULTIVARIATE(*SC, *site_view.supersite, site_view.supersite_index, rel_idx, supersite_sc_offset);
					supersite_handled = true;
				}
			}
			if (!supersite_handled) {
				IMPUTE(missing_probabilities);
			}
			curr_abs_missing--;
		}

		curr_segment_locus--;
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
		const bool can_retreat_amb = data_amb && has_amb_range && (curr_abs_ambiguous >= ambiguous_first);
		const int expected_delta = can_retreat_amb ? -1 : 0;
		const int cursor_before_bwd = curr_abs_ambiguous;
		if (can_retreat_amb) curr_abs_ambiguous--;
		trace_ambiguous_cursor("bwd_post", curr_abs_locus, is_sibling, expected_delta);
		const int lower_bound = ambiguous_first - 1;
		if (has_amb_range && (curr_abs_ambiguous < lower_bound || curr_abs_ambiguous > ambiguous_last)) {
			if (supersite_trace_enabled_d()) {
				int seg_len = (curr_segment_index >= 0 && curr_segment_index < (int)G->Lengths.size()) ? G->Lengths[curr_segment_index] : -1;
				std::fprintf(stderr,
					"[ss-amb-oob][double] stage=bwd_post locus=%d before=%d after=%d expected_delta=%d range=[%d,%d] seg_idx=%d seg_len=%d seg_loc=%d is_sibling=%d\n",
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
	trace_backward_active = false;

	// DIAGNOSTIC: Check if transitions_stored matches expected
	if (debug_track_transitions && debug_transitions_written != G->n_transitions) {
		std::fprintf(stderr, "[TRANS_MISMATCH_ERROR] Sample=%s: Backward pass stored %u transitions but expected %u\n",
		             G->name.c_str(), debug_transitions_written, G->n_transitions);
		std::fprintf(stderr, "  Vector size=%zu, locus_range=[%d,%d], n_segments=%u\n",
		             transition_probabilities.size(), locus_first, locus_last, G->n_segments);
		std::fprintf(stderr, "  This will cause out-of-bounds access in sampleBackward!\n");
	}
	return n_underflow_recovered;
}

void haplotype_segment_double::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	const unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	std::vector<double> cprobs(n_transitions, 0.0);

	// Guard: caller must provision at least n_transitions slots.
	if (transition_probabilities.size() < n_transitions) {
		if (supersite_trace_enabled_d()) {
			std::fprintf(stdout,
			             "SET_FIRST_TRANS double: skip write (insufficient buffer) n_transitions=%u buf_size=%zu\n",
			             n_transitions, transition_probabilities.size());
		}
		curr_abs_transition -= (n_transitions - 1);
		curr_abs_transition --;
		return;
	}

	double lane_probs[HAP_NUMBER];
	bool use_outer = supersites_enabled_flag && !AlphaSum.empty();

	if (use_outer) {
		const int rel_seg = 0;
		const double* alpha_lane = AlphaSum[rel_seg].data();
		double lane_total = 0.0;
		for (int h = 0; h < HAP_NUMBER; ++h) {
			double weight = alpha_lane[h] * probSumH[h];
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
		for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] = probSumH[h] * inv_total;
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

int haplotype_segment_double::SET_OTHER_TRANS(vector < double > & transition_probabilities) {
    int underflow_recovered = 0;
    if (TRANS_HAP()) {
        if (underflow_trace_enabled_d()) {
            std::lock_guard<std::mutex> lk(g_underflow_log_mutex_d);
            ensure_logs_dir_d();
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
                std::fprintf(f, "%s\tdouble\t%s\t%d\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g",
                             G->name.c_str(), "SET_OTHER_TRANS", curr, prev, cm_curr, cm_prev,
                             yt, nt, probSumT, sumHProbs, alphaSumSum_prev);
                if (rel_prev_seg < AlphaSum.size()) {
                    for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(f, "\t%.8g", AlphaSum[rel_prev_seg][h]);
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
        bool under = TRANS_DIP_ADD();
        if (under) {
            if (trans_parity_trace_enabled_d()) {
                std::fprintf(stderr, "[TRANS_DIP_ADD debug][double] locus=%d seg=%d curr_abs_transition=%d sumHProbs=%g sumDProbs=%g\n",
                             curr_abs_locus, curr_segment_index, curr_abs_transition, sumHProbs, sumDProbs);
                // Dump first few HProbs rows to compare with single path
                for (int h = 0; h < std::min(4, HAP_NUMBER); ++h) {
                    std::fprintf(stderr, "  HProbs[%d]:", h);
                    for (int hh = 0; hh < HAP_NUMBER; ++hh) std::fprintf(stderr, " %.6g", HProbs[h*HAP_NUMBER + hh]);
                    std::fprintf(stderr, "\n");
                }
            }
            return -2;
        } else underflow_recovered = 1;
    }
	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;

	// Bounds check before writing into transition_probabilities.
	const long long start_idx = static_cast<long long>(curr_abs_transition) - static_cast<long long>(n_transitions) + 1;
	const size_t buf_size = transition_probabilities.size();
	if (transition_probabilities.empty() ||
	    start_idx < 0 ||
	    (static_cast<size_t>(start_idx) + n_transitions) > buf_size) {
		if (supersite_trace_enabled_d()) {
			std::fprintf(stdout, "SET_OTHER_TRANS double: skip write (insufficient buffer) start_idx=%lld n_transitions=%u buf_size=%zu\n",
			             start_idx, n_transitions, buf_size);
		}
		curr_abs_transition -= (n_transitions - 1);
		curr_abs_transition --;
		return underflow_recovered;
	}

	curr_abs_transition -= (n_transitions - 1);
	for (int t = 0 ; t < n_transitions ; t ++) transition_probabilities[curr_abs_transition + t] = DProbs[t] * scaleDip;
	curr_abs_transition --;
	return underflow_recovered;
}

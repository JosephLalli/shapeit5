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
#include <models/supersite_trace_utils.h>
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

static bool trans_trace_enabled_s() {
	static int flag = -1;
	if (flag < 0) {
		const char* env = std::getenv("SHAPEIT5_TRANS_TRACE");
		flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
	}
	return flag == 1;
}

static const char* trans_trace_sample_s() {
	static const char* sample = std::getenv("SHAPEIT5_TRANS_TRACE_SAMPLE");
	return sample && sample[0] ? sample : nullptr;
}

} // namespace

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
	// Only log ambiguous cursor for a sibling locus using a neutral stage label
	// to avoid interfering with the standard fwd_pre/fwd_post or bwd_pre/bwd_post
	// diagnostics that run once per locus outside of the sibling path.
	trace_ambiguous_cursor("sib", curr_abs_locus, true, 0);
	// Optionally, add asserts here to ensure ambiguous/missing cursor is in bounds
	if (ambiguous_first <= ambiguous_last) {
		const bool backward_stage = trace_backward_active && !trace_forward_active;
		const int lower_inclusive = backward_stage ? (ambiguous_first - 1) : ambiguous_first;
		const int upper_exclusive = ambiguous_last + 1;
		assert(curr_abs_ambiguous >= lower_inclusive && curr_abs_ambiguous <= upper_exclusive);
	}
	if (missing_first <= missing_last) {
		const bool backward_stage = trace_backward_active && !trace_forward_active;
		const int lower_inclusive = backward_stage ? (missing_first - 1) : missing_first;
		const int upper_exclusive = missing_last + 1;
		assert(curr_abs_missing >= lower_inclusive && curr_abs_missing <= upper_exclusive);
	}
	if (supersite_trace_enabled()) {
		supersite_trace_log("[SupersiteSibling] stage=%s locus=%d kind=%d curr_abs_amb=%d curr_abs_mis=%d\n",
		                    trace_forward_active ? "FWD" : "BWD",
		                    curr_abs_locus,
		                    static_cast<int>(site_view.kind),
		                    curr_abs_ambiguous,
		                    curr_abs_missing);
	}
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

	const int lower_base = ambiguous_first;
	const int upper_base = ambiguous_last;
	const bool allow_lower_exclusive = backward_stage && trace_backward_active;
	const bool allow_upper_exclusive = forward_stage && trace_forward_active;
	const int lower_limit = allow_lower_exclusive ? (lower_base - 1) : lower_base;
	const int upper_limit = allow_upper_exclusive ? (upper_base + 1) : upper_base;
	if (curr_abs_ambiguous < lower_limit || curr_abs_ambiguous > upper_limit) {
		std::fprintf(stdout,
				 "[ss-amb-drift][single] stage=%s locus=%d curr_abs_ambiguous=%d lower_limit=%d upper_limit=%d rel=%d seg_idx=%d is_sibling=%d\n",
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
			 "[ss-amb-cursor][single] stage=%s locus=%d curr_abs_ambiguous=%d range=[%d,%d] lower_limit=%d upper_limit=%d actual_delta=%s expected_delta=%d seg_idx=%d is_sibling=%d\n",
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

void haplotype_segment_single::trace_log_forward_state(int locus,
                                                       int prev_before,
                                                       int prev_after,
                                                       double yt_val,
                                                       double nt_val,
                                                       bool update_prev,
                                                       bool is_anchor,
                                                       bool is_sibling,
                                                       bool hmm_amb,
                                                       bool hmm_mis,
                                                       bool hmm_hom) {
	if (!supersite_trace_enabled()) return;
	if (!trace_forward_table_header_emitted) {
		std::fprintf(stdout,
		             "# Forward trace - sample=%s loci=[%d,%d] segments=[%d,%d] n_cond_haps=%u\n",
		             G ? G->name.c_str() : "<null>",
		             locus_first,
		             locus_last,
		             segment_first,
		             segment_last,
		             n_cond_haps);
		std::fprintf(stdout,
		             "locus\tprev_before\tprev_after\tyt\tnt\tupdate_prev\tis_anchor\tis_sibling\tamb\tmis\thom");
		for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, "\tprobSumH[%d]", h);
		std::fprintf(stdout, "\tprobSumT\n");
		trace_forward_table_header_emitted = true;
	}
	std::fprintf(stdout,
	             "%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\t%d",
	             locus,
	             prev_before,
	             prev_after,
	             yt_val,
	             nt_val,
	             update_prev ? 1 : 0,
	             is_anchor ? 1 : 0,
	             is_sibling ? 1 : 0,
	             hmm_amb ? 1 : 0,
	             hmm_mis ? 1 : 0,
	             hmm_hom ? 1 : 0);
	for (int h = 0; h < HAP_NUMBER; ++h) {
		std::fprintf(stdout, "\t%.6f", static_cast<double>(probSumH[h]));
	}
	std::fprintf(stdout, "\t%.6f\n", static_cast<double>(probSumT));
}

void haplotype_segment_single::trace_log_backward_state(int locus,
                                                        int prev_before,
                                                        int prev_after,
                                                        double yt_val,
                                                        double nt_val,
                                                        bool update_prev,
                                                        bool is_anchor,
                                                        bool is_sibling,
                                                        bool hmm_amb,
                                                        bool hmm_mis,
                                                        bool hmm_hom) {
	if (!supersite_trace_enabled()) return;
	if (!trace_backward_table_header_emitted) {
		std::fprintf(stdout,
		             "# Backward trace - sample=%s loci=[%d,%d] segments=[%d,%d] n_cond_haps=%u\n",
		             G ? G->name.c_str() : "<null>",
		             locus_first,
		             locus_last,
		             segment_first,
		             segment_last,
		             n_cond_haps);
		std::fprintf(stdout,
		             "locus\tprev_before\tprev_after\tyt\tnt\tupdate_prev\tis_anchor\tis_sibling\tamb\tmis\thom");
		for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, "\tprobSumH[%d]", h);
		std::fprintf(stdout, "\tprobSumT\n");
		trace_backward_table_header_emitted = true;
	}
	std::fprintf(stdout,
	             "%d\t%d\t%d\t%.6f\t%.6f\t%d\t%d\t%d\t%d\t%d\t%d",
	             locus,
	             prev_before,
	             prev_after,
	             yt_val,
	             nt_val,
	             update_prev ? 1 : 0,
	             is_anchor ? 1 : 0,
	             is_sibling ? 1 : 0,
	             hmm_amb ? 1 : 0,
	             hmm_mis ? 1 : 0,
	             hmm_hom ? 1 : 0);
	for (int h = 0; h < HAP_NUMBER; ++h) {
		std::fprintf(stdout, "\t%.6f", static_cast<double>(probSumH[h]));
	}
	std::fprintf(stdout, "\t%.6f\n", static_cast<double>(probSumT));
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
	// Tests may forget to attach supersite context; do it defensively so emissions
	// always see immutable base classes (c0/c1) even without explicit setup.
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

        // Populate static panel code matrix (like Hvar for biallelic)
        // Take a snapshot of current cond_idx to avoid dynamic dependency
        if (super_sites && cond_idx) {
            ss_panel_matrix.resize(super_sites->size());

            // MATRIX_POPULATE_TRACE: Log when and how matrix is populated
            const char* matrix_trace = std::getenv("SHAPEIT5_MATRIX_POPULATE_TRACE");
            if (matrix_trace && matrix_trace[0] != '\0' && matrix_trace[0] != '0') {
                static int constructor_count = 0;
                constructor_count++;
                std::fprintf(stderr, "\n[MATRIX_POPULATE] ===== Constructor #%d called =====\n", constructor_count);
                std::fprintf(stderr, "[MATRIX_POPULATE] Segment=[%d,%d] n_cond_haps=%u n_supersites=%zu\n",
                             segment_first, segment_last, n_cond_haps, super_sites->size());
                std::fprintf(stderr, "[MATRIX_POPULATE] cond_idx pointer=%p\n", (void*)cond_idx);
                std::fprintf(stderr, "[MATRIX_POPULATE] cond_idx values (first 16): ");
                for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
                    std::fprintf(stderr, "%u ", (*cond_idx)[k]);
                }
                std::fprintf(stderr, "\n");
            }

            for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
                const SuperSite& ss = (*super_sites)[ss_idx];
                ss_panel_matrix[ss_idx].resize(n_cond_haps);

                // Unpack panel codes for all conditioning haplotypes at this supersite
                for (unsigned int k = 0; k < n_cond_haps; ++k) {
                    unsigned int gh = (*cond_idx)[k];  // Snapshot: read cond_idx ONCE here
                    ss_panel_matrix[ss_idx][k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
                }

                // Log panel codes for ss_idx 3 (locus 6, the diverging one)
                if (matrix_trace && matrix_trace[0] != '\0' && matrix_trace[0] != '0' && ss_idx == 3) {
                    std::fprintf(stderr, "[MATRIX_POPULATE] ss_idx=3 (locus 6) panel_codes: ");
                    for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
                        std::fprintf(stderr, "%u ", (unsigned)ss_panel_matrix[ss_idx][k]);
                    }
                    std::fprintf(stderr, "\n");
                }
            }
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
	trace_forward_table_header_emitted = false;
	trace_backward_table_header_emitted = false;

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
	trace_forward_table_header_emitted = false;

	// Diagnostics: track ambiguous-site bookkeeping counts to verify cursor correctness
	int diag_expected_amb_sites = 0; // number of data_amb sites (excluding siblings) encountered
	int diag_advanced_amb = 0;      // number of times we actually advanced the ambiguous cursor

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
    // Reset per-window trace budget for anchor match summaries
    trace_anchor_match_logs_remaining = 2;

    for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
        if (supersite_trace_enabled()) std::fprintf(stderr, "FWD3 loop enter abs=%d rel_off=%d\n", curr_abs_locus, curr_rel_locus_offset);
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		const int prev_before = prev_abs_locus;
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
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
                    supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
		const int prev_after = prev_abs_locus;
		trace_log_forward_state(curr_abs_locus,
		                        prev_before,
		                        prev_after,
		                        static_cast<double>(yt),
		                        static_cast<double>(nt),
		                        update_prev_locus,
		                        is_anchor,
		                        is_sibling,
		                        hmm_amb,
		                        hmm_mis,
		                        hmm_hom);

		if (curr_segment_locus == (G->Lengths_bio[curr_segment_index] - 1)) SUMK();
		if (curr_segment_locus == G->Lengths_bio[curr_segment_index] - 1) {
			const int rel_seg = curr_segment_index - segment_first;
			// Diagnostic guard: print sizes/indices when trace enabled to catch OOB
			if (supersite_trace_enabled()) {
				std::fprintf(stdout,
					"DBG Alpha write: rel_seg=%d curr_segment_index=%d segment_first=%d Alpha.size=%zu AlphaSum.size=%zu AlphaLaneSum.size=%zu AlphaSumSum.size=%zu G.Lengths[curr_segment_index]=%d prob.size=%zu n_cond_haps=%u\n",
					rel_seg, curr_segment_index, segment_first,
					Alpha.size(), AlphaSum.size(), AlphaLaneSum.size(), AlphaSumSum.size(),
					(curr_segment_index >= 0 && curr_segment_index < (int)G->Lengths.size()) ? G->Lengths[curr_segment_index] : -1,
					prob.size(), (unsigned)n_cond_haps);
			}
			// Safety assert to capture the exact failure during testing
			assert(rel_seg >= 0 && rel_seg < static_cast<int>(Alpha.size()) && "rel_seg out of range when writing Alpha");
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
						// Record the rel-missing index for backward SC (supersites) or IMPUTE (biallelic)
						int idx = curr_abs_locus - locus_first;
						if (idx >= 0 && idx < (int)missing_index_by_locus.size()) missing_index_by_locus[idx] = curr_rel_missing;
						curr_abs_missing ++;
					}
		// Only increment segment locus for non-sibling loci (siblings don't count toward segment length)
		if (!is_sibling) {
			curr_segment_locus ++;
		}
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
		const int cursor_before = curr_abs_ambiguous;
		const bool can_advance_amb = data_amb && has_amb_range && (curr_abs_ambiguous <= ambiguous_last);
		// CRITICAL FIX: Siblings do not advance the ambiguous cursor
		// Sibling variants do not have entries in the Ambiguous array (only anchors do).
		// If siblings advanced curr_abs_ambiguous, it would cause out-of-bounds access and indexing errors.
		// This fix prevents cursor drift and maintains proper alignment with the Ambiguous data structure.
		const bool sib_advance_amb = false; //is_sibling && hmm_amb && has_amb_range && (curr_abs_ambiguous < ambiguous_last);
		const int expected_delta = (can_advance_amb) ? 1 : 0;

		// Diagnostic accounting: count data_amb sites and actual advances
		if (expected_delta) {
			diag_expected_amb_sites++;
			if (supersite_trace_enabled()) {
				const char* amb_reason = is_anchor ? "super_anchor" : "biallelic";
				std::fprintf(stdout,
					"[AmbEncounter] locus=%d reason=%s curr_idx=%d range=[%d,%d] diag_expected=%d diag_advanced=%d\n",
					curr_abs_locus,
					amb_reason,
					curr_abs_ambiguous,
					ambiguous_first,
					ambiguous_last,
					diag_expected_amb_sites,
					diag_advanced_amb);
			}
		}
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
		if (can_advance_amb || sib_advance_amb) {
			curr_abs_ambiguous++;
			diag_advanced_amb++;
		}
		trace_ambiguous_cursor("fwd_post", curr_abs_locus, is_sibling, expected_delta);
		if (has_amb_range && (curr_abs_ambiguous < ambiguous_first || curr_abs_ambiguous > ambiguous_last + 1)) {
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
		if (curr_segment_locus >= G->Lengths_bio[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}

		// End-of-window diagnostic: verify ambiguous-site bookkeeping parity
		if (supersite_trace_enabled()) {
			const bool has_amb_range = (ambiguous_first <= ambiguous_last);
			int slots = has_amb_range ? (ambiguous_last - ambiguous_first + 1) : 0;
			if (diag_expected_amb_sites > slots) {
				std::fprintf(stderr, "[ss-amb-diag][single] expected_amb_sites=%d slots=%d\n", diag_expected_amb_sites, slots);
				assert(false && "more ambiguous data sites than available slots");
			}
			if (diag_expected_amb_sites != diag_advanced_amb) {
				std::fprintf(stderr, "[ss-amb-diag][single] expected_amb_sites=%d advanced_amb=%d slots=%d\n",
					 diag_expected_amb_sites, diag_advanced_amb, slots);
				assert(false && "ambiguous cursor advanced count mismatch");
			}
			else {
				if (supersite_trace_enabled()) {
					std::fprintf(stdout, "[ss-amb-diag-PASS][single] expected_amb_sites=%d advanced_amb=%d slots=%d\n",
						 diag_expected_amb_sites, diag_advanced_amb, slots);
				}
			}
		}
	}
		if (curr_segment_index > segment_last) {
			curr_segment_index = segment_last;
			curr_segment_locus = (segment_last >= segment_first) ? G->Lengths_bio[segment_last] - 1 : 0;
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
	curr_segment_locus = G->Lengths_bio[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;
	prev_abs_locus = locus_last;
	trace_forward_active = false;
	trace_backward_active = true;
	const bool has_transition_buffer = !transition_probabilities.empty();
	const size_t expected_transitions = static_cast<size_t>(G->n_transitions);
	if (!has_transition_buffer) {
		assert(expected_transitions == 0);
	} else {
		assert(transition_probabilities.size() >= expected_transitions);
	}
	trace_backward_table_header_emitted = false;
	const bool trans_trace = trans_trace_enabled_s();
	const char* trans_trace_sample = trans_trace_sample_s();

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx, panel_codes_size);

	// Flag: backward pass always starts with initialization; if locus_last is a sibling, defer until first non-sibling
	bool need_init = true;
	bool pending_collapse = false;
	yt = 0.0f;
	nt = 1.0f;

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		if (curr_segment_locus == G->Lengths_bio[curr_segment_index] - 1) pending_collapse = true;
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		const int prev_before = prev_abs_locus;
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
			if (supersite_trace_enabled()) {
				std::fprintf(stdout, "[DEFERRED_INIT] locus=%d is_sibling=1, advancing counters and continue\n", curr_abs_locus);
			}
			// Minimal interaction: just advance the counters as if this locus was processed, but do no HMM math.
			// This ensures that when the first anchor is reached, prev_abs_locus is equal to it, making yt=0,
			// which is correct for the start of a backward pass.
			INIT_SIB(site_view);
			// CRITICAL: Check for segment boundary BEFORE decrementing, since we'll skip the normal check with continue
			if (curr_segment_locus == 0 && curr_abs_locus != locus_first && curr_segment_index > segment_first) {
				if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
					unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
					unsigned int prev_dipcount = (curr_segment_index > 0) ? G->countDiplotypes(G->Diplotypes[curr_segment_index-1]) : 0;
					unsigned int n_trans = curr_dipcount * prev_dipcount;
					std::fprintf(stdout,
					             "[TRANS_TRACE][single][skip-sib] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d buf_size=%zu\n",
					             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, transition_probabilities.size());
					std::fflush(stdout);
				}
				int ret = SET_OTHER_TRANS(transition_probabilities);
				if (ret < 0) return ret;
				else n_underflow_recovered += ret;
			}
			// Siblings don't count toward segment length, so don't decrement curr_segment_locus
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
		if (supersite_trace_enabled()) {
			std::fprintf(stderr, "[YT_DEBUG] locus=%d is_sibling=%d is_anchor=%d prev_abs=%d locus_last=%d yt_before=%.10f\n",
			             curr_abs_locus, is_sibling, is_anchor, prev_abs_locus, locus_last, yt);
		}
		if (!is_sibling) {
			yt = (curr_abs_locus == locus_last || curr_abs_locus == prev_abs_locus) ? 0.0f : M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
			nt = 1.0f - yt;
			if (supersite_trace_enabled()) {
				std::fprintf(stderr, "[YT_DEBUG] locus=%d COMPUTED yt=%.10f (curr==last:%d curr==prev:%d)\n",
				             curr_abs_locus, yt, (curr_abs_locus == locus_last), (curr_abs_locus == prev_abs_locus));
			}
		} else {
			if (supersite_trace_enabled()) {
				std::fprintf(stderr, "[YT_DEBUG] locus=%d SKIPPED (is_sibling=1), yt unchanged=%.10f\n",
				             curr_abs_locus, yt);
			}
		}
		trace_ambiguous_cursor("bwd_pre", curr_abs_locus, is_sibling, 0);

		// Deferred initialization: first non-sibling after starting on a sibling - use INIT
		if (need_init && is_anchor) {
			if (supersite_trace_enabled()) {
				std::fprintf(stdout, "[DEFERRED_INIT] locus=%d is_anchor=1, need_init=true, using INIT_FROM_MASK\n", curr_abs_locus);
			}
			need_init = false;
			pending_collapse = false;
			// FIX: Set prev_abs_locus to curr_abs_locus to ensure yt=0 for the first anchor
			prev_abs_locus = curr_abs_locus;
			if (M.ss_anchor_split_emissions) {
				supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
		} else if (need_init && !is_sibling) {
			// Deferred initialization: first biallelic after starting on a sibling - use INIT
			need_init = false;
			pending_collapse = false;
			if (hmm_mis) {
				INIT_MIS();
			} else {
				bial_adapter.build_match_mask(site_view, n_cond_haps, curr_rel_locus + curr_rel_locus_offset, init_match_mask);
				INIT_FROM_MASK(init_match_mask, static_cast<float>(M.ed/M.ee));
			}
		} else if (!pending_collapse) {
			if (is_anchor) {
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
				if (M.ss_anchor_split_emissions) {
					supersite_adapter.build_match_mask(site_view, n_cond_haps, /*use_anchor_split_semantics*/false, init_match_mask);
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
		if (supersite_trace_enabled()) {
			fprintf(stderr, "[SIBLING_SKIP_TRACE] loc=%d, prev_loc_in=%d, kind=%d, is_sib=%d, update_prev=%d, yt=%.8f, prev_loc_out=%d\n",
					curr_abs_locus,
					prev_before,
					(int)site_view.kind,
					(int)is_sibling,
					(int)update_prev_locus,
					yt,
					prev_after);

			// Calculate biological anchor index (for supersites, count only anchors, not siblings)
			int bio_anchor_idx = curr_abs_locus;
			if (super_sites && locus_to_super_idx) {
				int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
				if (ss_idx >= 0) {
					// Map to biological anchor index (count only anchors, not siblings)
					bio_anchor_idx = 0;
					for (int i = 0; i < curr_abs_locus; i++) {
						int check_ss = (*locus_to_super_idx)[i];
						if (check_ss < 0 || i == (*super_sites)[check_ss].global_site_id) {
							bio_anchor_idx++;
						}
					}
				}
			}

			// Enhanced BWD_PROB_DETAIL logging with bio_anchor, segment info, and first 8 prob[] values
			fprintf(stderr, "[BWD_PROB_DETAIL] sample=%s locus=%d bio_anchor=%d is_sib=%d seg=%d seg_locus=%d probSumT=%.15f\n",
			        G->name.c_str(), curr_abs_locus, bio_anchor_idx, (int)is_sibling,
			        curr_segment_index, curr_segment_locus, (double)probSumT);
			fprintf(stderr, "  prob[0][0-7]= ");
			for (int i = 0; i < std::min(8, (int)prob.size()); i++) {
				fprintf(stderr, "%.8e ", (double)prob[i]);
			}
			fprintf(stderr, "\n");
		}
		prev_abs_locus = prev_after;
		trace_log_backward_state(curr_abs_locus,
		                         prev_before,
		                         prev_after,
		                         static_cast<double>(yt),
		                         static_cast<double>(nt),
		                         update_prev_locus,
		                         is_anchor,
		                         is_sibling,
		                         hmm_amb,
		                         hmm_mis,
		                         hmm_hom);

			// Only emit transition probabilities if caller provided storage
			if (!transition_probabilities.empty()) {
				if (curr_abs_locus == 0) {
					if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
						unsigned int n_trans = G->countDiplotypes(G->Diplotypes[0]);
						std::fprintf(stdout,
						             "[TRANS_TRACE][single][first] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d buf_size=%zu\n",
						             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, transition_probabilities.size());
						std::fflush(stdout);
					}
					SET_FIRST_TRANS(transition_probabilities);
				}
				if (!is_sibling && curr_segment_locus == 0 && curr_abs_locus != locus_first && curr_segment_index > segment_first) {
					if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
						unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
						unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
						unsigned int n_trans = curr_dipcount * prev_dipcount;
						std::fprintf(stdout,
						             "[TRANS_TRACE][single][other] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d buf_size=%zu\n",
						             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, transition_probabilities.size());
						std::fflush(stdout);
					}
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
					if (supersite_trace_enabled()) {
						std::fprintf(stdout, "SCTrace anchor=%d rel_missing=%d\n", curr_abs_locus, idx);
						std::fprintf(stdout, "AlphaSumMissing:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", AlphaSumMissing[idx][h]);
						std::fprintf(stdout, "\nBeta:");
						for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", prob[h + idx * HAP_NUMBER]);
						std::fprintf(stdout, "\n");
					}
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
				if (supersite_trace_enabled()) std::fprintf(stdout, "IMPUTE-guard locus=%d map_i=%d rel_missing_idx=%d n_missing=%d\n", curr_abs_locus, map_i, rel_missing_idx, n_missing);
				if (rel_missing_idx >= 0) {
					// Forward recorded a missing slot for this locus -> perform imputation
					IMPUTE(missing_probabilities);
				} else {
					// Nothing recorded in forward (e.g., sibling or no-missing window) -> skip IMPUTE
					if (supersite_trace_enabled()) std::fprintf(stdout, "Skipping IMPUTE for locus=%d rel_missing_idx=%d n_missing=%d\n", curr_abs_locus, rel_missing_idx, n_missing);
				}
			}
			curr_abs_missing--;
		}


		// Only decrement segment locus for non-sibling loci (siblings don't count toward segment length)
		if (!is_sibling) {
			curr_segment_locus--;
		}
		const bool has_amb_range = (ambiguous_first <= ambiguous_last);
	const bool can_retreat_amb = data_amb && has_amb_range && (curr_abs_ambiguous >= ambiguous_first);
	// Sibling variants do not have entries in the Ambiguous array, so curr_abs_ambiguous should not retreat for them.
	const bool sib_retreat_amb = false; //is_sibling && hmm_amb && has_amb_range && (curr_abs_ambiguous > ambiguous_first);
	const int expected_delta = (can_retreat_amb) ? -1 : 0;
		const int cursor_before_bwd = curr_abs_ambiguous;
	if (expected_delta) curr_abs_ambiguous--;
		trace_ambiguous_cursor("bwd_post", curr_abs_locus, is_sibling, expected_delta);
		const int lower_bound = ambiguous_first - 1;
		if (has_amb_range && (curr_abs_ambiguous < lower_bound || curr_abs_ambiguous > ambiguous_last)) {
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
			curr_segment_locus = G->Lengths_bio[curr_segment_index] - 1;
		}
	}
		curr_segment_index = segment_first;
		curr_segment_locus = 0;
		curr_abs_locus = locus_first;
	trace_backward_active = false;
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

void haplotype_segment_single::build_lane_priors_first(const SiteView& site_view, double lane_probs[HAP_NUMBER], bool use_outer_prod) const {
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

	if (use_outer_prod && supersites_enabled_flag && !AlphaSum.empty()) {
		const int rel_seg = 0;
		double lane_total = 0.0;
		for (int h = 0; h < HAP_NUMBER; ++h) {
			double w = lane_probs[h] * static_cast<double>(AlphaSum[rel_seg][h]);
			lane_probs[h] = w;
			lane_total += w;
		}
		if (lane_total > std::numeric_limits<double>::min()) {
			const double inv = 1.0 / lane_total;
			for (int h = 0; h < HAP_NUMBER; ++h) lane_probs[h] *= inv;
		}
	}
}

void haplotype_segment_single::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	const unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	std::vector<double> cprobs(n_transitions, 0.0);

	const bool trace_trans = []() {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
		return env && env[0] != '\0' && env[0] != '0';
	}();
	if (supersite_trace_enabled()) {
		std::fprintf(stdout, "SET_FIRST_TRANS single: n_transitions=%u buf_size=%zu\n",
					 n_transitions, transition_probabilities.size());
	}

	// Loudly enforce buffer sizing instead of skipping writes.
	assert(n_transitions > 0);
	assert(transition_probabilities.size() >= n_transitions);

	double lane_probs[HAP_NUMBER];
	const bool allow_outer = []() {
		const char* env = std::getenv("SHAPEIT5_SS_USE_OUTER");
		return env && env[0] != '\0' && env[0] != '0';
	}();

	// Build SiteView for the first locus to derive lane_class
	SiteView sv{};
	bool has_supersite = false;
	if (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx) {
		SupersiteEmissionAdapter sup_adapt(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx, panel_codes_size);
		has_supersite = sup_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	if (!has_supersite) {
		BiallelicEmissionAdapter bial_adapt(G, &Hvar);
		bial_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	build_lane_priors_first(sv, lane_probs, allow_outer);

	if (trace_trans && G && G->n_segments > 1 && (segment_first == 0)) {
		std::fprintf(stderr, "[SET_FIRST_TRANS_DEBUG] sample=%s use_outer=%d probSumT=%.15f\n",
		             G->name.c_str(), (int)allow_outer, (double)probSumT);
		std::fprintf(stderr, "  lane_probs:");
		for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " %.9f", lane_probs[h]);
		std::fprintf(stderr, "\n");
		std::fprintf(stderr, "  probSumH:");
		for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " %.9f", (double)probSumH[h]);
		std::fprintf(stderr, "\n");
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
    // LOG HMM STATE - Check values before and after TRANS_DIP_MULT
    bool is_test_sample = (G->name.find("_sample") != std::string::npos);
    bool is_seg2 = (curr_segment_index == 2);
    if (is_test_sample && is_seg2) {
        std::fprintf(stderr, "[HMM_BEFORE_TRANS_DIP_MULT] sample=%s locus=%d seg=%d sumHProbs=%.20g (hex:%a)\n",
                     G->name.c_str(), curr_abs_locus, curr_segment_index, (double)sumHProbs, (double)sumHProbs);
    }

    if (TRANS_DIP_MULT()) {
        if (TRANS_DIP_ADD()) return -2;
        else underflow_recovered = 1;
    }

    // LOG DPROBS AFTER TRANS_DIP_MULT
    if (is_test_sample && is_seg2) {
        std::fprintf(stderr, "[HMM_AFTER_TRANS_DIP_MULT] sample=%s locus=%d seg=%d sumDProbs=%.20g (hex:%a)\n",
                     G->name.c_str(), curr_abs_locus, curr_segment_index, sumDProbs, sumDProbs);
    }

	unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
	unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index-1]);
	unsigned int n_transitions = curr_dipcount * prev_dipcount;
	double scaleDip = 1.0 / sumDProbs;

	if (trans_parity_trace_enabled_s()) {
		std::fprintf(stderr, "[TRANS_DIP_ADD debug][single] locus=%d seg=%d curr_abs_transition=%d sumHProbs=%g sumDProbs=%g\n",
		             curr_abs_locus, curr_segment_index, curr_abs_transition, sumHProbs, sumDProbs);
		for (int h = 0; h < std::min(4, HAP_NUMBER); ++h) {
			std::fprintf(stderr, "  HProbs[%d]:", h);
			for (int hh = 0; hh < HAP_NUMBER; ++hh) std::fprintf(stderr, " %.6g", HProbs[h*HAP_NUMBER + hh]);
			std::fprintf(stderr, "\n");
		}
		std::fflush(stderr);
	}
	bool enable_debug = (!debug::SUPERDEBUG_SAMPLENAME.empty() && G->name == debug::SUPERDEBUG_SAMPLENAME) ||
						(G->name.find("_sample") != std::string::npos);
	if (enable_debug) {
		// Determine biological anchor index (for supersites, map to anchor; for biallelic, locus = anchor)
		int bio_anchor_idx = curr_abs_locus;
		bool is_supersite_anchor = false;
		if (super_sites && locus_to_super_idx) {
			int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
			if (ss_idx >= 0) {
				const SuperSite& ss = (*super_sites)[ss_idx];
				is_supersite_anchor = (curr_abs_locus == ss.global_site_id);
				// Map to biological anchor index (count only anchors, not siblings)
				bio_anchor_idx = 0;
				for (int i = 0; i < curr_abs_locus; i++) {
					int check_ss = (*locus_to_super_idx)[i];
					if (check_ss < 0 || i == (*super_sites)[check_ss].global_site_id) {
						bio_anchor_idx++;
					}
				}
			}
		}
		std::fprintf(stderr, "[SET_OTHER_TRANS] sample=%s locus=%d bio_anchor=%d seg=%d is_ss_anchor=%d n_trans=%u sumHProbs=%.15f sumDProbs=%.15f\n",
					 G->name.c_str(), curr_abs_locus, bio_anchor_idx, curr_segment_index, (int)is_supersite_anchor, n_transitions, sumHProbs, sumDProbs);
		// Show first 8 transition probabilities
		for (unsigned int t = 0; t < std::min(8u, n_transitions); t++) {
			std::fprintf(stderr, "  trans[%u] DProb=%.15f final=%.15f\n", t, DProbs[t], DProbs[t] * scaleDip);
		}
	}
	if (supersite_trace_enabled()) {
		std::fprintf(stdout, "SET_OTHER_TRANS single: curr_abs_transition=%d n_transitions=%u buf_size=%zu\n",
					 curr_abs_transition, n_transitions, transition_probabilities.size());
	}

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

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
		if (G->ss_observed_gts.empty()) {
			G->snapshotSupersiteObservedGts(*super_sites, *super_site_var_index);
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
                if (matrix_trace && matrix_trace[0] != '\0' && matrix_trace[0] != '0' &&
                    (ss_idx == 1 || ss_idx == 6 || ss_idx == 3)) {
                std::fprintf(stderr, "[MATRIX_POPULATE] ss_idx=%zu (anchor=%u) panel_codes first 16: ",
                                 ss_idx, ss.global_site_id);
                for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
                    std::fprintf(stderr, "%u ", (unsigned)ss_panel_matrix[ss_idx][k]);
                }
                std::fprintf(stderr, "\n");

                // Cross-check packed supersite codes against Hvar (biallelic view)
                const bool anchor_in_window = (ss.global_site_id >= (unsigned)locus_first) &&
                                              (ss.global_site_id <= (unsigned)locus_last);
                const int rel_row = anchor_in_window
                                      ? static_cast<int>(ss.global_site_id - locus_first + curr_rel_locus_offset)
                                      : -1;
                std::fprintf(stderr,
                             "[MATRIX_POPULATE] anchor_in_window=%d rel_row=%d locus_first=%d locus_last=%d offset=%d\n",
                             anchor_in_window ? 1 : 0, rel_row, locus_first, locus_last, curr_rel_locus_offset);
                const unsigned int dump_limit = std::min(32u, n_cond_haps);
                for (unsigned int k = 0; k < dump_limit; ++k) {
                    const unsigned int gh = (*cond_idx)[k];
                    const uint8_t packed_code = ss_panel_matrix[ss_idx][k];
                    const int hvar = anchor_in_window ? static_cast<int>(Hvar.get(rel_row, k)) : -1;
                    std::fprintf(stderr,
                                 "[PACK_CHK] ss_idx=%zu anchor=%u k=%u cond_idx=%u packed=%u hvar=%d\n",
                                 ss_idx, ss.global_site_id, k, gh, (unsigned)packed_code, hvar);
                }
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
	AlphaLaneSum.clear();
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
			__m256 lane_vec = _mm256_load_ps(&probSumH[0]);
			_mm256_store_ps(AlphaLaneSum[rel_seg].lane, lane_vec);
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
	const bool trans_trace = trans_trace_enabled_s();
	const char* trans_trace_sample = trans_trace_sample_s();
	const int seg_count = segment_last - segment_first + 1;
	std::vector<uint8_t> seg_has_anchor(seg_count, 0);
	std::vector<uint8_t> seg_has_sibling(seg_count, 0);
	std::vector<uint8_t> seg_wrote_transition(seg_count, 0);
	std::vector<std::vector<int>> seg_anchor_ss(seg_count);
	std::vector<std::vector<int>> seg_sibling_ss(seg_count);
	std::vector<unsigned int> seg_non_sib_count(seg_count, 0);
	std::vector<int> seg_first_non_sib_locus(seg_count, -1);
	std::vector<int> seg_last_non_sib_locus(seg_count, -1);
	bool wrote_first_transition = false;
	// TEMP DEBUG (TRANS_MISMATCH_ERROR): one-off trace to capture skipped segment boundaries.
	// Remove after confirming why transition counts diverge in supersite windows.
	const bool trans_mismatch_trace = std::getenv("SHAPEIT5_TRACE_TRANS_MISMATCH");
	static bool trans_mismatch_trace_dumped = false;
	struct trans_mismatch_event {
		int locus;
		int seg_idx;
		int seg_locus;
		unsigned int n_trans;
		const char* reason;
	};
	std::vector<trans_mismatch_event> trans_mismatch_events;
	if (trans_mismatch_trace && !trans_mismatch_trace_dumped) trans_mismatch_events.reserve(16);

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
		const int rel_seg = curr_segment_index - segment_first;
		if (rel_seg >= 0 && rel_seg < seg_count) {
			if (is_sibling) seg_has_sibling[rel_seg] = 1;
			else seg_has_anchor[rel_seg] = 1;
			if (!is_sibling) {
				seg_non_sib_count[rel_seg] += 1u;
				if (seg_first_non_sib_locus[rel_seg] < 0) seg_first_non_sib_locus[rel_seg] = curr_abs_locus;
				seg_last_non_sib_locus[rel_seg] = curr_abs_locus;
			}
			const int ss_idx = site_view.supersite_index;
			auto add_unique = [](std::vector<int>& vec, int val) {
				if (val < 0) return;
				for (int existing : vec) {
					if (existing == val) return;
				}
				vec.push_back(val);
			};
			if (is_anchor) add_unique(seg_anchor_ss[rel_seg], ss_idx);
			else if (is_sibling) add_unique(seg_sibling_ss[rel_seg], ss_idx);
		}

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
				const bool at_segment_boundary = (curr_segment_locus == 0 && curr_segment_index > segment_first);
				if (at_segment_boundary && trans_mismatch_trace && !trans_mismatch_trace_dumped) {
					const unsigned int curr_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index]);
					const unsigned int prev_dipcount = G->countDiplotypes(G->Diplotypes[curr_segment_index - 1]);
					const unsigned int n_trans = curr_dipcount * prev_dipcount;
					const char* reason = nullptr;
					if (curr_abs_locus == locus_first) reason = "window_start";
					else if (is_sibling) reason = "sibling_boundary";
					if (reason && trans_mismatch_events.size() < 16) {
						trans_mismatch_events.push_back({curr_abs_locus, curr_segment_index, curr_segment_locus, n_trans, reason});
					}
				}
				if (curr_abs_locus == 0) {
					if (trans_trace && (!trans_trace_sample || G->name == trans_trace_sample)) {
						unsigned int n_trans = G->countDiplotypes(G->Diplotypes[0]);
						std::fprintf(stdout,
						             "[TRANS_TRACE][single][first] sample=%s curr_abs_locus=%d seg_idx=%d n_trans=%u curr_abs_transition=%d buf_size=%zu\n",
						             G->name.c_str(), curr_abs_locus, curr_segment_index, n_trans, curr_abs_transition, transition_probabilities.size());
					std::fflush(stdout);
				}
					SET_FIRST_TRANS(transition_probabilities);
					wrote_first_transition = true;
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
					if (rel_seg >= 0 && rel_seg < seg_count) seg_wrote_transition[rel_seg] = 1;
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
	if (!transition_probabilities.empty()) {
		bool error = false;
		if (segment_first == 0 && locus_first == 0 && !wrote_first_transition) {
			std::fprintf(stderr, "[TRANS_MISMATCH_ERROR] Sample=%s: Missing first-segment transition write in window\n",
			             G->name.c_str());
			error = true;
		}
		for (int s = segment_first; s <= segment_last; ++s) {
			const int rel = s - segment_first;
			if (rel < 0 || rel >= seg_count) continue;
			if (!seg_sibling_ss[rel].empty()) {
				for (int ss_idx : seg_sibling_ss[rel]) {
					bool has_anchor = false;
					for (int anchor_idx : seg_anchor_ss[rel]) {
						if (anchor_idx == ss_idx) {
							has_anchor = true;
							break;
						}
					}
					if (!has_anchor) {
						std::fprintf(stderr,
						             "[TRANS_MISMATCH_ERROR] Sample=%s: Segment %d has sibling ss_idx=%d without its anchor\n",
						             G->name.c_str(), s, ss_idx);
						error = true;
					}
				}
			}
		}
		for (int s = segment_first + 1; s <= segment_last; ++s) {
			const int rel = s - segment_first;
			if (rel < 0 || rel >= seg_count) continue;
			if (seg_has_anchor[rel] && !seg_wrote_transition[rel]) {
				const unsigned int len_bio = (s >= 0 && s < (int)G->Lengths_bio.size()) ? G->Lengths_bio[s] : 0u;
				std::fprintf(stderr,
				             "[TRANS_MISMATCH_ERROR] Sample=%s: Missing transition write for segment %d in window (len_bio=%u non_sib=%u first_non_sib=%d last_non_sib=%d)\n",
				             G->name.c_str(), s, len_bio, seg_non_sib_count[rel], seg_first_non_sib_locus[rel], seg_last_non_sib_locus[rel]);
				error = true;
			}
		}
		if (error) {
			std::fprintf(stderr, "  segments=[%d,%d] locus_range=[%d,%d]\n",
			             segment_first, segment_last, locus_first, locus_last);
			std::fprintf(stderr, "  Whole-sample transitions=%u, buffer size=%zu\n",
			             G->n_transitions, transition_probabilities.size());
			// TEMP DEBUG (TRANS_MISMATCH_ERROR): dump the first few skipped boundaries once.
			if (trans_mismatch_trace && !trans_mismatch_trace_dumped) {
				std::fprintf(stderr, "[TRANS_MISMATCH_TRACE][single] sample=%s skipped_boundaries=%zu (max 16)\n",
				             G->name.c_str(), trans_mismatch_events.size());
				for (const auto& ev : trans_mismatch_events) {
					std::fprintf(stderr,
					             "  [TRANS_MISMATCH_TRACE][single] reason=%s locus=%d seg=%d seg_locus=%d n_trans=%u\n",
					             ev.reason, ev.locus, ev.seg_idx, ev.seg_locus, ev.n_trans);
				}
				trans_mismatch_trace_dumped = true;
			}
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
		SupersiteEmissionAdapter sup_adapt(G, super_sites, locus_to_super_idx, super_site_var_index);
		has_supersite = sup_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	if (!has_supersite) {
		BiallelicEmissionAdapter bial_adapt(G, &Hvar);
		bial_adapt.build_view(locus_first, ambiguous_first, sv);
	}
	build_lane_priors_first(sv, lane_probs, allow_outer);

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

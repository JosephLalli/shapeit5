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

haplotype_segment_single::haplotype_segment_single(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M,
    const std::vector<SuperSite>* _super_sites,
    const std::vector<bool>* _is_super_site,
    const std::vector<int>* _locus_to_super_idx,
    const uint8_t* _panel_codes,
    const std::vector<int>* _super_site_var_index) :
    G(_G), M(_M), super_sites(_super_sites), is_super_site(_is_super_site),
    locus_to_super_idx(_locus_to_super_idx), panel_codes(_panel_codes), super_site_var_index(_super_site_var_index), cond_idx(&idxH),
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
	AlphaLaneSum = vector < aligned_vector32 < float > > (segment_last - segment_first + 1, aligned_vector32 < float > (HAP_NUMBER, 0.0f));
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

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx);
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
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb;
		yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		if (supersite_trace_enabled()) {
			std::fprintf(stdout,
					"FWD locus=%d rel=%d seg_locus=%d kind=%d emit=%d amb=%d mis=%d\n",
					curr_abs_locus, curr_rel_locus, curr_segment_locus,
					(int)site_view.kind, (int)emit, (int)hmm_amb, (int)hmm_mis);
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
                // Sibling within window: true no-op (preserve probability state)
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
                // Sibling at segment boundary: true no-op (preserve probability state)
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
			Alpha[curr_segment_index - segment_first] = prob;
			AlphaSum[curr_segment_index - segment_first] = probSumH;
			AlphaLaneSum[curr_segment_index - segment_first] = probSumH;
			AlphaSumSum[curr_segment_index - segment_first] = probSumT;
			AlphaLocus[curr_segment_index - segment_first] = prev_abs_locus;
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
		curr_abs_ambiguous += data_amb;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
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

	const bool supersites_enabled = (super_sites && locus_to_super_idx && super_site_var_index && panel_codes && cond_idx);
	BiallelicEmissionAdapter bial_adapter(G, &Hvar);
	SupersiteEmissionAdapter supersite_adapter(G, super_sites, locus_to_super_idx, super_site_var_index, panel_codes, cond_idx);

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
		const bool data_mis = hmm_mis && !is_sibling;
		const bool data_amb = hmm_amb;
		yt = (curr_abs_locus == locus_last)?0.0:M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

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
					// Sibling within window (backward): true no-op (preserve probability state)
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
					// Sibling at segment boundary (backward): true no-op (preserve probability state)
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
		curr_abs_ambiguous -= data_amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
	return n_underflow_recovered;
}

void haplotype_segment_single::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
	double scale = 1.0f / probSumT, scaleDip = 0.0f;
	unsigned int n_transitions = G->countDiplotypes(G->Diplotypes[0]);
	vector < double > cprobs = vector < double > (n_transitions, 0.0);
	for (unsigned int d = 0, t = 0 ; d < 64 ; ++d) {
		if (DIP_GET(G->Diplotypes[0], d)) {
			cprobs[t] = (double)(probSumH[DIP_HAP0(d)]*scale) * (double)(probSumH[DIP_HAP1(d)]*scale);
			scaleDip += cprobs[t++];
		}
	}
	scaleDip = 1.0f / scaleDip;
	for (unsigned int t = 0 ; t < n_transitions ; t ++) transition_probabilities[t] = cprobs[t] * scaleDip;
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

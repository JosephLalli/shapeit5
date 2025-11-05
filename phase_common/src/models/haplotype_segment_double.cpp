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
#include <mutex>
#include <cstdio>
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

// EXPERIMENTAL: COLLAPSE transition normalization testing (SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION)
// TODO: Remove before release after determining optimal behavior (see Bug #4)
static bool normalize_collapse_transition_enabled_d() {
    static int flag = -1;
    if (flag < 0) {
        const char* env = std::getenv("SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION");
        flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return flag == 1;
}
} // namespace

haplotype_segment_double::haplotype_segment_double(genotype * _G, bitmatrix & H, vector < unsigned int > & idxH, window & W, hmm_parameters & _M,
    const std::vector<SuperSite>* _super_sites,
    const std::vector<bool>* _is_super_site,
    const std::vector<int>* _locus_to_super_idx,
    const uint8_t* _panel_codes,
    const std::vector<int>* _super_site_var_index) :
    G(_G), M(_M), super_sites(_super_sites), is_super_site(_is_super_site),
    locus_to_super_idx(_locus_to_super_idx), panel_codes(_panel_codes), super_site_var_index(_super_site_var_index), cond_idx(&idxH) {
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
	if (super_sites && locus_to_super_idx && panel_codes && super_site_var_index) {
		ss_cond_codes = aligned_vector32<uint8_t>(n_cond_haps, 0);
		ss_emissions = aligned_vector32<double>(n_cond_haps, 1.0);
		ss_emissions_h1 = aligned_vector32<double>(n_cond_haps, 1.0);
		ss_cached.resize(super_sites->size(), false);  // Phase 3: initialize cache flags
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

	for (curr_abs_locus = locus_first ; curr_abs_locus <= locus_last ; curr_abs_locus++) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		bool update_prev_locus = true;
		char rare_allele = M.rare_allele[curr_abs_locus];
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);
		genotype::SuperSiteContext ss_ctx = G->getSuperSiteContext(curr_abs_locus);
		if (ss_ctx.is_member) {
			if (!ss_ctx.is_anchor) {
				amb = false;
				mis = false;
				hom = true;
			} else {
				amb = ss_ctx.has_het || ss_ctx.has_sca;
				mis = ss_ctx.all_missing;
				hom = !(amb || mis);
			}
		}
		yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		if (curr_rel_locus == 0) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != 0) {
			if (hom) update_prev_locus = RUN_HOM(rare_allele);
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else  COLLAPSE_MIS();
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
            Alpha[curr_segment_index - segment_first] = prob;
            AlphaSum[curr_segment_index - segment_first] = probSumH;
            AlphaSumSum[curr_segment_index - segment_first] = probSumT;
            AlphaLocus[curr_segment_index - segment_first] = prev_abs_locus;
        }
		if (mis) {
			AlphaMissing[curr_rel_missing] = prob;
			AlphaSumMissing[curr_rel_missing] = probSumH;
			curr_abs_missing ++;
		}

		curr_segment_locus ++;
		curr_abs_ambiguous += amb;
		if (curr_segment_locus >= G->Lengths[curr_segment_index]) {
			curr_segment_index++;
			curr_segment_locus = 0;
		}
	}
}

int haplotype_segment_double::backward(vector < double > & transition_probabilities, vector < float > & missing_probabilities,
                                       vector < float > * SC, const vector < bool > * anchor_has_missing) {
	int n_underflow_recovered = 0;
	curr_segment_index = segment_last;
	curr_segment_locus = G->Lengths[segment_last] - 1;
	curr_abs_ambiguous = ambiguous_last;
	curr_abs_missing = missing_last;
	curr_abs_transition = transition_last;
	prev_abs_locus = locus_last;

	for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
		curr_rel_locus = curr_abs_locus - locus_first;
		curr_rel_missing = curr_abs_missing - missing_first;
		char rare_allele = M.rare_allele[curr_abs_locus];
		bool update_prev_locus = true;
		bool amb = VAR_GET_AMB(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool mis = VAR_GET_MIS(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
		bool hom = !(amb || mis);
		genotype::SuperSiteContext ss_ctx = G->getSuperSiteContext(curr_abs_locus);
		if (ss_ctx.is_member) {
			if (!ss_ctx.is_anchor) {
				amb = false;
				mis = false;
				hom = true;
			} else {
				amb = ss_ctx.has_het || ss_ctx.has_sca;
				mis = ss_ctx.all_missing;
				hom = !(amb || mis);
			}
		}
		yt = (curr_abs_locus == locus_last)?0.0:M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
		nt = 1.0f - yt;

		if (curr_abs_locus == locus_last) {
			if (hom) INIT_HOM();
			else if (amb) INIT_AMB();
			else INIT_MIS();
		} else if (curr_segment_locus != G->Lengths[curr_segment_index] - 1) {
			if (hom) update_prev_locus = RUN_HOM(rare_allele);
			else if (amb) RUN_AMB();
			else RUN_MIS();
		} else {
			if (hom) COLLAPSE_HOM();
			else if (amb) COLLAPSE_AMB();
			else COLLAPSE_MIS();
		}
		if (curr_segment_locus == 0) SUMK();
		prev_abs_locus=update_prev_locus?curr_abs_locus:prev_abs_locus;

		if (curr_abs_locus == 0) SET_FIRST_TRANS(transition_probabilities);
		if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
			int ret = SET_OTHER_TRANS(transition_probabilities);
			if (ret < 0) return ret;
			else n_underflow_recovered += ret;
		}

		if (mis) {
			// Phase 3: Check if this is a supersite with all members missing
			int ss_idx_here = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[curr_abs_locus] : -1;
			
			if (ss_idx_here >= 0 && anchor_has_missing && (*anchor_has_missing)[ss_idx_here] && SC) {
				// Part of a missing supersite
				const SuperSite& ss = (*super_sites)[ss_idx_here];
				
				if (curr_abs_locus == (int)ss.global_site_id) {
					// Anchor: compute multivariant for entire supersite
					IMPUTE_SUPERSITE_MULTIVARIATE(*SC, ss, ss_idx_here);
				}
				// Else: sibling, skip (no IMPUTE call)
				
				curr_abs_missing--;  // Still decrement counter
			} else {
				// Normal biallelic missing site
				IMPUTE(missing_probabilities);
				curr_abs_missing--;
			}
		}


		curr_segment_locus--;
		curr_abs_ambiguous -= amb;
		if (curr_segment_locus < 0 && curr_segment_index > 0) {
			curr_segment_index--;
			curr_segment_locus = G->Lengths[curr_segment_index] - 1;
		}
	}
	return n_underflow_recovered;
}

void haplotype_segment_double::SET_FIRST_TRANS(vector < double > & transition_probabilities) {
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

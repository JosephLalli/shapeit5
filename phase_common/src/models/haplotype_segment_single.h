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

#ifndef _HAPLOTYPE_SEGMENT_SINGLE_H
#define _HAPLOTYPE_SEGMENT_SINGLE_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>
#include <models/site_emission_types.h>
#include <models/super_site_emissions.h>
#include <models/supersite_trace_utils.h>
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <cstring>

#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>


inline bool trans_parity_trace_enabled_s() {
    static int flag = -1;
    if (flag < 0) {
        const char* env = std::getenv("SHAPEIT5_DEBUG_TRANS_PARITY");
        flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    return flag == 1;
}

// Structured imputation trace (TSV) helpers, gated by SHAPEIT5_IMPUTE_TSV
inline FILE* impute_tsv_file() {
    static bool init = false;
    static FILE* file = nullptr;
    if (!init) {
        init = true;
        const char* path = std::getenv("SHAPEIT5_IMPUTE_TSV");
        if (path && path[0]) {
            if (std::strcmp(path, "-") == 0) {
                file = stdout;
            } else {
                file = std::fopen(path, "w");
                if (!file) {
                    std::fprintf(stderr, "[IMPUTE_TSV] failed to open %s\n", path);
                }
            }
        }
    }
    return file;
}

inline const char* impute_tsv_sample_filter() {
    static bool init = false;
    static const char* sample = nullptr;
    if (!init) {
        init = true;
        const char* env = std::getenv("SHAPEIT5_IMPUTE_TSV_SAMPLE");
        sample = (env && env[0]) ? env : nullptr;
    }
    return sample;
}

inline bool impute_tsv_should_log(const genotype* G) {
    FILE* f = impute_tsv_file();
    if (!f) return false;
    const char* filter = impute_tsv_sample_filter();
    if (!filter) return true;
    return (G != nullptr && G->name == filter);
}

inline void impute_tsv_log_row(const char* mode,
                               const genotype* G,
                               int locus,
                               int bio_locus,
                               int ss_idx,
                               int rel_mis,
                               int C,
                               unsigned int n_cond,
                               int lane,
                               float alpha_sum,
                               float denom,
                               const float* class_nums,
                               const float* sc_row) {
    FILE* f = impute_tsv_file();
    if (!f || !class_nums || !sc_row) return;
    const char* sample = G ? G->name.c_str() : "?";
    std::fprintf(f,
                 "IMPUTE_TSV\tsample=%s\tmode=%s\tlocus=%d\tbio_locus=%d\tss_idx=%d\trel_mis=%d\tC=%d\tn_cond=%u\tlane=%d\talpha_sum=%.8e\tdenom=%.8e\tnum=[",
                 sample, mode, locus, bio_locus, ss_idx, rel_mis, C, n_cond, lane, alpha_sum, denom);
    for (int c = 0; c < C; ++c) {
        if (c) std::fputc(',', f);
        std::fprintf(f, "%.8e", class_nums[c]);
    }
    std::fprintf(f, "]\tsc=[");
    for (int c = 0; c < C; ++c) {
        if (c) std::fputc(',', f);
        std::fprintf(f, "%.8e", sc_row[c]);
    }
    std::fprintf(f, "]\n");
    std::fflush(f);
}

class haplotype_segment_single {
private:
    //EXTERNAL DATA
    hmm_parameters & M;
    genotype * G;
    bitmatrix Hhap, Hvar;

    //COORDINATES & CONSTANTS
    int segment_first;
    int segment_last;
    int locus_first;
    int locus_last;
    int ambiguous_first;
    int ambiguous_last;
    int missing_first;
    int missing_last;
    int transition_first;
    int transition_last;
    unsigned int n_cond_haps;
    unsigned int n_missing;

    //CURSORS
    // === SUPERSITE INDEXING SEMANTICS ===
    // PER-VARIANT cursors: Include all biallelic splits (anchor + siblings)
    // PER-BASEPAIR cursors: Only biological positions (anchors=1, siblings excluded)
    int curr_segment_index;
    int curr_segment_locus;     // PER-BASEPAIR: biological position within segment
    int curr_abs_locus;         // PER-VARIANT: absolute variant index [0, n_variants)
    int prev_abs_locus;
    int curr_rel_locus;
    int curr_rel_locus_offset;
    int curr_abs_ambiguous;     // PER-BASEPAIR: index into G->Ambiguous[], siblings excluded
                                // Controlled by data_amb = hmm_amb && !is_sibling
    int curr_abs_transition;
    int curr_abs_missing;
    int curr_rel_missing;


    //DYNAMIC ARRAYS
    float probSumT;
    aligned_vector32 < float > prob;
    aligned_vector32 < float > probSumK;
    // Store both donor (row) and lane (column) marginals for outer-product seeding at segment boundaries
    aligned_vector32 < float > probSumH;
    struct alignas(32) LaneMarginal {
        float lane[HAP_NUMBER];
        LaneMarginal() : lane{0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f} {}
    };
    std::vector < aligned_vector32 < float > > Alpha;
    std::vector < aligned_vector32 < float > > AlphaSum;
    std::vector<LaneMarginal> AlphaLaneSum;
    std::vector < int > AlphaLocus;
    aligned_vector32 < float > AlphaSumSum;
    std::vector < aligned_vector32 < float > > AlphaMissing;
    std::vector < aligned_vector32 < float > > AlphaSumMissing;
    float HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
    double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));

    //STATIC ARRAYS
    float sumHProbs;
    double sumDProbs;
    float g0[HAP_NUMBER], g1[HAP_NUMBER];
    float nt, yt;

    // Lane weights helpers for first-transition seeding
    struct LaneWeights {
        float match[HAP_NUMBER];
        float neither;
    };
    LaneWeights compute_lane_match_weights(const SiteView& site_view) const;
    void build_lane_priors_first(const SiteView& site_view, double lane_probs[HAP_NUMBER], bool use_outer_prod) const;

    //SUPER-SITE SUPPORT
    const std::vector<SuperSite>* super_sites;
    const std::vector<bool>* is_super_site;
    const std::vector<int>* locus_to_super_idx;
    const uint8_t* panel_codes;
    size_t panel_codes_size;
    const std::vector<int>* super_site_var_index;
    const std::vector<unsigned int>* cond_idx;
    const std::vector<uint32_t>* supersite_sc_offset;  // Thread-local SC offsets (set during backward)

    // Static panel code matrix: [supersite_idx][cond_hap_relative_idx] -> panel_code
    // Populated once during initialization using cond_idx snapshot
    // Similar to Hvar for biallelic, avoids dynamic cond_idx dependency
    std::vector<aligned_vector32<uint8_t>> ss_panel_matrix;

    aligned_vector32<uint8_t> ss_cond_codes;  // Temporary workspace for current supersite
    aligned_vector32<float> ss_emissions;
    aligned_vector32<float> ss_emissions_h1;
    const bool supersites_enabled_flag;

    // Anchor MIS mapping: record which relative-missing index belongs to a given locus
    std::vector<int> missing_index_by_locus; // size = locus_last - locus_first + 1, init -1

    // EMISSION HELPERS
    bool prepare_outer_product_mix(int rel_prev_segment, __m256& col_mix, float& row_stay, float& row_switch, bool allow_outer = true);
    void trace_ambiguous_cursor(const char* stage, int locus, bool is_sibling, int expected_delta) const;
    void trace_log_forward_state(int locus,
                                 int prev_before,
                                 int prev_after,
                                 double yt_val,
                                 double nt_val,
                                 bool update_prev,
                                 bool is_anchor,
                                 bool is_sibling,
                                 bool hmm_amb,
                                 bool hmm_mis,
                                 bool hmm_hom);
    void trace_log_backward_state(int locus,
                                  int prev_before,
                                  int prev_after,
                                  double yt_val,
                                  double nt_val,
                                  bool update_prev,
                                  bool is_anchor,
                                  bool is_sibling,
                                  bool hmm_amb,
                                  bool hmm_mis,
                                  bool hmm_hom);
    mutable int trace_forward_pre_cursor;
    mutable int trace_forward_pre_locus;
    mutable bool trace_forward_pre_valid;
    mutable int trace_backward_pre_cursor;
    mutable int trace_backward_pre_locus;
    mutable bool trace_backward_pre_valid;
    mutable bool trace_forward_active;
    mutable bool trace_backward_active;
    bool trace_forward_table_header_emitted;
    bool trace_backward_table_header_emitted;

    // Trace noise guard: limit anchor emission logs per window
    int trace_anchor_match_logs_remaining = 0;

    //INLINED AND UNROLLED ROUTINES
    void INIT_HOM();
    void INIT_AMB();
    void INIT_MIS();
    // Sibling-specific no-op DP helpers (bookkeeping only)
    void INIT_SIB(const SiteView& site_view);
    void RUN_SIB(const SiteView& site_view);
    void COLLAPSE_SIB(const SiteView& site_view);
    void handle_sibling_bookkeeping(const SiteView& site_view);
    bool RUN_HOM(char);
    void RUN_AMB();
    void RUN_MIS();
    void COLLAPSE_HOM();
    void COLLAPSE_AMB();
    void COLLAPSE_MIS();
    
    //SUPERSITE HELPER ROUTINES
    void SS_INIT_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code);
    void SS_INIT_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1);
    void SS_INIT_MIS();
    bool SS_RUN_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code);
    bool SS_RUN_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1);
    bool SS_RUN_MIS();
    void SS_COLLAPSE_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code);
    void SS_COLLAPSE_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1);
    void SS_COLLAPSE_MIS();
    
    void ss_load_cond_codes(const SuperSite& ss, int ss_idx); // Cache helper
    
    void SUMK();
    void IMPUTE(std::vector < float > & );
    void IMPUTE_SUPERSITE_MULTIVARIATE(std::vector<float>& SC, const SuperSite& ss, int ss_idx, int rel_missing_index, const std::vector<uint32_t>* supersite_sc_offset);
    bool TRANS_HAP();
    bool TRANS_DIP_MULT();
    bool TRANS_DIP_ADD();
    void SET_FIRST_TRANS(std::vector < double > & );
    int SET_OTHER_TRANS(std::vector < double > & );

public:
    //CONSTRUCTOR/DESTRUCTOR
    haplotype_segment_single(genotype *, bitmatrix &, std::vector < unsigned int > &, window &, hmm_parameters &,
        const std::vector<SuperSite>* _super_sites = nullptr,
        const std::vector<bool>* _is_super_site = nullptr,
        const std::vector<int>* _locus_to_super_idx = nullptr,
        const uint8_t* _panel_codes = nullptr,
        size_t _panel_codes_size = 0,
        const std::vector<int>* _super_site_var_index = nullptr);
    ~haplotype_segment_single();

    //void fetch();
    void forward();
    int backward(std::vector < double > &, std::vector < float > &, 
                 std::vector<float>* SC = nullptr, 
                 const std::vector<bool>* anchor_has_missing = nullptr,
                 const std::vector<uint32_t>* supersite_sc_offset = nullptr);
    
};

/*******************************************************************************/
/*****************         SUPERSITE HELPER FUNCTIONS      ********************/
/*******************************************************************************/

/*
 * EMISSION PATTERN DESIGN RATIONALE (BUG #6 DOCUMENTATION):
 * 
 * Biallelic and supersite code use different emission computation patterns,
 * but both implement identical Li & Stephens emission probabilities:
 *   emit[h] = (donor_matches_sample[h]) ? 1.0 : (ed/ee)
 * 
 * BIALLELIC PATTERN (optimized for binary alleles):
 *   Uses inline conditional multiplication:
 *     if (ag != ah) _prob *= mismatch;
 *   Advantages:
 *     - Minimal overhead for simple boolean comparison (0 vs 1)
 *     - Direct scalar-vector multiplication
 *     - Cache-friendly for dense biallelic data
 * 
 * SUPERSITE PATTERN (required for per-lane multi-allelic semantics):
 *   Uses precomputed emission vectors or SIMD blend:
 *     emit = _mm256_blendv_ps(mis_f, match_f, match_mask);
 *   Advantages:
 *     - Handles per-lane expected class (each lane may expect different ALT)
 *     - Supports 4-bit allele codes (0-15) vs binary (0-1)
 *     - AMB sites need different emissions per lane based on amb_code mask
 *   Requirements:
 *     - Must build per-lane expected class array from amb_code
 *     - Must compare donor 4-bit code against per-lane expected codes
 *     - Cannot use simple scalar broadcast (each lane has different expectation)
 * 
 * MATHEMATICAL EQUIVALENCE:
 *   Biallelic: if (donor_allele != sample_allele) → multiply by (ed/ee)
 *   Supersite: if (donor_code != expected_class[lane]) → multiply by (ed/ee)
 *   Both: P(observe sample | donor state) = match ? 1.0 : error_rate
 * 
 * DECISION: Keep both patterns (each optimized for its use case).
 * The divergence is intentional and necessary, not a bug to fix.
 */

// Cache donor codes for a supersite to avoid repeated unpacking
inline
void haplotype_segment_single::ss_load_cond_codes(const SuperSite& ss, int ss_idx) {
    // NEW APPROACH: Use static panel matrix populated during initialization
    // This avoids dynamic dependency on cond_idx, fixing the forward/backward divergence bug

    if (ss_idx < 0 || ss_idx >= (int)ss_panel_matrix.size()) {
        std::fprintf(stderr, "ERROR: ss_load_cond_codes called with invalid ss_idx=%d (matrix size=%zu)\n",
                     ss_idx, ss_panel_matrix.size());
        std::abort();
    }

    // LOAD_TRACE: Log when codes are loaded and what we're reading
    const char* load_trace = std::getenv("SHAPEIT5_LOAD_TRACE");
    bool trace_this_call = false;
    if (load_trace && load_trace[0] != '\0' && load_trace[0] != '0' &&
        (ss_idx == 1 || ss_idx == 3 || ss_idx == 6)) {
        trace_this_call = true;
        std::fprintf(stderr, "[LOAD_TRACE] ss_load_cond_codes called: ss_idx=%d locus=%d\n",
                     ss_idx, curr_abs_locus);
        std::fprintf(stderr, "[LOAD_TRACE] ss_panel_matrix[%d] address=%p\n",
                     ss_idx, (void*)&ss_panel_matrix[ss_idx]);
        std::fprintf(stderr, "[LOAD_TRACE] Matrix panel codes (first 16): ");
        for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
            std::fprintf(stderr, "%u ", (unsigned)ss_panel_matrix[ss_idx][k]);
        }
        std::fprintf(stderr, "\n");

        // Cross-check packed codes against Hvar for the current anchor
        const bool anchor_in_window = (curr_abs_locus >= locus_first) && (curr_abs_locus <= locus_last);
        const int rel_row = anchor_in_window ? (curr_abs_locus - locus_first + curr_rel_locus_offset) : -1;
        std::fprintf(stderr,
                     "[LOAD_TRACE] anchor_in_window=%d rel_row=%d locus_first=%d locus_last=%d offset=%d\n",
                     anchor_in_window ? 1 : 0, rel_row, locus_first, locus_last, curr_rel_locus_offset);
        const unsigned int dump_limit = std::min(32u, n_cond_haps);
        for (unsigned int k = 0; k < dump_limit; ++k) {
            const unsigned int gh = (*cond_idx)[k];
            const uint8_t packed_code = ss_panel_matrix[ss_idx][k];
            const int hvar = anchor_in_window ? static_cast<int>(Hvar.get(rel_row, k)) : -1;
            std::fprintf(stderr,
                         "[LOAD_CHK] ss_idx=%d locus=%d k=%u cond_idx=%u packed=%u hvar=%d\n",
                         ss_idx, curr_abs_locus, k, gh, (unsigned)packed_code, hvar);
        }
    }

    // Simply copy from static matrix using relative indices (like biallelic Hvar)
    ss_cond_codes.resize(n_cond_haps);
    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        ss_cond_codes[k] = ss_panel_matrix[ss_idx][k];
    }

    // Verify what was loaded
    if (trace_this_call) {
        std::fprintf(stderr, "[LOAD_TRACE] Loaded ss_cond_codes (first 16): ");
        for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
            std::fprintf(stderr, "%u ", (unsigned)ss_cond_codes[k]);
        }
        std::fprintf(stderr, "\n");
    }

    // The static matrix is our "cache" and it's always valid.
}

inline
void haplotype_segment_single::SS_INIT_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code) {
    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);
    
    // Precompute emissions
    precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, M.error_ratio[curr_abs_locus], ss_emissions);
    
    // Initialize probabilities
    __m256 _sum = _mm256_set1_ps(0.0f);
 for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
    __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
    __m256 _prob = _emit;
    _sum = _mm256_add_ps(_sum, _prob);
    _mm256_store_ps(&prob[i], _prob);
}
_mm256_store_ps(&probSumH[0], _sum);
probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

 // Trace: matches per lane at anchor (HOM)
 if (supersite_trace_enabled() && trace_anchor_match_logs_remaining > 0) {
     int lane_matches[HAP_NUMBER] = {0};
     for (unsigned int k = 0; k < n_cond_haps; ++k) {
         const uint8_t dc = ss_cond_codes[k];
         for (int h = 0; h < HAP_NUMBER; ++h) {
             if (dc == sample_code) lane_matches[h]++;
         }
     }
     supersite_trace_log("[SupersiteEmit] matches_per_lane locus=%d ss_idx=%d n_cond=%u",
                         curr_abs_locus, ss_idx, n_cond_haps);
     for (int h = 0; h < HAP_NUMBER; ++h) supersite_trace_log(" %d:%d", h, lane_matches[h]);
     supersite_trace_log("\n");
     trace_anchor_match_logs_remaining--;
 }
}

inline
void haplotype_segment_single::SS_INIT_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1) {
    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);

    // Build per-lane expected class vector (8 lanes)
    uint8_t expected_class[HAP_NUMBER];
    unsigned char amb_mask = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                            ? G->Ambiguous[curr_abs_ambiguous] : 0u;
    if (c0 == c1) {
        for (int h = 0; h < HAP_NUMBER; ++h) expected_class[h] = c0;
    } else {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            bool use_c1 = ((amb_mask >> h) & 1U);
            expected_class[h] = use_c1 ? c1 : c0;
        }
    }
    __m256 match_f = _mm256_set1_ps(1.0f);
    __m256 mis_f   = _mm256_set1_ps(M.error_ratio[curr_abs_locus]);
    alignas(32) int expv[HAP_NUMBER];
    for (int h = 0; h < HAP_NUMBER; ++h) expv[h] = (int)expected_class[h];
    __m256i exp_vec = _mm256_load_si256((__m256i*)expv);

    // Determine anchor ALT code (1..n_alts) for current anchor locus
    int anchor_code = 0; 
    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
        if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { anchor_code = (int)ai + 1; break; }
    }

    // Strict multi-allelic semantics: donor must match per-lane expected class
    __m256 _sum = _mm256_set1_ps(0.0f);
    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256i donor_vec = _mm256_set1_epi32((int)ss_cond_codes[k]);
        __m256i match_vec_i = _mm256_cmpeq_epi32(donor_vec, exp_vec);
        __m256 emit = _mm256_blendv_ps(mis_f, match_f, _mm256_castsi256_ps(match_vec_i));
        _sum = _mm256_add_ps(_sum, emit);
        _mm256_store_ps(&prob[i], emit);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

    // Trace: matches per lane at anchor (AMB)
    if (supersite_trace_enabled() && trace_anchor_match_logs_remaining > 0) {
        int lane_matches[HAP_NUMBER] = {0};
        for (unsigned int k = 0; k < n_cond_haps; ++k) {
            const uint8_t dc = ss_cond_codes[k];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                if ((int)dc == (int)expected_class[h]) lane_matches[h]++;
            }
        }
        supersite_trace_log("[SupersiteEmit] matches_per_lane locus=%d ss_idx=%d n_cond=%u",
                            curr_abs_locus, ss_idx, n_cond_haps);
        for (int h = 0; h < HAP_NUMBER; ++h) supersite_trace_log(" %d:%d", h, lane_matches[h]);
        supersite_trace_log("\n");
        trace_anchor_match_logs_remaining--;
    }
}

inline
void haplotype_segment_single::SS_INIT_MIS() {
    // Missing data: uniform probabilities
    fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
    fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
    probSumT = 1.0f;
}

inline
bool haplotype_segment_single::SS_RUN_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code) {
    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);

    // DETAILED BACKWARD TRACING
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    bool trace_enabled = (tr && tr[0] != '\0' && tr[0] != '0');
    if (trace_enabled) {
        std::fprintf(stdout, "[SS_RUN_HOM_TRACE] locus=%d ss_idx=%d sample_code=%u n_cond=%u\n",
                     curr_abs_locus, ss_idx, (unsigned)sample_code, n_cond_haps);
        std::fprintf(stdout, "  Input probSumT=%.10f\n", (double)probSumT);
        std::fprintf(stdout, "  Input probSumH=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                     (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3],
                     (double)probSumH[4], (double)probSumH[5], (double)probSumH[6], (double)probSumH[7]);
        std::fprintf(stdout, "  Transition: yt=%.10f nt=%.10f\n", yt, nt);
        std::fprintf(stdout, "  Panel codes: ");
        for (unsigned int k = 0; k < n_cond_haps; ++k) {
            std::fprintf(stdout, "%u ", (unsigned)ss_cond_codes[k]);
        }
        std::fprintf(stdout, "\n  Match counts: ");
        unsigned int matches = 0, mismatches = 0;
        for (unsigned int k = 0; k < n_cond_haps; ++k) {
            if (ss_cond_codes[k] == sample_code) matches++; else mismatches++;
        }
        std::fprintf(stdout, "match=%u mismatch=%u\n", matches, mismatches);
        // Emit per-lane classes (all lanes want sample_code for HOM)
        std::fprintf(stdout, "[SS_RUN_HOM_CLASSES] locus=%d ss_idx=%d class=%u lanes:", curr_abs_locus, ss_idx, (unsigned)sample_code);
        for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %u", (unsigned)sample_code);
        std::fprintf(stdout, "\n");
    }

    // Precompute emissions
    precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, M.error_ratio[curr_abs_locus], ss_emissions);

    if (trace_enabled) {
        std::fprintf(stdout, "  Panel codes and emissions:\n");
        std::fprintf(stdout, "  SS_DETAILED: _factor=%.10f _nt=%.10f\n",
                     yt / (n_cond_haps * probSumT), nt / probSumT);
        alignas(32) float tfreq_arr[HAP_NUMBER];
        _mm256_store_ps(tfreq_arr, _mm256_load_ps(&probSumH[0]));
        std::fprintf(stdout, "  SS_DETAILED: probSumH_before=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                     tfreq_arr[0], tfreq_arr[1], tfreq_arr[2], tfreq_arr[3],
                     tfreq_arr[4], tfreq_arr[5], tfreq_arr[6], tfreq_arr[7]);
    }

    // Update probabilities with transitions
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);

    if (trace_enabled) {
        alignas(32) float tfreq_arr2[HAP_NUMBER];
        _mm256_store_ps(tfreq_arr2, _tFreq);
        std::fprintf(stdout, "  SS_DETAILED: _tFreq=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                     tfreq_arr2[0], tfreq_arr2[1], tfreq_arr2[2], tfreq_arr2[3],
                     tfreq_arr2[4], tfreq_arr2[5], tfreq_arr2[6], tfreq_arr2[7]);
    }

    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
        __m256 _prob = _mm256_load_ps(&prob[i]);
        if (trace_enabled && k == 0) {
            alignas(32) float prob_before[HAP_NUMBER];
            _mm256_store_ps(prob_before, _prob);
            std::fprintf(stdout, "  SS_DETAILED k=0 prob_before=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                         prob_before[0], prob_before[1], prob_before[2], prob_before[3],
                         prob_before[4], prob_before[5], prob_before[6], prob_before[7]);
        }
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        if (trace_enabled && k == 0) {
            alignas(32) float prob_after_fma[HAP_NUMBER];
            _mm256_store_ps(prob_after_fma, _prob);
            std::fprintf(stdout, "  SS_DETAILED k=0 prob_after_FMA=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                         prob_after_fma[0], prob_after_fma[1], prob_after_fma[2], prob_after_fma[3],
                         prob_after_fma[4], prob_after_fma[5], prob_after_fma[6], prob_after_fma[7]);
        }
        if (ss_emissions[k] != 1.0f) _prob = _mm256_mul_ps(_prob, _emit);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);

        if (trace_enabled) {
            std::fprintf(stdout, "    k=%u code=%u match=%d emission=%.6f\n",
                         k, (unsigned)ss_cond_codes[k], (ss_cond_codes[k] == sample_code), ss_emissions[k]);
        }
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

    // DETAILED BACKWARD TRACING - OUTPUT
    if (trace_enabled) {
        std::fprintf(stdout, "  Output probSumT=%.10f\n", (double)probSumT);
        std::fprintf(stdout, "  Output probSumH=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                     (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3],
                     (double)probSumH[4], (double)probSumH[5], (double)probSumH[6], (double)probSumH[7]);
    }

    return true;
}

inline
bool haplotype_segment_single::SS_RUN_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1) {
    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);

    // DEBUG: Check if we're in trace mode
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    bool trace_enabled = (tr && tr[0] != '\0' && tr[0] != '0');

    // Build per-lane expected class vector
    uint8_t expected_class[HAP_NUMBER];
    unsigned char amb_mask = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                            ? G->Ambiguous[curr_abs_ambiguous] : 0u;
    if (c0 == c1) {
        for (int h = 0; h < HAP_NUMBER; ++h) expected_class[h] = c0;
    } else {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            bool use_c1 = ((amb_mask >> h) & 1U);
            expected_class[h] = use_c1 ? c1 : c0;
        }
    }

    if (trace_enabled) {
        std::fprintf(stderr, "DEBUG SS_RUN_AMB: trace_enabled is TRUE. locus=%d\n", curr_abs_locus);
    }

	if (trace_enabled) {
		std::fprintf(stderr, "[SS_LANE_DEBUG] amb_mask=%d (0x%02x) c0=%d c1=%d\n", (int)amb_mask, (int)amb_mask, (int)c0, (int)c1);
		for (int h = 0; h < HAP_NUMBER; ++h) {
			bool bit = ((amb_mask >> h) & 1U);
			std::fprintf(stderr, "  Lane %d: bit=%d expected=%d\n", h, (int)bit, (int)expected_class[h]);
		}
		// Emit per-lane expected classes to stdout for easier capture
		std::fprintf(stdout, "[SS_RUN_AMB_CLASSES] locus=%d ss_idx=%d c0=%u c1=%u amb_mask=0x%02x lanes:", curr_abs_locus, ss_idx, (unsigned)c0, (unsigned)c1, (unsigned)amb_mask);
		for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %u", (unsigned)expected_class[h]);
		std::fprintf(stdout, "\n");
		// Emit first few donor codes and match mask
		std::fprintf(stdout, "[SS_RUN_AMB_CODES] locus=%d ss_idx=%d n_cond=%u codes:", curr_abs_locus, ss_idx, (unsigned)n_cond_haps);
		for (unsigned int k = 0; k < std::min(16u, n_cond_haps); ++k) {
			std::fprintf(stdout, " %u", (unsigned)ss_cond_codes[k]);
		}
		std::fprintf(stdout, "\n");
	}

    // HYPOTHESIS 1 DEBUGGING
    if (!debug::SUPERDEBUG_SAMPLENAME.empty() && G->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
        std::cout << "[SUPERDEBUG] H1: ENTER SS_RUN_AMB" << std::endl;
        std::cout << "[SUPERDEBUG] H1: Sample=" << G->name << " Pos=" << ss.global_site_id << " c0=" << (int)c0 << " c1=" << (int)c1 << " amb_mask=" << (int)amb_mask << std::endl;
        std::string expected_str = "  Expected classes: ";
        for (int h=0; h<HAP_NUMBER; ++h) expected_str += std::to_string(expected_class[h]) + " ";
        std::cout << expected_str << std::endl;
    }
    
    // DEBUG: Print emission setup
    if (trace_enabled) {
        std::fprintf(stdout, "SS_RUN_AMB: locus=%d c0=%u c1=%u amb_mask=0x%02x\n", 
                     curr_abs_locus, (unsigned)c0, (unsigned)c1, (unsigned)amb_mask);
        std::fprintf(stdout, "  Expected classes per lane:");
        for (int h = 0; h < HAP_NUMBER; ++h) {
            std::fprintf(stdout, " %u", (unsigned)expected_class[h]);
        }
        std::fprintf(stdout, "\n");
        std::fprintf(stdout, "  Panel codes:");
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            std::fprintf(stdout, " %u", (unsigned)ss_cond_codes[k]);
        }
        std::fprintf(stdout, "\n");
    }
    
    alignas(32) int expv[HAP_NUMBER];
    for (int h = 0; h < HAP_NUMBER; ++h) expv[h] = (int)expected_class[h];
    __m256i exp_vec = _mm256_load_si256((__m256i*)expv);
    __m256 match_f = _mm256_set1_ps(1.0f);
    __m256 mis_f   = _mm256_set1_ps(M.error_ratio[curr_abs_locus]);

    // Update probabilities with transitions and per-lane emissions
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);

    // Strict multi-allelic semantics: donor must match per-lane expected class
    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256 _prob = _mm256_load_ps(&prob[i]);
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);

        __m256i donor_vec = _mm256_set1_epi32((int)ss_cond_codes[k]);
        __m256i match_vec_i = _mm256_cmpeq_epi32(donor_vec, exp_vec);
        __m256 emit = _mm256_blendv_ps(mis_f, match_f, _mm256_castsi256_ps(match_vec_i));

        if (!debug::SUPERDEBUG_SAMPLENAME.empty() && G->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP && k < 5) {
            std::cout << "[SUPERDEBUG] H1: k=" << k << " donor_code=" << (int)ss_cond_codes[k] << " (strict) first-lane emit=" << ((float*)&emit)[0] << std::endl;
        }

        if (trace_enabled) {
            std::fprintf(stdout, "    SS_donor_%d dc=%d strict per-lane emit:", k, (int)ss_cond_codes[k]);
            alignas(32) float emit_arr_dbg[HAP_NUMBER];
            _mm256_store_ps(emit_arr_dbg, emit);
            for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", emit_arr_dbg[h]);
            std::fprintf(stdout, "\n");
        }

        _prob = _mm256_mul_ps(_prob, emit);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
    return true;
}

inline
bool haplotype_segment_single::SS_RUN_MIS() {
    // Missing data: only apply transitions, no emissions
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);
    
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        __m256 _prob = _mm256_load_ps(&prob[i]);
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
    return true;
}

inline
void haplotype_segment_single::SS_COLLAPSE_HOM(const SuperSite& ss, int ss_idx, uint8_t sample_code) {
    // DEBUG: Log entry state
    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[SS_COLLAPSE_HOM_ENTER] locus=%d ss_idx=%d probSumT_before=%.10f yt=%.10f nt=%.10f n_cond=%u\n",
                     curr_abs_locus, ss_idx, (double)probSumT, (double)yt, (double)nt, n_cond_haps);
    }

    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);

    // Precompute emissions
    precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, M.error_ratio[curr_abs_locus], ss_emissions);

    // Collapse from probSumK with optional outer-product seeding
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);

    // DEBUG: Log flags and computed values
    if (supersite_trace_enabled()) {
        float nt_val = use_outer ? 0.0f : (nt / probSumT);
        float tFreq_val = use_outer ? 0.0f : (yt / n_cond_haps);
        std::fprintf(stderr, "[SS_COLLAPSE_HOM_FLAGS] use_outer=%d rel_prev_seg=%d nt_val=%.10f tFreq_val=%.10f\n",
                     (int)use_outer, rel_prev_seg, (double)nt_val, (double)tFreq_val);
    }
    
    // DEBUG: Log first few probSumK values
    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[SS_COLLAPSE_HOM_PROBSUMK] first_4:");
        for (int k = 0; k < std::min(4, (int)n_cond_haps); ++k) {
            std::fprintf(stderr, " %.10f", (double)probSumK[k]);
        }
        std::fprintf(stderr, "\n");
    }

    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }
        if (ss_emissions[k] != 1.0f) _prob = _mm256_mul_ps(_prob, _emit);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

    // DEBUG: Log exit state
    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[SS_COLLAPSE_HOM_EXIT] locus=%d probSumT_after=%.10f probSumH[0-3]=%.6f %.6f %.6f %.6f\n",
                     curr_abs_locus, (double)probSumT, (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3]);
    }

    // Trace: matches per lane at anchor (HOM)
    if (supersite_trace_enabled() && trace_anchor_match_logs_remaining > 0) {
        int lane_matches[HAP_NUMBER] = {0};
        for (unsigned int k = 0; k < n_cond_haps; ++k) {
            const uint8_t dc = ss_cond_codes[k];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                if (dc == sample_code) lane_matches[h]++;
            }
        }
        supersite_trace_log("[SupersiteEmit] matches_per_lane locus=%d ss_idx=%d n_cond=%u",
                            curr_abs_locus, ss_idx, n_cond_haps);
        for (int h = 0; h < HAP_NUMBER; ++h) supersite_trace_log(" %d:%d", h, lane_matches[h]);
        supersite_trace_log("\n");
        trace_anchor_match_logs_remaining--;
    }
}

inline
void haplotype_segment_single::SS_COLLAPSE_AMB(const SuperSite& ss, int ss_idx, uint8_t c0, uint8_t c1) {
    // Load cached conditioning haplotype codes
    ss_load_cond_codes(ss, ss_idx);

    // DEBUG: Check if we're in trace mode
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    bool trace_enabled = (tr && tr[0] != '\0' && tr[0] != '0');

    // Build per-lane expected class vector
    uint8_t expected_class[HAP_NUMBER];
    unsigned char amb_mask = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                            ? G->Ambiguous[curr_abs_ambiguous] : 0u;
    if (c0 == c1) {
        for (int h = 0; h < HAP_NUMBER; ++h) expected_class[h] = c0;
    } else {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            bool use_c1 = ((amb_mask >> h) & 1U);
            expected_class[h] = use_c1 ? c1 : c0;
        }
    }
    
    // DEBUG: Print emission setup
    if (trace_enabled) {
        std::fprintf(stdout, "SS_COLLAPSE_AMB: locus=%d c0=%u c1=%u amb_mask=0x%02x\n", 
                     curr_abs_locus, (unsigned)c0, (unsigned)c1, (unsigned)amb_mask);
        std::fprintf(stdout, "  Expected classes per lane:");
        for (int h = 0; h < HAP_NUMBER; ++h) {
            std::fprintf(stdout, " %u", (unsigned)expected_class[h]);
        }
        std::fprintf(stdout, "\n");
    }
    
    // Precompute emission vectors based on mode (BUG FIX #5: unified code path)
    __m256 match_f = _mm256_set1_ps(1.0f);
    __m256 mis_f   = _mm256_set1_ps(M.error_ratio[curr_abs_locus]);
    __m256i exp_vec;
    
    // Strict 4-bit class equality semantics
    alignas(32) int expv[HAP_NUMBER];
    for (int h = 0; h < HAP_NUMBER; ++h)
        expv[h] = (int)expected_class[h];
    exp_vec = _mm256_load_si256((__m256i*)expv);

    // Collapse from previous segment with per-lane emissions
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);

    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }

        __m256i donor_vec = _mm256_set1_epi32((int)ss_cond_codes[k]);
        __m256i match_vec_i = _mm256_cmpeq_epi32(donor_vec, exp_vec);
        __m256 emit = _mm256_blendv_ps(mis_f, match_f, _mm256_castsi256_ps(match_vec_i));
        if (trace_enabled) {
            alignas(32) float emit_arr_dbg[HAP_NUMBER];
            _mm256_store_ps(emit_arr_dbg, emit);
            std::fprintf(stdout, "  Donor %d: raw_code=%u strict emit:", k, (unsigned)ss_cond_codes[k]);
            for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", emit_arr_dbg[h]);
            std::fprintf(stdout, "\n");
        }
        _prob = _mm256_mul_ps(_prob, emit);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
    if (trace_enabled) {
        std::fprintf(stdout, "  BIAL sums lanes:");
        for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", probSumH[h]);
        std::fprintf(stdout, " probSumT=%.6f\n", probSumT);
    }
}

inline
void haplotype_segment_single::SS_COLLAPSE_MIS() {
    // Missing data: only apply transitions, no emissions
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);
    
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

}

/*******************************************************************************/
/*****************            HOMOZYGOUS GENOTYPE            ************************/
/*******************************************************************************/


inline
bool haplotype_segment_single::prepare_outer_product_mix(int rel_prev_segment, __m256& col_mix, float& row_stay, float& row_switch, bool allow_outer) {
    // Outer product is DISABLED by default due to divergence bugs
    // Set SHAPEIT5_ENABLE_OUTER_PRODUCT=1 to enable (for testing/comparison)
    static int enable_outer = -1;
    if (enable_outer < 0) {
        const char* env = std::getenv("SHAPEIT5_ENABLE_OUTER_PRODUCT");
        enable_outer = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
    }
    if (!enable_outer) return false;

    if (!allow_outer) return false;
    if (!supersites_enabled_flag) return false;
    if (rel_prev_segment < 0 || rel_prev_segment >= static_cast<int>(AlphaLaneSum.size())) return false;
    if (!n_cond_haps) return false;

    const float prev_total = AlphaSumSum[rel_prev_segment];
    if (prev_total <= std::numeric_limits<float>::min()) return false;

    const __m256 prev_cols = _mm256_load_ps(AlphaLaneSum[rel_prev_segment].lane);
    const __m256 stay_factor = _mm256_set1_ps(nt / prev_total);
    const __m256 switch_vec = _mm256_set1_ps(yt / static_cast<float>(HAP_NUMBER));
    col_mix = _mm256_fmadd_ps(prev_cols, stay_factor, switch_vec);

    row_stay = nt / prev_total;
    row_switch = yt / static_cast<float>(n_cond_haps);
    return true;
}

inline
void haplotype_segment_single::INIT_HOM() {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling at window boundary: initialize neutrally to avoid underflow
            INIT_MIS();
            return;
        }

        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;
            case SSClass::HOM: SS_INIT_HOM(ss, ss_idx, c0); return;
            case SSClass::AMB: SS_INIT_AMB(ss, ss_idx, c0, c1); return;
        }
    }

    // Biallelic path (original simple loop)
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _mm256_set1_ps((ag==ah)?1.0f:M.error_ratio[curr_abs_locus]);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
bool haplotype_segment_single::RUN_HOM(char rare_allele) {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: true no-op (preserve probability state, no DP)
            return true;
        }
        
        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return true;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: return SS_RUN_HOM(ss, ss_idx, c0);
            case SSClass::AMB: return SS_RUN_AMB(ss, ss_idx, c0, c1);
        }
    }
    
    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    if (rare_allele < 0 || ag == rare_allele) {
        // EMISSION TRACING
        const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
        bool trace_enabled = (tr && tr[0] != '\0' && tr[0] != '0');
        if (trace_enabled) {
            std::fprintf(stdout, "[BIAL_RUN_HOM_TRACE] locus=%d sample_allele=%u n_cond=%u\n",
                         curr_abs_locus, (unsigned)ag, n_cond_haps);
            std::fprintf(stdout, "  Input probSumT=%.10f\n", (double)probSumT);
            std::fprintf(stdout, "  Input probSumH=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                         (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3],
                         (double)probSumH[4], (double)probSumH[5], (double)probSumH[6], (double)probSumH[7]);
            std::fprintf(stdout, "  Transition: yt=%.10f nt=%.10f\n", yt, nt);
        }

        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
        __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
        _tFreq = _mm256_mul_ps(_tFreq, _factor);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        __m256 _mismatch = _mm256_set1_ps(M.error_ratio[curr_abs_locus]);

        if (trace_enabled) {
            std::fprintf(stdout, "  Panel alleles and emissions:\n");
            std::fprintf(stdout, "  DETAILED: _factor=%.10f _nt=%.10f\n",
                         yt / (n_cond_haps * probSumT), nt / probSumT);
            alignas(32) float tfreq_arr[HAP_NUMBER];
            _mm256_store_ps(tfreq_arr, _tFreq);
            std::fprintf(stdout, "  DETAILED: _tFreq=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                         tfreq_arr[0], tfreq_arr[1], tfreq_arr[2], tfreq_arr[3],
                         tfreq_arr[4], tfreq_arr[5], tfreq_arr[6], tfreq_arr[7]);
        }
        for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            if (trace_enabled && k == 0) {
                alignas(32) float prob_before[HAP_NUMBER];
                _mm256_store_ps(prob_before, _prob);
                std::fprintf(stdout, "  DETAILED k=0 prob_before=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                             prob_before[0], prob_before[1], prob_before[2], prob_before[3],
                             prob_before[4], prob_before[5], prob_before[6], prob_before[7]);
            }
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            if (trace_enabled && k == 0) {
                alignas(32) float prob_after_fma[HAP_NUMBER];
                _mm256_store_ps(prob_after_fma, _prob);
                std::fprintf(stdout, "  DETAILED k=0 prob_after_FMA=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                             prob_after_fma[0], prob_after_fma[1], prob_after_fma[2], prob_after_fma[3],
                             prob_after_fma[4], prob_after_fma[5], prob_after_fma[6], prob_after_fma[7]);
            }
            float emission = (ag != ah) ? M.error_ratio[curr_abs_locus] : 1.0f;
            if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);

            if (trace_enabled) {
                std::fprintf(stdout, "    k=%u ah=%u match=%d emission=%.6f\n",
                             k, (unsigned)ah, (ag==ah), emission);
            }
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

        if (trace_enabled) {
            std::fprintf(stdout, "  Output probSumT=%.10f\n", (double)probSumT);
            std::fprintf(stdout, "  Output probSumH=[%.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f, %.10f]\n",
                         (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3],
                         (double)probSumH[4], (double)probSumH[5], (double)probSumH[6], (double)probSumH[7]);
        }
        return true;
    }
    return false;
}

inline
void haplotype_segment_single::COLLAPSE_HOM() {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: true no-op (preserve probability state, no DP)
            return;
        }
        
        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_COLLAPSE_HOM(ss, ss_idx, c0); return;
            case SSClass::AMB: SS_COLLAPSE_AMB(ss, ss_idx, c0, c1); return;
        }
    }
    
    // Biallelic path
    // DEBUG: Log entry state
    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[BIAL_COLLAPSE_HOM_ENTER] locus=%d probSumT_before=%.10f yt=%.10f nt=%.10f n_cond=%u\n",
                     curr_abs_locus, (double)probSumT, (double)yt, (double)nt, n_cond_haps);
    }

    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch, true);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);

    // DEBUG: Log flags and computed values
    if (supersite_trace_enabled()) {
        float nt_val = use_outer ? 0.0f : (nt / probSumT);
        float tFreq_val = use_outer ? 0.0f : (yt / n_cond_haps);
        std::fprintf(stderr, "[BIAL_COLLAPSE_HOM_FLAGS] use_outer=%d rel_prev_seg=%d nt_val=%.10f tFreq_val=%.10f\n",
                     (int)use_outer, rel_prev_seg, (double)nt_val, (double)tFreq_val);
    }
    __m256 _mismatch = _mm256_set1_ps(M.error_ratio[curr_abs_locus]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }
        if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

    // DEBUG: Log exit state
    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[BIAL_COLLAPSE_HOM_EXIT] locus=%d probSumT_after=%.10f probSumH[0-3]=%.6f %.6f %.6f %.6f\n",
                     curr_abs_locus, (double)probSumT, (double)probSumH[0], (double)probSumH[1], (double)probSumH[2], (double)probSumH[3]);
    }

}

/*******************************************************************************/
/*****************            HETEROZYGOUS GENOTYPE            ********************/
/*******************************************************************************/

inline
void haplotype_segment_single::INIT_AMB() {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling at window boundary: initialize neutrally to avoid underflow
            INIT_MIS();  // BUG FIX #1: Use biallelic MIS
            return;
        }

        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_INIT_HOM(ss, ss_idx, c0); return;
            case SSClass::AMB: SS_INIT_AMB(ss, ss_idx, c0, c1); return;
        }
    }
    
    // Biallelic path
    unsigned char amb_code = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                             ? G->Ambiguous[curr_abs_ambiguous] : 0u;

    if (std::getenv("SHAPEIT5_TEST_TRACE")) {
        // Prepare g0, g1 emissions for tracing
        float g0_trace[HAP_NUMBER], g1_trace[HAP_NUMBER];
        for (int h = 0 ; h < HAP_NUMBER ; h ++) {
            g0_trace[h] = HAP_GET(amb_code,h)?M.error_ratio[curr_abs_locus]:1.0f;
            g1_trace[h] = HAP_GET(amb_code,h)?1.0f:M.error_ratio[curr_abs_locus];
        }
        std::fprintf(stdout, "[BIAL_INIT_AMB] locus=%d BIALLELIC PATH amb_code=0x%02x\n", curr_abs_locus, (unsigned)amb_code);
        std::fprintf(stdout, "  Biallelic emission g0: %f %f %f %f %f %f %f %f\n",
                     g0_trace[0], g0_trace[1], g0_trace[2], g0_trace[3],
                     g0_trace[4], g0_trace[5], g0_trace[6], g0_trace[7]);
        std::fprintf(stdout, "  Biallelic emission g1: %f %f %f %f %f %f %f %f\n",
                     g1_trace[0], g1_trace[1], g1_trace[2], g1_trace[3],
                     g1_trace[4], g1_trace[5], g1_trace[6], g1_trace[7]);
    }

    // Biallelic INIT_AMB: simple AVX2 loop (original algorithm)
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.error_ratio[curr_abs_locus]:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.error_ratio[curr_abs_locus];
    }
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _emit[ah];
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_single::RUN_AMB() {
    // DEBUG: Check if we're in trace mode
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    bool trace_enabled = (tr && tr[0] != '\0' && tr[0] != '0');
    
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: true no-op (preserve probability state, no DP)
            if (trace_enabled) {
                std::fprintf(stdout, "BIAL_RUN_AMB: locus=%d SIBLING (no-op)\n", curr_abs_locus);
            }
            return;
        }

        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        if (trace_enabled) {
            std::fprintf(stdout, "BIAL_RUN_AMB: locus=%d SUPERSITE DISPATCH c0=%u c1=%u cls=%d\n", 
                         curr_abs_locus, (unsigned)c0, (unsigned)c1, (int)cls);
        }
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_RUN_HOM(ss, ss_idx, c0); return;
            case SSClass::AMB: SS_RUN_AMB(ss, ss_idx, c0, c1); return;
        }
    }
    
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    if (trace_enabled) {
        std::fprintf(stdout, "BIAL_RUN_AMB: locus=%d BIALLELIC PATH amb_code=0x%02x\n", 
                     curr_abs_locus, (unsigned)amb_code);
    }
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.error_ratio[curr_abs_locus]:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.error_ratio[curr_abs_locus];
    }
    if (trace_enabled) {
        std::fprintf(stderr, "[BIAL_LANE_DEBUG] amb_code=%d (0x%02x)\n", (int)amb_code, (int)amb_code);
        for (int h = 0; h < HAP_NUMBER; ++h) {
             bool bit = HAP_GET(amb_code, h);
             std::fprintf(stderr, "  Lane %d: bit=%d g0=%.4f g1=%.4f\n", h, (int)bit, g0[h], g1[h]);
        }
    }
    if (supersite_trace_enabled()) {
        std::fprintf(stdout, "[BIAL_AMB_TRACE] loc=%d amb_code=0x%02x\n", curr_abs_locus, (unsigned)amb_code);
        std::fprintf(stdout, "  g0:"); for(int h=0; h<HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", g0[h]); std::fprintf(stdout, "\n");
        std::fprintf(stdout, "  g1:"); for(int h=0; h<HAP_NUMBER; ++h) std::fprintf(stdout, " %.6f", g1[h]); std::fprintf(stdout, "\n");
    }
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _mm256_load_ps(&prob[i]);
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        _prob = _mm256_mul_ps(_prob, _emit[ah]);
        
        // DEBUG: Log detailed emission calculations for all donors
        if (trace_enabled) {
            alignas(32) float emit_vals[HAP_NUMBER];
            _mm256_store_ps(emit_vals, _emit[ah]);
            std::fprintf(stdout, "    Bial_donor_%d ah=%d per-lane emit:", k, (int)ah);
            for (int h = 0; h < HAP_NUMBER; ++h) {
                std::fprintf(stdout, " %.6f", emit_vals[h]);
            }
            std::fprintf(stdout, "\n");
        }
        
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*
inline
void haplotype_segment_single::RUN_AMB() {
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.error_ratio[curr_abs_locus]:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.error_ratio[curr_abs_locus];
    }
    __m256d _sum0 = _mm256_set1_pd(0.0f);
    __m256d _sum1 = _mm256_set1_pd(0.0f);
    __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
    __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
    __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
    _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
    _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    __m256d _emit0[2], _emit1[2];
    _emit0[0] = _mm256_loadu_pd(&g0[0]);
    _emit0[1] = _mm256_loadu_pd(&g1[0]);
    _emit1[0] = _mm256_loadu_pd(&g0[4]);
    _emit1[1] = _mm256_loadu_pd(&g1[4]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0 = _mm256_load_pd(&prob[i+0]);
        __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
        _prob0 = _mm256_mul_pd(_prob0, _emit0[ah]);
        _prob1 = _mm256_mul_pd(_prob1, _emit1[ah]);
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i+0], _prob0);
        _mm256_store_pd(&prob[i+4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
*/

inline
void haplotype_segment_single::COLLAPSE_AMB() {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: true no-op (preserve probability state, no DP)
            return;
        }

        // Classify and dispatch to supersite logic
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_COLLAPSE_HOM(ss, ss_idx, c0); return;
            case SSClass::AMB: SS_COLLAPSE_AMB(ss, ss_idx, c0, c1); return;
        }
    }
    
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.error_ratio[curr_abs_locus]:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.error_ratio[curr_abs_locus];
    }
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch, true);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);
    __m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }
        _prob = _mm256_mul_ps(_prob, _emit[ah]);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************            MISSING GENOTYPE            ************************/
/*******************************************************************************/

inline
void haplotype_segment_single::INIT_MIS() {
    fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
    fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
    probSumT = 1.0f;
}

inline
void haplotype_segment_single::RUN_MIS() {
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
    __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
    _tFreq = _mm256_mul_ps(_tFreq, _factor);
    __m256 _nt = _mm256_set1_ps(nt / probSumT);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        __m256 _prob = _mm256_load_ps(&prob[i]);
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_single::COLLAPSE_MIS() {
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch, true);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        __m256 _prob;
        if (use_outer) {
            float row_mix = row_stay * probSumK[k] + row_switch;
            _prob = _mm256_mul_ps(col_mix, _mm256_set1_ps(row_mix));
        } else {
            _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        }
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************                    SUM Ks                ************************/
/*******************************************************************************/

inline
void haplotype_segment_single::SUMK() {
        for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
            probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
        }
}

/*******************************************************************************/
/*****************        TRANSITION COMPUTATIONS            ************************/
/*******************************************************************************/

inline
bool haplotype_segment_single::TRANS_HAP() {
    sumHProbs = 0.0f;
    unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
    yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], prev_abs_locus);
    nt = 1.0f - yt;
    
    // Guard against zero/near-zero AlphaSumSum from previous segment.
    // Root cause: If prev segment has Lengths[s] = 0 (empty segment), the forward pass storage
    // code (haplotype_segment_single.cpp:291-295) never executes because the condition
    // `curr_segment_locus == Lengths[s] - 1` becomes `0 == MAX_UINT` (unsigned wraparound),
    // leaving AlphaSumSum[s] at its initialized value of 0.0f.
    // This can occur with genotypes having zero variants or other edge cases in segment building.
    // Signal underflow to trigger double-precision retry in production pipeline.
    const float prev_total = AlphaSumSum[curr_rel_segment_index - 1];
    const float epsilon = 1e-30f;
    if (prev_total < epsilon) {
        // Return true to signal underflow and trigger double-precision fallback
        sumHProbs = 0.0f;
        return true;
    }
    
    // Log inputs to TRANS_HAP for first seg=2 call
    bool is_test_sample = (G->name.find("_sample") != std::string::npos);
    bool is_seg2 = (curr_segment_index == 2);
    bool is_first_iter = (curr_abs_locus == 15 || curr_abs_locus == 30);
    if (is_test_sample && is_seg2 && is_first_iter) {
        std::fprintf(stderr, "[TRANS_HAP_INPUTS] sample=%s locus=%d seg=%d\n",
                     G->name.c_str(), curr_abs_locus, curr_segment_index);
        std::fprintf(stderr, "  prev_total=%.20g yt=%.20g nt=%.20g n_cond=%u\n",
                     (double)prev_total, (double)yt, (double)nt, n_cond_haps);
        std::fprintf(stderr, "  AlphaSum[prev_seg][0-3]: %.20g %.20g %.20g %.20g\n",
                     (double)AlphaSum[curr_rel_segment_index-1][0],
                     (double)AlphaSum[curr_rel_segment_index-1][1],
                     (double)AlphaSum[curr_rel_segment_index-1][2],
                     (double)AlphaSum[curr_rel_segment_index-1][3]);
        std::fprintf(stderr, "  prob[0][0-7]: %.20g %.20g %.20g %.20g %.20g %.20g %.20g %.20g\n",
                     (double)prob[0], (double)prob[1], (double)prob[2], (double)prob[3],
                     (double)prob[4], (double)prob[5], (double)prob[6], (double)prob[7]);
        std::fprintf(stderr, "  Alpha[prev_seg][0*8+0-3]: %.20g %.20g %.20g %.20g\n",
                     (double)Alpha[curr_rel_segment_index-1][0],
                     (double)Alpha[curr_rel_segment_index-1][1],
                     (double)Alpha[curr_rel_segment_index-1][2],
                     (double)Alpha[curr_rel_segment_index-1][3]);
    }

    float fact1 = nt / prev_total;
    for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
        __m256 _sum = _mm256_set1_ps(0.0f);
        float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/prev_total) * yt / n_cond_haps;
        for (int k = 0 ; k < n_cond_haps ; k ++) {
            __m256 _alpha = _mm256_set1_ps(Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1] * fact1 + fact2);
            __m256 _beta = _mm256_load_ps(&prob[k*HAP_NUMBER]);
            _sum = _mm256_add_ps(_sum, _mm256_mul_ps(_alpha, _beta));
        }
        _mm256_store_ps(&HProbs[h1*HAP_NUMBER], _sum);
        sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
    }
    bool enable_debug = (!debug::SUPERDEBUG_SAMPLENAME.empty() && G->name == debug::SUPERDEBUG_SAMPLENAME) ||
                        (G->name.find("_sample") != std::string::npos);
    if (enable_debug && curr_segment_index >= 1 && curr_segment_index <= 5) {
        std::fprintf(stderr, "[TRANS_HAP] sample=%s locus=%d seg=%d prev_total=%.15f yt=%.15f nt=%.15f n_cond=%u sumHProbs=%.15f\n",
                     G->name.c_str(), curr_abs_locus, curr_segment_index, (double)prev_total, (double)yt, (double)nt, n_cond_haps, (double)sumHProbs);
        // Show first 4 prob values
        for (int k = 0; k < std::min(4, (int)n_cond_haps); k++) {
            std::fprintf(stderr, "  prob[%d][0-3]=%.6g %.6g %.6g %.6g\n", k,
                         (double)prob[k*HAP_NUMBER+0], (double)prob[k*HAP_NUMBER+1],
                         (double)prob[k*HAP_NUMBER+2], (double)prob[k*HAP_NUMBER+3]);
        }
    }
    return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<float>::min());
}

inline
bool haplotype_segment_single::TRANS_DIP_MULT() {
    sumDProbs= 0.0f;
    double scaling = 1.0 / sumHProbs;
    int pd_hits = 0, nd_hits = 0, t_hits = 0;

    // Log HProbs array and diplotype iteration for first seg=2 call per sample
    bool is_test_sample = (G->name.find("_sample") != std::string::npos);
    bool is_seg2 = (curr_segment_index == 2);
    bool is_first_iter = (curr_abs_locus == 15 || curr_abs_locus == 30); // bial locus 15, supersite locus 30
    if (is_test_sample && is_seg2 && is_first_iter) {
        std::fprintf(stderr, "[HPROBS_ARRAY] sample=%s locus=%d seg=%d sumHProbs=%.20g\n",
                     G->name.c_str(), curr_abs_locus, curr_segment_index, (double)sumHProbs);
        std::fprintf(stderr, "  HProbs[0-15]:\n");
        for (int i = 0; i < 16; i++) {
            std::fprintf(stderr, "    [%d] %.20g (hex:%a)\n", i, (double)HProbs[i], (double)HProbs[i]);
        }
        std::fprintf(stderr, "  HProbs[16-31]:\n");
        for (int i = 16; i < 32; i++) {
            std::fprintf(stderr, "    [%d] %.20g (hex:%a)\n", i, (double)HProbs[i], (double)HProbs[i]);
        }
        std::fprintf(stderr, "  Diplotypes prev[%d]=0x%016llx curr[%d]=0x%016llx\n",
                     curr_segment_index-1, (unsigned long long)G->Diplotypes[curr_segment_index-1],
                     curr_segment_index, (unsigned long long)G->Diplotypes[curr_segment_index]);
    }

    for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
        if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
            pd_hits++;
            for (int nd = 0 ; nd < 64 ; ++nd) {
                if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
                    if (pd_hits == 1) nd_hits++;
                    DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
                    if (is_test_sample && is_seg2 && is_first_iter && t < 8) {
                        int h0_pd = DIP_HAP0(pd);
                        int h1_pd = DIP_HAP1(pd);
                        int h0_nd = DIP_HAP0(nd);
                        int h1_nd = DIP_HAP1(nd);
                        std::fprintf(stderr, "    t=%d pd=%d(h%d,h%d) nd=%d(h%d,h%d) DProbs[%d]=%.20g\n",
                                     t, pd, h0_pd, h1_pd, nd, h0_nd, h1_nd, t, DProbs[t]);
                    }
                    sumDProbs += DProbs[t];
                    t_hits++;
                    t++;
                }
            }
        }
    }

    if (is_test_sample && is_seg2 && is_first_iter) {
        std::fprintf(stderr, "  Total: pd_hits=%d nd_hits=%d t_hits=%d sumDProbs=%.20g\n",
                     pd_hits, nd_hits, t_hits, sumDProbs);
    }

    if (trans_parity_trace_enabled_s()) {
        std::fprintf(stderr, "[TRANS_DIP_MULT debug][single] seg=%d sumHProbs=%g scaling=%g pd_hits=%d nd_hits_first_pd=%d t_hits=%d sumDProbs=%g firstD=%g\n",
                     curr_segment_index, sumHProbs, scaling, pd_hits, nd_hits, t_hits, sumDProbs, (t_hits>0?DProbs[0]:-1.0));
        std::fflush(stderr);
    }
    return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

inline
bool haplotype_segment_single::TRANS_DIP_ADD() {
    sumDProbs = 0.0f;
    double scaling = 1.0 / sumHProbs;
    for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
        if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
            for (int nd = 0 ; nd < 64 ; ++nd) {
                if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
                    DProbs[t] = DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) + ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
                    sumDProbs += DProbs[t];
                    t++;
                }
            }
        }
    }
    return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

inline
void haplotype_segment_single::IMPUTE(std::vector < float > & missing_probabilities) {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];

    auto deep_trace_enabled = []() {
        static int flag = -1;
        if (flag < 0) {
            const char* env = std::getenv("SHAPEIT5_IMPUTE_DEEP");
            flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
        }
        return flag == 1;
    };

    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[IMPUTE_TRACE] locus=%d rel_mis=%d ss_idx=%d n_cond=%u\n",
                     curr_abs_locus, curr_rel_missing, ss_idx, n_cond_haps);
        std::fprintf(stderr, "  AlphaSum: ");
        for (int h=0; h<HAP_NUMBER; ++h) std::fprintf(stderr, "%.6e ", AlphaSumMissing[curr_rel_missing][h]);
        std::fprintf(stderr, "\n");
    }

    if (deep_trace_enabled()) {
        std::fprintf(stderr,
                     "[IMPUTE_DEEP] locus=%d ss_idx=%d curr_abs_missing=%d curr_rel_missing=%d locus_first=%d\n",
                     curr_abs_locus, ss_idx, curr_abs_missing, curr_rel_missing, locus_first);
        std::fprintf(stderr, "  AlphaSumMissing per lane:");
        for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " %.8e", AlphaSumMissing[curr_rel_missing][h]);
        std::fprintf(stderr, "\n");
        // Dump AlphaMissing/Beta for all donors
        for (unsigned int k = 0, i = 0; k < n_cond_haps; ++k, i += HAP_NUMBER) {
            std::fprintf(stderr, "  [A/B] k=%u", k);
            for (int h = 0; h < HAP_NUMBER; ++h) {
                std::fprintf(stderr, " A%.8e", AlphaMissing[curr_rel_missing][i + h]);
            }
            for (int h = 0; h < HAP_NUMBER; ++h) {
                std::fprintf(stderr, " B%.8e", prob[i + h]); // prob holds Beta in backward
            }
            std::fprintf(stderr, "\n");
        }
    }

    __m256 _sum = _mm256_set1_ps(0.0f);
    // NOTE: We intentionally do NOT divide by AlphaSumMissing here.
    // The 1/AlphaSum factor cancels out in the final ratio normalization.
    // Removing it avoids overflow when AlphaSum is subnormal (< FLT_MIN).

    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Determine class for current split variant (0=REF, 1..n_alts)
        int target_class = 0;
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { target_class = (int)ai + 1; break; }
        }
        ss_load_cond_codes(ss, ss_idx);
        // Prepare accumulators per class
        __m256 sum_classes[SUPERSITE_MAX_ALTS + 1];
        for (int c = 0; c <= (int)ss.var_count; ++c) sum_classes[c] = _mm256_set1_ps(0.0f);
        __m256 denom = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int code = (int)ss_cond_codes[k]; if (code > (int)ss.var_count) code = 0;
            __m256 _prob = _mm256_load_ps(&prob[i]);
            __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
            _sum = _mm256_mul_ps(_alpha, _prob);  // No AlphaSum division needed
            sum_classes[code] = _mm256_add_ps(sum_classes[code], _sum);
            denom = _mm256_add_ps(denom, _sum);
        }
        __m256 p = sum_classes[target_class];
        alignas(32) float pv[8], dv[8];
        alignas(32) float denom_lanes[HAP_NUMBER];
        alignas(32) float class_lanes[SUPERSITE_MAX_ALTS + 1][HAP_NUMBER];
        _mm256_store_ps(pv, p);
        _mm256_store_ps(dv, denom);
        _mm256_store_ps(denom_lanes, denom);
        for (int c = 0; c <= (int)ss.var_count; ++c) _mm256_store_ps(class_lanes[c], sum_classes[c]);
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float d = dv[h];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = (d > 0.0f) ? (pv[h] / d) : 0.0f;
        }
        // Structured TSV logging for supersite splits (non-anchor)
        if (impute_tsv_should_log(G)) {
            int C = static_cast<int>(ss.var_count) + 1;
            for (int h = 0; h < HAP_NUMBER; ++h) {
                float lane_nums[SUPERSITE_MAX_ALTS + 1];
                float lane_sc[SUPERSITE_MAX_ALTS + 1];
                const float denom_val = denom_lanes[h];
                float uniform = 1.0f / static_cast<float>(C);
                for (int c = 0; c < C; ++c) {
                    lane_nums[c] = class_lanes[c][h];
                    lane_sc[c] = (denom_val > 0.0f) ? (class_lanes[c][h] / denom_val) : uniform;
                }
                impute_tsv_log_row("supersite-split",
                                   G,
                                   curr_abs_locus,
                                   ss.global_site_id,
                                   ss_idx,
                                   curr_rel_missing,
                                   C,
                                   n_cond_haps,
                                   h,
                                   AlphaSumMissing[curr_rel_missing][h],
                                   denom_val,
                                   lane_nums,
                                   lane_sc);
            }
        }
    } else {
        __m256 _sumA[2]; _sumA[0] = _mm256_set1_ps(0.0f); _sumA[1] = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
            _sum = _mm256_mul_ps(_alpha, _prob);  // No AlphaSum division needed
            _sumA[ah] = _mm256_add_ps(_sumA[ah], _sum);

            if (supersite_trace_enabled()) {
                alignas(32) float av[8], bv[8], sv[8];
                _mm256_store_ps(av, _alpha);
                _mm256_store_ps(bv, _prob);
                _mm256_store_ps(sv, _sum);
                std::fprintf(stderr, "  k=%d ah=%d Alpha=%.2e Beta=%.2e Sum=%.2e\n", k, (int)ah, av[0], bv[0], sv[0]);
            }
        }
        float* prob0 = (float*)&_sumA[0];
        float* prob1 = (float*)&_sumA[1];
        for (int h = 0; h < HAP_NUMBER; h++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h] + prob1[h]);
        
        if (supersite_trace_enabled()) {
            std::fprintf(stderr, "  Result Probs (Alt):");
            for (int h=0; h<HAP_NUMBER; ++h) std::fprintf(stderr, " %.6f", missing_probabilities[curr_abs_missing * HAP_NUMBER + h]);
            std::fprintf(stderr, "\n");
        }

        // Structured TSV logging for biallelic missing imputations
        if (impute_tsv_should_log(G)) {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                float denom = prob0[h] + prob1[h];
                float lane_nums[2] = {prob0[h], prob1[h]};
                float lane_sc[2];
                if (denom > 0.0f) {
                    lane_sc[0] = prob0[h] / denom;
                    lane_sc[1] = prob1[h] / denom;
                } else {
                    lane_sc[0] = 0.5f;
                    lane_sc[1] = 0.5f;
                }
                impute_tsv_log_row("bial",
                                   G,
                                   curr_abs_locus,
                                   curr_abs_locus,
                                   -1,
                                   curr_rel_missing,
                                   2,
                                   n_cond_haps,
                                   h,
                                   AlphaSumMissing[curr_rel_missing][h],
                                   denom,
                                   lane_nums,
                                   lane_sc);
            }
        }

        // Optional side-by-side SC vs bial probability dump for comparison debugging
        static int sc_compare = -1;
        if (sc_compare < 0) {
            const char* env = std::getenv("SHAPEIT5_SC_COMPARE");
            sc_compare = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
        }
        if (sc_compare == 1) {
            std::fprintf(stderr, "[SC_COMPARE_BIAL] sample=%s locus=%d rel_mis=%d ss_idx=%d\n",
                         G ? G->name.c_str() : "<nil>", curr_abs_locus, curr_rel_missing, ss_idx);
            for (int h = 0; h < HAP_NUMBER; ++h) {
                float p_alt = missing_probabilities[curr_abs_missing * HAP_NUMBER + h];
                float p_ref = 1.0f - p_alt;
                std::fprintf(stderr, "  lane=%d p_ref=%.6f p_alt=%.6f\n", h, p_ref, p_alt);
            }
        }
    }
}

// Phase 3: Impute multivariant posteriors for missing supersite
// Computes P(class_c | Alpha, Beta) for each class c ∈ {REF, ALT1, ..., ALTn}
// Storage: SC[offset + h*C + c] where h=lane, C=num_classes, c=class_index
inline
void haplotype_segment_single::IMPUTE_SUPERSITE_MULTIVARIATE(
    std::vector<float>& SC,
    const SuperSite& ss,
    int ss_idx,
    int rel_missing_index,
    const std::vector<uint32_t>* supersite_sc_offset)
{
    // Return early if panel codes not available (no supersites built)
    if (panel_codes == nullptr) return;

    ss_load_cond_codes(ss, ss_idx);

    int C = (int)ss.n_classes;  // 1 + n_alts

    // Calculate offset - use thread-local offset if available, otherwise fallback to simple calculation
    uint32_t offset;
    if (supersite_sc_offset != nullptr) {
        offset = (*supersite_sc_offset)[ss_idx];  // Use thread-local offset to prevent race conditions
    } else {
        // Fallback for tests: assume simple layout with supersite index * (lanes * classes)
        const bool has_guard = (SC.size() >= 2 && std::isnan(SC[0]) && std::isnan(SC[SC.size() - 1]));
        offset = static_cast<uint32_t>(ss_idx) * HAP_NUMBER * C + (has_guard ? 1u : 0u);
    }

    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "[IMPUTE_SS_TRACE] locus=%d rel_mis=%d ss_idx=%d n_cond=%u offset=%u C=%d\n",
                     curr_abs_locus, rel_missing_index, ss_idx, n_cond_haps, offset, C);
        std::fprintf(stderr, "  AlphaSum: ");
        for (int h=0; h<HAP_NUMBER; ++h) std::fprintf(stderr, "%.6e ", AlphaSumMissing[rel_missing_index][h]);
        std::fprintf(stderr, "\n");
        // Show donor mapping for the first few conditioning states to confirm parity
        unsigned int dbg_k = std::min(16u, n_cond_haps);
        for (unsigned int k = 0; k < dbg_k; ++k) {
            int cond_row = (cond_idx && k < cond_idx->size()) ? (int)(*cond_idx)[k] : -1;
            bool bial_alt = Hvar.get(curr_abs_locus, k);
            std::fprintf(stderr, "  [SS_DONOR_MAP] k=%u cond_idx=%d code=%u bial_alt=%d\n",
                         k, cond_row, (unsigned)ss_cond_codes[k], (int)bial_alt);
        }
        // Dump the biallelic panel allele for each donor to check cond_idx consistency
        fprintf(stderr, "  [SS_DONOR_ALLELES] locus=%d n_cond=%u\n", curr_abs_locus, (unsigned)n_cond_haps);
        for (unsigned int k = 0; k < dbg_k; ++k) {
            bool allele = Hvar.get(curr_abs_locus, k);
            std::fprintf(stderr, "    k=%u allele=%d\n", k, (int)allele);
        }
    }

    const size_t payload = static_cast<size_t>(HAP_NUMBER) * static_cast<size_t>(C);
    const bool guard_sentinels_present = supersite_debug::guard_checks_enabled() &&
                                         SC.size() >= 2 &&
                                         std::isnan(SC.front()) &&
                                         std::isnan(SC.back());
    const size_t base_offset = guard_sentinels_present ? 1u : 0u;  // first data slot
    const size_t usable = guard_sentinels_present ? (SC.size() - 2) : SC.size();

    if (supersite_debug::guard_checks_enabled()) {
        // Validate sentinel integrity and bounds before writing
        if (guard_sentinels_present) {
            if (!std::isnan(SC.front()) || !std::isnan(SC.back())) {
                std::fprintf(stderr,
                    "[supersite-guard] Sentinel corrupt before write sample=%s ss_idx=%d front=%g back=%g\n",
                    G ? G->name.c_str() : "?", ss_idx, SC.front(), SC.back());
                assert(std::isnan(SC.front()) && std::isnan(SC.back()));
                return;
            }
        }
        const size_t rel_offset = (offset >= base_offset) ? static_cast<size_t>(offset - base_offset) : static_cast<size_t>(-1);
        if (offset < base_offset || payload > usable || rel_offset > usable || rel_offset + payload > usable) {
            std::fprintf(stderr,
                "[supersite-guard] IMPUTE_SUPERSITE_MULTIVARIATE bounds: offset=%u base=%zu rel=%zu C=%d size=%zu usable=%zu payload=%zu sample=%s ss_idx=%d\n",
                offset, base_offset, rel_offset, C, SC.size(), usable, payload, G ? G->name.c_str() : "?", ss_idx);
            assert(offset >= base_offset);
            assert(payload <= usable);
            assert(rel_offset <= usable);
            assert(rel_offset + payload <= usable);
            return;
        }
    }
    
    // Initialize per-class accumulators (8 lanes = 8 samples)
    // sum[c] = Σ_k [Alpha_k × Beta_k × 1{donor_k carries class_c}]
    __m256 sum[SUPERSITE_MAX_ALTS + 1];
    for (int c = 0; c < C; ++c) {
        sum[c] = _mm256_set1_ps(0.0f);
    }
    __m256 denom = _mm256_set1_ps(0.0f);
    
    if (rel_missing_index < 0 || rel_missing_index >= (int)AlphaSumMissing.size()) {
        std::fprintf(stderr,
            "[supersite-guard] AlphaSumMissing OOB sample=%s ss_idx=%d rel_missing_index=%d size=%zu\n",
            G ? G->name.c_str() : "?", ss_idx, rel_missing_index, AlphaSumMissing.size());
        assert(rel_missing_index >= 0 && rel_missing_index < (int)AlphaSumMissing.size());
        return;
    }
    if (rel_missing_index >= (int)AlphaMissing.size()) {
        std::fprintf(stderr,
            "[supersite-guard] AlphaMissing OOB sample=%s ss_idx=%d rel_missing_index=%d size=%zu\n",
            G ? G->name.c_str() : "?", ss_idx, rel_missing_index, AlphaMissing.size());
        assert(rel_missing_index < (int)AlphaMissing.size());
        return;
    }
    if (AlphaSumMissing[rel_missing_index].size() < HAP_NUMBER ||
        AlphaMissing[rel_missing_index].size() < HAP_NUMBER) {
        std::fprintf(stderr,
            "[supersite-guard] Missing arrays under-sized sample=%s ss_idx=%d rel_missing_index=%d AlphaSumMissing_size=%zu AlphaMissing_size=%zu\n",
            G ? G->name.c_str() : "?", ss_idx, rel_missing_index,
            AlphaSumMissing[rel_missing_index].size(), AlphaMissing[rel_missing_index].size());
        assert(AlphaSumMissing[rel_missing_index].size() >= HAP_NUMBER);
        assert(AlphaMissing[rel_missing_index].size() >= HAP_NUMBER);
        return;
    }

    if (supersite_trace_enabled()) {
        supersite_trace_log("[SupersiteImpute] enter ss_idx=%d offset=%u C=%d rel_missing_index=%d n_cond_haps=%u\n",
                            ss_idx, offset, C, rel_missing_index, (unsigned)n_cond_haps);
    }

    // NOTE: We intentionally do NOT divide by AlphaSumMissing here.
    // The 1/AlphaSum factor would be a per-lane constant that multiplies every term,
    // so it cancels out in the final SC = sum[c] / denom normalization.
    // Removing it avoids overflow when AlphaSum is subnormal (< FLT_MIN).

    // Accumulate: for each conditioning donor haplotype k
    for (int k = 0, i = 0; k < (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        uint8_t code = ss_cond_codes[k];  // 0=REF, 1..n_alts
        if (code >= C) code = 0;  // Safety: invalid code -> REF

        // term = Alpha[k] x Beta[k] (8 lanes) - no AlphaSum division needed
        __m256 _alpha = _mm256_load_ps(&AlphaMissing[rel_missing_index][i]);
        __m256 _beta  = _mm256_load_ps(&prob[i]);  // prob=Beta in backward pass
        __m256 term = _mm256_mul_ps(_alpha, _beta);

        // Add to class bucket
        sum[code] = _mm256_add_ps(sum[code], term);
        denom = _mm256_add_ps(denom, term);

        if (supersite_trace_enabled()) {
            alignas(32) float av[8], bv[8], tv[8];
            _mm256_store_ps(av, _alpha);
            _mm256_store_ps(bv, _beta);
            _mm256_store_ps(tv, term);
            std::fprintf(stderr, "  k=%d code=%u Alpha=%.2e Beta=%.2e Term=%.2e\n", k, code, av[0], bv[0], tv[0]);
        }
    }
    if (supersite_trace_enabled()) supersite_trace_log("[SupersiteImpute] accumulation complete\n");

    // Normalize: SC[offset + h*C + c] = sum[c][h] / denom[h]
    // Extract 8 lanes and write to SC buffer
    float sum_lanes[SUPERSITE_MAX_ALTS + 1][HAP_NUMBER];
    float denom_lanes[HAP_NUMBER];

        for (int c = 0; c < C; ++c) {
        _mm256_storeu_ps(sum_lanes[c], sum[c]);
    }
    _mm256_storeu_ps(denom_lanes, denom);

    // Write normalized posteriors to SC
    for (int h = 0; h < HAP_NUMBER; ++h) {
        float d = denom_lanes[h];
        if (d > 0.0f) {
            for (int c = 0; c < C; ++c) {
                SC[offset + h * C + c] = sum_lanes[c][h] / d;
            }
        } else {
            // Fallback: uniform distribution if denominator is zero
            float uniform = 1.0f / (float)C;
            for (int c = 0; c < C; ++c) {
                SC[offset + h * C + c] = uniform;
            }
        }
    }

    // Structured TSV logging for supersite anchor missing imputations
    if (impute_tsv_should_log(G)) {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float lane_nums[SUPERSITE_MAX_ALTS + 1];
            float lane_sc[SUPERSITE_MAX_ALTS + 1];
            float denom_val = denom_lanes[h];
            float uniform = 1.0f / (float)C;
            for (int c = 0; c < C; ++c) {
                lane_nums[c] = sum_lanes[c][h];
                lane_sc[c] = (denom_val > 0.0f) ? SC[offset + h * C + c] : uniform;
            }
            impute_tsv_log_row("supersite-anchor",
                               G,
                               curr_abs_locus,
                               ss.global_site_id,
                               ss_idx,
                               rel_missing_index,
                               C,
                               n_cond_haps,
                               h,
                               AlphaSumMissing[rel_missing_index][h],
                               denom_val,
                               lane_nums,
                               lane_sc);
        }
    }

    if (supersite_trace_enabled()) {
        std::fprintf(stderr, "  Result SC (h0):");
        for (int c=0; c<C; ++c) std::fprintf(stderr, " %.6f", SC[offset + 0*C + c]);
        std::fprintf(stderr, "\n");
        }

        // Guard: row-normalization and finiteness check for SC
        if (supersite_debug::guard_checks_enabled()) {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                float sum_row = 0.0f;
                bool finite = true;
                for (int c = 0; c < C; ++c) {
                    const float val = SC[offset + h * C + c];
                    sum_row += val;
                    if (!std::isfinite(val)) {
                        finite = false;
                    }
                }
                if (!finite || std::fabs(sum_row - 1.0f) > 1e-3f) {
                    std::fprintf(stderr,
                        "[supersite-guard] SC row anomaly sample=%s ss_idx=%d hap=%d sum=%.6f finite=%d\n",
                        G ? G->name.c_str() : "?", ss_idx, h, sum_row, finite ? 1 : 0);
                }
            }
        }

        // Deep trace: dump all inputs/outputs for manual parity math
        auto deep_trace_enabled = []() {
            static int flag = -1;
            if (flag < 0) {
                const char* env = std::getenv("SHAPEIT5_IMPUTE_DEEP");
                flag = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
            }
            return flag == 1;
        };
        if (deep_trace_enabled()) {
            std::fprintf(stderr,
                         "[IMPUTE_DEEP_SS] locus=%d ss_idx=%d rel_mis=%d C=%d curr_abs_missing=%d locus_first=%d\n",
                         curr_abs_locus, ss_idx, rel_missing_index, C, curr_abs_missing, locus_first);
            std::fprintf(stderr, "  AlphaSumMissing per lane:");
            for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " %.8e", AlphaSumMissing[rel_missing_index][h]);
            std::fprintf(stderr, "\n");
            for (unsigned int k = 0, i = 0; k < n_cond_haps; ++k, i += HAP_NUMBER) {
                std::fprintf(stderr, "  [SS_A/B] k=%u code=%u", k, (unsigned)ss_cond_codes[k]);
                for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " A%.8e", AlphaMissing[rel_missing_index][i + h]);
                for (int h = 0; h < HAP_NUMBER; ++h) std::fprintf(stderr, " B%.8e", prob[i + h]);
                std::fprintf(stderr, "\n");
            }
            for (int h = 0; h < HAP_NUMBER; ++h) {
                std::fprintf(stderr, "  lane=%d denom=%.8e", h, denom_lanes[h]);
                for (int c = 0; c < C; ++c) {
                    std::fprintf(stderr, " num_c%d=%.8e", c, sum_lanes[c][h]);
                }
                std::fprintf(stderr, "\n");
            }
        }

        // Optional side-by-side SC dump for comparison debugging
        static int sc_compare = -1;
        if (sc_compare < 0) {
            const char* env = std::getenv("SHAPEIT5_SC_COMPARE");
            sc_compare = (env && env[0] != '\0' && env[0] != '0') ? 1 : 0;
        }
        if (sc_compare == 1) {
            std::fprintf(stderr, "[SC_COMPARE_SS] sample=%s locus=%d ss_idx=%d rel_mis=%d C=%d\n",
                         G ? G->name.c_str() : "<nil>", curr_abs_locus, ss_idx, rel_missing_index, C);
            for (int h = 0; h < HAP_NUMBER; ++h) {
                std::fprintf(stderr, "  lane=%d", h);
                for (int c = 0; c < C; ++c) {
                    std::fprintf(stderr, " c%d=%.6f", c, SC[offset + h * C + c]);
                }
                std::fprintf(stderr, "\n");
            }
        }


    if (supersite_trace_enabled()) supersite_trace_log("[SupersiteImpute] wrote SC block\n");

    if (supersite_debug::guard_checks_enabled()) {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float row_sum = 0.0f;
            for (int c = 0; c < C; ++c) {
                float value = SC[offset + h * C + c];
                if (!std::isfinite(value)) {
                    std::fprintf(stderr,
                        "[supersite-guard] Non-finite SC entry sample=%s ss_idx=%d lane=%d class=%d value=%g\n",
                        G ? G->name.c_str() : "?", ss_idx, h, c, value);
                    assert(std::isfinite(value));
                    return;
                }
                if (value < -1e-6f) {
                    std::fprintf(stderr,
                        "[supersite-guard] Negative SC entry sample=%s ss_idx=%d lane=%d class=%d value=%g\n",
                        G ? G->name.c_str() : "?", ss_idx, h, c, value);
                    assert(value >= -1e-6f);
                    return;
                }
                row_sum += value;
            }
            if (std::fabs(row_sum - 1.0f) > 1e-3f) {
                std::fprintf(stderr,
                    "[supersite-guard] Row sum drift sample=%s ss_idx=%d lane=%d sum=%g\n",
                    G ? G->name.c_str() : "?", ss_idx, h, row_sum);
                assert(std::fabs(row_sum - 1.0f) <= 1e-3f);
                return;
            }
        }
        const size_t guard_hi = SC.size() ? SC.size() - 1 : 0;
        if (guard_hi > 0 && guard_sentinels_present) {
            const float pre_guard = SC.front();
            const float post_guard = SC[guard_hi];
            if (!std::isnan(pre_guard) || !std::isnan(post_guard)) {
                std::fprintf(stderr,
                    "[supersite-guard] Guard corrupt after write sample=%s ss_idx=%d pre=%g post=%g\n",
                    G ? G->name.c_str() : "?", ss_idx, pre_guard, post_guard);
                assert(std::isnan(pre_guard));
                assert(std::isnan(post_guard));
            }
        }
    }

    if (supersite_trace_enabled()) {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float row_sum = 0.0f;
            supersite_trace_log("[SupersiteImpute] lane=%d", h);
            for (int c = 0; c < C; ++c) {
                float value = SC[offset + h * C + c];
                row_sum += value;
                supersite_trace_log(" %d:%.4f", c, value);
            }
            supersite_trace_log(" | sum=%.4f\n", row_sum);
        }
    }
}

#endif

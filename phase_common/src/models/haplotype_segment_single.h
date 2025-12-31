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
#include <cmath>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <cstring>

#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>


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
    const std::vector<aligned_vector32<uint8_t>>* get_ss_panel_matrix() const {
        return ss_panel_matrix.empty() ? nullptr : &ss_panel_matrix;
    }

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

    // Simply copy from static matrix using relative indices (like biallelic Hvar)
    ss_cond_codes.resize(n_cond_haps);
    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        ss_cond_codes[k] = ss_panel_matrix[ss_idx][k];
    }

    // The static matrix is our "cache" and it's always valid.
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

        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;
            case SSClass::AMB: INIT_AMB(); return;
            case SSClass::HOM:
                ss_load_cond_codes(ss, ss_idx);
                {
                    const float error_ratio = M.ed / M.ee;
                    __m256 _sum = _mm256_set1_ps(0.0f);
                    for (int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
                        float emission = (c0 == ss_cond_codes[k]) ? 1.0f : error_ratio;
                        __m256 _prob = _mm256_set1_ps(emission);
                        _sum = _mm256_add_ps(_sum, _prob);
                        _mm256_store_ps(&prob[i], _prob);
                    }
                    _mm256_store_ps(&probSumH[0], _sum);
                    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                }

                // Trace: matches per lane at anchor (HOM)
                return;
        }
    }

    // Biallelic path (original simple loop)
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _mm256_set1_ps((ag==ah)?1.0f:M.ed / M.ee);
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

        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return true;
            case SSClass::AMB: RUN_AMB(); return true;
            case SSClass::HOM: {
                const uint8_t sample_code = c0;
                if (ss.rare_code_mask != 0 && sample_code <= ss.n_alts) {
                    if ((ss.rare_code_mask & static_cast<uint16_t>(1u << sample_code)) == 0) return false;
                }

                ss_load_cond_codes(ss, ss_idx);

                // DETAILED BACKWARD TRACING


                // Update probabilities with transitions
                __m256 _sum = _mm256_set1_ps(0.0f);
                __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
                __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
                _tFreq = _mm256_mul_ps(_tFreq, _factor);
                __m256 _nt = _mm256_set1_ps(nt / probSumT);
                const float error_ratio = M.ed / M.ee;
                __m256 _mismatch = _mm256_set1_ps(error_ratio);


                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    __m256 _prob = _mm256_load_ps(&prob[i]);
                    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
                    const bool match = (ss_cond_codes[k] == sample_code);
                    if (!match) _prob = _mm256_mul_ps(_prob, _mismatch);
                    _sum = _mm256_add_ps(_sum, _prob);
                    _mm256_store_ps(&prob[i], _prob);

                }
                _mm256_store_ps(&probSumH[0], _sum);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

                // DETAILED BACKWARD TRACING - OUTPUT

                return true;
            }
        }
    }
    
    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    if (rare_allele < 0 || ag == rare_allele) {
        // EMISSION TRACING

        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
        __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
        _tFreq = _mm256_mul_ps(_tFreq, _factor);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        __m256 _mismatch = _mm256_set1_ps(M.ed / M.ee);

        for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            float emission = (ag != ah) ? M.ed / M.ee : 1.0f;
            if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);

        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];

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
        
        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;
            case SSClass::AMB: COLLAPSE_AMB(); return;
            case SSClass::HOM: {
                const uint8_t sample_code = c0;
                // Load cached conditioning haplotype codes
                ss_load_cond_codes(ss, ss_idx);

                // Precompute emissions
                precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, M.ed / M.ee, ss_emissions);

                // Collapse from probSumK with optional outer-product seeding
                __m256 _sum = _mm256_set1_ps(0.0f);
                __m256 col_mix;
                float row_stay = 0.0f, row_switch = 0.0f;
                int rel_prev_seg = (curr_segment_index - segment_first) - 1;
                const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch);
                __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
                __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);

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

                return;
            }
        }
    }

    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 col_mix;
    float row_stay = 0.0f, row_switch = 0.0f;
    int rel_prev_seg = (curr_segment_index - segment_first) - 1;
    const bool use_outer = prepare_outer_product_mix(rel_prev_seg, col_mix, row_stay, row_switch, true);
    __m256 _tFreq = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(yt / n_cond_haps);
    __m256 _nt = use_outer ? _mm256_set1_ps(0.0f) : _mm256_set1_ps(nt / probSumT);

    __m256 _mismatch = _mm256_set1_ps(M.ed / M.ee);
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

        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: INIT_HOM(); return;
            case SSClass::AMB: {
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
                __m256 mis_f   = _mm256_set1_ps(M.ed / M.ee);
                alignas(32) int expv[HAP_NUMBER];
                for (int h = 0; h < HAP_NUMBER; ++h) expv[h] = (int)expected_class[h];
                __m256i exp_vec = _mm256_load_si256((__m256i*)expv);

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
                return;
            }
        }
    }
    
    // Biallelic path
    unsigned char amb_code = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                             ? G->Ambiguous[curr_abs_ambiguous] : 0u;

    // Biallelic INIT_AMB: simple AVX2 loop (original algorithm)
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed / M.ee;
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
    const char rare_allele = M.rare_allele[curr_abs_locus];
    
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

        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: RUN_HOM(rare_allele); return;
            case SSClass::AMB: {
                // Load cached conditioning haplotype codes
                ss_load_cond_codes(ss, ss_idx);

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



                alignas(32) int expv[HAP_NUMBER];
                for (int h = 0; h < HAP_NUMBER; ++h) expv[h] = (int)expected_class[h];
                __m256i exp_vec = _mm256_load_si256((__m256i*)expv);
                __m256 match_f = _mm256_set1_ps(1.0f);
                __m256 mis_f   = _mm256_set1_ps(M.ed / M.ee);

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

                    _prob = _mm256_mul_ps(_prob, emit);
                    _sum = _mm256_add_ps(_sum, _prob);
                    _mm256_store_ps(&prob[i], _prob);
                }
                _mm256_store_ps(&probSumH[0], _sum);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed / M.ee;
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
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed / M.ee;
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

        // Classify using IMMUTABLE observed genotype (c0, c1)
        uint8_t c0, c1;
        G->getSupersiteObservedGt(ss_idx, c0, c1);
        SSClass cls = classifyObservedGt(c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: COLLAPSE_HOM(); return;
            case SSClass::AMB: {
                // Load cached conditioning haplotype codes
                ss_load_cond_codes(ss, ss_idx);

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
                // Precompute emission vectors based on mode (BUG FIX #5: unified code path)
                __m256 match_f = _mm256_set1_ps(1.0f);
                __m256 mis_f   = _mm256_set1_ps(M.ed / M.ee);
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
                    _prob = _mm256_mul_ps(_prob, emit);
                    _sum = _mm256_add_ps(_sum, _prob);
                    _mm256_store_ps(&prob[i], _prob);
                }
                
                _mm256_store_ps(&probSumH[0], _sum);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed / M.ee;
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
    return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<float>::min());
}

inline
bool haplotype_segment_single::TRANS_DIP_MULT() {
    sumDProbs= 0.0f;
    double scaling = 1.0 / sumHProbs;
    for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
        if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
            for (int nd = 0 ; nd < 64 ; ++nd) {
                if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
                    DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
                    sumDProbs += DProbs[t];
                    t++;
                }
            }
        }
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
        _mm256_store_ps(pv, p);
        _mm256_store_ps(dv, denom);
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float d = dv[h];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = (d > 0.0f) ? (pv[h] / d) : 0.0f;
        }
    } else {
        __m256 _sumA[2]; _sumA[0] = _mm256_set1_ps(0.0f); _sumA[1] = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
            _sum = _mm256_mul_ps(_alpha, _prob);  // No AlphaSum division needed
            _sumA[ah] = _mm256_add_ps(_sumA[ah], _sum);

        }
        float* prob0 = (float*)&_sumA[0];
        float* prob1 = (float*)&_sumA[1];
        for (int h = 0; h < HAP_NUMBER; h++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h] + prob1[h]);
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

    uint32_t offset = 0;
    if (supersite_sc_offset != nullptr) {
        offset = (*supersite_sc_offset)[ss_idx];  // Use thread-local offset to prevent race conditions
    } else {
        // Fallback for tests: assume simple layout with supersite index * (lanes * classes)
        offset = static_cast<uint32_t>(ss_idx) * HAP_NUMBER * C;
    }
    
    // Initialize per-class accumulators (8 lanes = 8 samples)
    // sum[c] = Σ_k [Alpha_k × Beta_k × 1{donor_k carries class_c}]
    __m256 sum[SUPERSITE_MAX_ALTS + 1];
    for (int c = 0; c < C; ++c) {
        sum[c] = _mm256_set1_ps(0.0f);
    }
    __m256 denom = _mm256_set1_ps(0.0f);
    
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

    }

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

}

#endif

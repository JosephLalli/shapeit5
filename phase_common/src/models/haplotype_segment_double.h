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

#ifndef _HAPLOTYPE_SEGMENT_DOUBLE_H
#define _HAPLOTYPE_SEGMENT_DOUBLE_H

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>
#include <models/site_emission_types.h>
#include <models/super_site_emissions.h>

#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>

class haplotype_segment_double {
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
	double probSumT;
	aligned_vector32 < double > prob;
    aligned_vector32 < double > probSumK;
    aligned_vector32 < double > probSumH;
	std::vector < aligned_vector32 < double > > Alpha;
	std::vector < aligned_vector32 < double > > AlphaSum;
	std::vector < int > AlphaLocus;
	aligned_vector32 < double > AlphaSumSum;
	std::vector < aligned_vector32 < double > > AlphaMissing;
	std::vector < aligned_vector32 < double > > AlphaSumMissing;
	double HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));

	//STATIC ARRAYS
	double sumHProbs;
	double sumDProbs;
	double g0[HAP_NUMBER], g1[HAP_NUMBER];
	double nt, yt;

	//SUPER-SITE SUPPORT
	const std::vector<SuperSite>* super_sites;
	const std::vector<int>* locus_to_super_idx;
	const uint8_t* panel_codes;
	const std::vector<int>* super_site_var_index;
	const std::vector<unsigned int>* cond_idx;
	const std::vector<uint32_t>* supersite_sc_offset;  // Thread-local SC offsets (set during backward)
	std::vector<aligned_vector32<uint8_t>> ss_panel_matrix;
	aligned_vector32<uint8_t> ss_cond_codes;
	aligned_vector32<double> ss_emissions;
	aligned_vector32<double> ss_emissions_h1;
	const bool supersites_enabled_flag;

	// Anchor MIS mapping: record rel-missing index per locus in window
	std::vector<int> missing_index_by_locus;

	//INLINED AND UNROLLED ROUTINES
	void INIT_HOM();
	void INIT_AMB();
	void INIT_MIS();
	void handle_sibling_bookkeeping(const SiteView& site_view);
	void INIT_SIB(const SiteView& site_view);
	void RUN_SIB(const SiteView& site_view);
	void COLLAPSE_SIB(const SiteView& site_view);
	bool RUN_HOM(char);
	void RUN_AMB();
	void RUN_MIS();
	void COLLAPSE_HOM();
	void COLLAPSE_AMB();
	void COLLAPSE_MIS();
	
	// Caching helper
	void ss_load_cond_codes(const SuperSite& ss, int ss_idx);
	
	void SUMK();
	void IMPUTE(std::vector < float > & );
	void IMPUTE_SUPERSITE_MULTIVARIATE(std::vector < float > & SC, const SuperSite& ss, int ss_idx, int rel_missing_index, const std::vector<uint32_t>* supersite_sc_offset);  // Phase 3
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(std::vector < double > & );
	int SET_OTHER_TRANS(std::vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment_double(genotype *, bitmatrix &, std::vector < unsigned int > &, window &, hmm_parameters &);
	~haplotype_segment_double();

	//void fetch();
	void forward();
	int backward(std::vector < double > &, std::vector < float > &, 
	            std::vector < float > * SC = nullptr,
	            const std::vector < bool > * anchor_has_missing = nullptr,
	            const std::vector<uint32_t>* supersite_sc_offset = nullptr);  // Phase 3: optional supersite posteriors
	
};

/*******************************************************************************/
/*****************     SUPERSITE HELPER FUNCTIONS (Phase 2)     ***************/
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

// Phase 3: Caching helper to load conditioning haplotype codes once per supersite
inline
void haplotype_segment_double::ss_load_cond_codes(const SuperSite& ss, int ss_idx) {
    if (ss_idx < 0 || ss_idx >= (int)ss_panel_matrix.size()) {
        std::fprintf(stderr, "ERROR: ss_load_cond_codes called with invalid ss_idx=%d (matrix size=%zu)\n",
                     ss_idx, ss_panel_matrix.size());
        std::abort();
    }

    ss_cond_codes.resize(n_cond_haps);
    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        ss_cond_codes[k] = ss_panel_matrix[ss_idx][k];
    }
}

// Note: pack_expected_codes_pd() helper removed - was causing half-lane split antipattern
// Now using proper per-lane emission logic throughout SS_*_AMB() functions

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/


inline
void haplotype_segment_double::INIT_HOM() {
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
            case SSClass::AMB: INIT_AMB(); return;
            case SSClass::HOM: {
                const uint8_t sample_code = c0;
                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                precomputeSuperSiteEmissions_AVX2(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0, M.ed / M.ee, ss_emissions);
                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);
                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    __m256d _emit = _mm256_set1_pd(ss_emissions[k]);
                    __m256d _prob0 = _emit;
                    __m256d _prob1 = _emit;
                    _sum0 = _mm256_add_pd(_sum0, _prob0);
                    _sum1 = _mm256_add_pd(_sum1, _prob1);
                    _mm256_store_pd(&prob[i+0], _prob0);
                    _mm256_store_pd(&prob[i+4], _prob1);
                }
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    // Default biallelic path - simple AVX2 loop (original algorithm)
    const bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _mismatch = _mm256_set1_pd(M.ed / M.ee);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0, _prob1;
        if (ag == ah) {
            _prob0 = _mm256_set1_pd(1.0);
            _prob1 = _mm256_set1_pd(1.0);
        } else {
            _prob0 = _mismatch;
            _prob1 = _mismatch;
        }
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i], _prob0);
        _mm256_store_pd(&prob[i+4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
bool haplotype_segment_double::RUN_HOM(char rare_allele) {
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
            case SSClass::MIS: RUN_MIS(); return true;  // BUG FIX #1: Use biallelic MIS
            case SSClass::AMB: RUN_AMB(); return true;
            case SSClass::HOM: {
                const uint8_t sample_code = c0;
                if (ss.rare_code_mask != 0 && sample_code <= ss.n_alts) {
                    if ((ss.rare_code_mask & static_cast<uint16_t>(1u << sample_code)) == 0) return false;
                }

                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                precomputeSuperSiteEmissions_AVX2(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0, M.ed / M.ee, ss_emissions);
                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);
                __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
                __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
                __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
                _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
                _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
                __m256d _nt = _mm256_set1_pd(nt / probSumT);
                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    __m256d _emit = _mm256_set1_pd(ss_emissions[k]);
                    __m256d _prob0 = _mm256_load_pd(&prob[i]);
                    __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
                    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
                    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
                    if (ss_emissions[k] != 1.0) {
                        _prob0 = _mm256_mul_pd(_prob0, _emit);
                        _prob1 = _mm256_mul_pd(_prob1, _emit);
                    }
                    _sum0 = _mm256_add_pd(_sum0, _prob0);
                    _sum1 = _mm256_add_pd(_sum1, _prob1);
                    _mm256_store_pd(&prob[i], _prob0);
                    _mm256_store_pd(&prob[i+4], _prob1);
                }
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return true;
            }
        }
    }
    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    if (rare_allele < 0 || ag == rare_allele) {
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
        __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
        __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
        __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
		_tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
		_tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
		__m256d _nt = _mm256_set1_pd(nt / probSumT);
		__m256d _mismatch = _mm256_set1_pd(M.ed / M.ee);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
			__m256d _prob0 = _mm256_load_pd(&prob[i]);
			__m256d _prob1 = _mm256_load_pd(&prob[i+4]);
			_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
			_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
			if (ag!=ah) {
				_prob0 = _mm256_mul_pd(_prob0, _mismatch);
				_prob1 = _mm256_mul_pd(_prob1, _mismatch);
			}
			_sum0 = _mm256_add_pd(_sum0, _prob0);
			_sum1 = _mm256_add_pd(_sum1, _prob1);
			_mm256_store_pd(&prob[i], _prob0);
			_mm256_store_pd(&prob[i+4], _prob1);
		}
		_mm256_store_pd(&probSumH[0], _sum0);
		_mm256_store_pd(&probSumH[4], _sum1);
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
		return true;
	}
	return false;
}

inline
void haplotype_segment_double::COLLAPSE_HOM() {
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
            case SSClass::AMB: COLLAPSE_AMB(); return;
            case SSClass::HOM: {
                const uint8_t sample_code = c0;
                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                precomputeSuperSiteEmissions_AVX2(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0, M.ed / M.ee, ss_emissions);
                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);
                __m256d _tFreq0 = _mm256_set1_pd(yt / n_cond_haps);
                __m256d _tFreq1 = _mm256_set1_pd(yt / n_cond_haps);
                __m256d _nt = _mm256_set1_pd(nt / probSumT);
                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    __m256d _emit = _mm256_set1_pd(ss_emissions[k]);
                    __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
                    __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
                    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
                    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
                    if (ss_emissions[k] != 1.0) {
                        _prob0 = _mm256_mul_pd(_prob0, _emit);
                        _prob1 = _mm256_mul_pd(_prob1, _emit);
                    }
                    _sum0 = _mm256_add_pd(_sum0, _prob0);
                    _sum1 = _mm256_add_pd(_sum1, _prob1);
                    _mm256_store_pd(&prob[i], _prob0);
                    _mm256_store_pd(&prob[i+4], _prob1);
                }
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _tFreq0 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _tFreq1 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    __m256d _mismatch = _mm256_set1_pd(M.ed / M.ee);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
        if (ag!=ah) {
            _prob0 = _mm256_mul_pd(_prob0, _mismatch);
            _prob1 = _mm256_mul_pd(_prob1, _mismatch);
        }
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i], _prob0);
        _mm256_store_pd(&prob[i+4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE			********************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_AMB() {
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
                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                unsigned char amb_mask = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                                         ? G->Ambiguous[curr_abs_ambiguous] : 0u;
                uint8_t expected_class[HAP_NUMBER];
                if (c0 == c1) {
                    for (int h = 0; h < HAP_NUMBER; ++h) expected_class[h] = c0;
                } else {
                    for (int h = 0; h < HAP_NUMBER; ++h) {
                        bool use_c1 = ((amb_mask >> h) & 1U);
                        expected_class[h] = use_c1 ? c1 : c0;
                    }
                }

                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);

                // Strict 4-bit class equality semantics
                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    alignas(32) double E8[HAP_NUMBER];
                    for (int h = 0; h < HAP_NUMBER; ++h) {
                        bool match = (expected_class[h] == ss_cond_codes[k]);
                        E8[h] = match ? 1.0 : (M.ed / M.ee);
                    }
                    __m256d emit_lo = _mm256_load_pd(&E8[0]);
                    __m256d emit_hi = _mm256_load_pd(&E8[4]);

                    _sum0 = _mm256_add_pd(_sum0, emit_lo);
                    _sum1 = _mm256_add_pd(_sum1, emit_hi);
                    _mm256_store_pd(&prob[i+0], emit_lo);
                    _mm256_store_pd(&prob[i+4], emit_hi);
                }
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    // Biallelic path - simple AVX2 loop (original algorithm)
    unsigned char amb_code = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                             ? G->Ambiguous[curr_abs_ambiguous] : 0u;
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0;
        g1[h] = HAP_GET(amb_code,h)?1.0:M.ed / M.ee;
    }
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _emit0[2], _emit1[2];
    _emit0[0] = _mm256_loadu_pd(&g0[0]);
    _emit0[1] = _mm256_loadu_pd(&g1[0]);
    _emit1[0] = _mm256_loadu_pd(&g0[4]);
    _emit1[1] = _mm256_loadu_pd(&g1[4]);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0 = _emit0[ah];
        __m256d _prob1 = _emit1[ah];
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i], _prob0);
        _mm256_store_pd(&prob[i+4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::RUN_AMB() {
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
                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                // Build expected class per lane
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

                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);
                __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
                __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
                __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
                _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
                _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
                __m256d _nt = _mm256_set1_pd(nt / probSumT);

                // Strict 4-bit class equality semantics
                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    __m256d _prob0 = _mm256_load_pd(&prob[i]);
                    __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
                    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
                    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);

                    alignas(32) double E8[HAP_NUMBER];
                    for (int h = 0; h < HAP_NUMBER; ++h) {
                        bool match = (expected_class[h] == ss_cond_codes[k]);
                        E8[h] = match ? 1.0 : (M.ed / M.ee);
                    }

                    __m256d emit_lo = _mm256_load_pd(&E8[0]);
                    __m256d emit_hi = _mm256_load_pd(&E8[4]);

                    _prob0 = _mm256_mul_pd(_prob0, emit_lo);
                    _prob1 = _mm256_mul_pd(_prob1, emit_hi);
                    _sum0  = _mm256_add_pd(_sum0, _prob0);
                    _sum1  = _mm256_add_pd(_sum1, _prob1);
                    _mm256_store_pd(&prob[i],   _prob0);
                    _mm256_store_pd(&prob[i+4], _prob1);
                }
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    // Biallelic path
    if (curr_abs_ambiguous < 0 || curr_abs_ambiguous >= (int)G->Ambiguous.size()) {
        fprintf(stderr, "RUN_AMB: Invalid curr_abs_ambiguous index %d, returning\n", curr_abs_ambiguous);
        return;
    }
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0;
        g1[h] = HAP_GET(amb_code,h)?1.0:M.ed / M.ee;
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

inline
void haplotype_segment_double::COLLAPSE_AMB() {
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
                // Load conditioning haplotype codes (cached after first call)
                ss_load_cond_codes(ss, ss_idx);

                // Build expected class per lane
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

                // Unified loop: collapse from probSumK with per-lane emissions
                // BUG FIX #5: Single code path with parameterized emission computation
                __m256d _sum0 = _mm256_set1_pd(0.0);
                __m256d _sum1 = _mm256_set1_pd(0.0);
                __m256d _tFreq0 = _mm256_set1_pd(yt / n_cond_haps);
                __m256d _tFreq1 = _mm256_set1_pd(yt / n_cond_haps);
                __m256d _nt = _mm256_set1_pd(nt / probSumT);

                for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
                    // Transition: collapse from previous segment boundary using donor/column marginals
                    double base = (k < probSumK.size()) ? probSumK[k] : 0.0;
                    __m256d _prob0 = _mm256_set1_pd(base);
                    __m256d _prob1 = _mm256_set1_pd(base);
                    _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
                    _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);

                    // Strict 4-bit class equality semantics
                    alignas(32) double E8[HAP_NUMBER];
                    for (int h = 0; h < HAP_NUMBER; ++h) {
                        bool match = (expected_class[h] == ss_cond_codes[k]);
                        E8[h] = match ? 1.0 : (M.ed / M.ee);
                    }
                    __m256d emit_lo = _mm256_load_pd(&E8[0]);
                    __m256d emit_hi = _mm256_load_pd(&E8[4]);
                    _prob0 = _mm256_mul_pd(_prob0, emit_lo);
                    _prob1 = _mm256_mul_pd(_prob1, emit_hi);
                    
                    _sum0 = _mm256_add_pd(_sum0, _prob0);
                    _sum1 = _mm256_add_pd(_sum1, _prob1);
                    _mm256_store_pd(&prob[i],   _prob0);
                    _mm256_store_pd(&prob[i+4], _prob1);
                }
                
                _mm256_store_pd(&probSumH[0], _sum0);
                _mm256_store_pd(&probSumH[4], _sum1);
                probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
                return;
            }
        }
    }
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed / M.ee:1.0;
        g1[h] = HAP_GET(amb_code,h)?1.0:M.ed / M.ee;
    }
    __m256d _sum0 = _mm256_set1_pd(0.0f);
    __m256d _sum1 = _mm256_set1_pd(0.0f);
    __m256d _tFreq0 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _tFreq1 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
	__m256d _emit0[2], _emit1[2];
	_emit0[0] = _mm256_loadu_pd(&g0[0]);
	_emit0[1] = _mm256_loadu_pd(&g1[0]);
	_emit1[0] = _mm256_loadu_pd(&g0[4]);
	_emit1[1] = _mm256_loadu_pd(&g1[4]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
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

/*******************************************************************************/
/*****************			MISSING GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_MIS() {
	fill(prob.begin(), prob.end(), 1.0/(HAP_NUMBER * n_cond_haps));
	fill(probSumH.begin(), probSumH.end(), 1.0/HAP_NUMBER);
	probSumT = 1.0;
}

inline
void haplotype_segment_double::handle_sibling_bookkeeping(const SiteView& site_view) {
	(void)site_view;
}

inline
void haplotype_segment_double::INIT_SIB(const SiteView& site_view) {
	INIT_MIS();
}

inline
void haplotype_segment_double::RUN_SIB(const SiteView& site_view) {
	handle_sibling_bookkeeping(site_view);
}

inline
void haplotype_segment_double::COLLAPSE_SIB(const SiteView& site_view) {
	handle_sibling_bookkeeping(site_view);
}

inline
void haplotype_segment_double::RUN_MIS() {
	__m256d _sum0 = _mm256_set1_pd(0.0);
	__m256d _sum1 = _mm256_set1_pd(0.0);
	__m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
	__m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
	__m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
	_tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
	_tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
	__m256d _nt = _mm256_set1_pd(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256d _prob0 = _mm256_load_pd(&prob[i]);
		__m256d _prob1 = _mm256_load_pd(&prob[i+4]);
		_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
		_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
		_sum0 = _mm256_add_pd(_sum0, _prob0);
		_sum1 = _mm256_add_pd(_sum1, _prob1);
		_mm256_store_pd(&prob[i], _prob0);
		_mm256_store_pd(&prob[i+4], _prob1);
	}
	_mm256_store_pd(&probSumH[0], _sum0);
	_mm256_store_pd(&probSumH[4], _sum1);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::COLLAPSE_MIS() {
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _tFreq0 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _tFreq1 = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    for (int k = 0, i = 0; k != n_cond_haps; ++k, i += HAP_NUMBER) {
        double base = (k < (int)probSumK.size()) ? probSumK[k] : 0.0;
        __m256d _prob0 = _mm256_set1_pd(base);
        __m256d _prob1 = _mm256_set1_pd(base);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i], _prob0);
        _mm256_store_pd(&prob[i + 4], _prob1);
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************					SUM Ks				************************/
/*******************************************************************************/

inline
void haplotype_segment_double::SUMK() {
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
	}
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS			************************/
/*******************************************************************************/

inline
bool haplotype_segment_double::TRANS_HAP() {
	sumHProbs = 0.0;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], prev_abs_locus);
	nt = 1.0 - yt;
	
	// Guard against zero/near-zero AlphaSumSum (can occur at window boundaries on supersite siblings,
	// or extreme underflow in forward pass). In double precision, use uniform prior as last resort.
	const double prev_total = AlphaSumSum[curr_rel_segment_index - 1];
	const double epsilon = 1e-60;  // Double precision epsilon
	double fact1, uniform_weight;
	bool use_uniform = (prev_total < epsilon);
	
	if (use_uniform) {
		fact1 = 0.0;
		uniform_weight = 1.0 / (n_cond_haps * HAP_NUMBER);
	} else {
		fact1 = nt / prev_total;
		uniform_weight = 0.0;  // Not used, but initialize to avoid uninitialized read
	}
	
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		__m256d _sum0 = _mm256_set1_pd(0.0);
		__m256d _sum1 = _mm256_set1_pd(0.0);
		double fact2;
		if (use_uniform) {
			fact2 = uniform_weight;
		} else {
			fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/prev_total) * yt / n_cond_haps;
		}
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			__m256d _alpha = _mm256_set1_pd(Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1] * fact1 + fact2);
			__m256d _beta0 = _mm256_load_pd(&prob[k*HAP_NUMBER+0]);
			__m256d _beta1 = _mm256_load_pd(&prob[k*HAP_NUMBER+4]);
			_sum0 = _mm256_add_pd(_sum0, _mm256_mul_pd(_alpha, _beta0));
			_sum1 = _mm256_add_pd(_sum1, _mm256_mul_pd(_alpha, _beta1));
		}
		_mm256_store_pd(&HProbs[h1*HAP_NUMBER+0], _sum0);
		_mm256_store_pd(&HProbs[h1*HAP_NUMBER+4], _sum1);
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<double>::min());
}

inline
bool haplotype_segment_double::TRANS_DIP_MULT() {
	sumDProbs= 0.0;
	double scaling = 1.0 / sumHProbs;
	int pd_hits = 0, nd_hits = 0, t_hits = 0;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			pd_hits++;
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					if (pd_hits == 1) nd_hits++;
					DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t_hits++;
					t++;
				}
			}
		}
	}
	return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

inline
bool haplotype_segment_double::TRANS_DIP_ADD() {
	sumDProbs = 0.0;
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
void haplotype_segment_double::IMPUTE(std::vector<float>& missing_probabilities) {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];

    __m256d _alphaSum_lo = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][0]);
    __m256d _alphaSum_hi = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][4]);
    __m256d _ones = _mm256_set1_pd(1.0);
    __m256d _alphaInv_lo = _mm256_div_pd(_ones, _alphaSum_lo);
    __m256d _alphaInv_hi = _mm256_div_pd(_ones, _alphaSum_hi);

    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];

        int target_class = 0;
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) {
                target_class = static_cast<int>(ai) + 1;
                break;
            }
        }

        ss_load_cond_codes(ss, ss_idx);

        __m256d sum_lo[SUPERSITE_MAX_ALTS + 1];
        __m256d sum_hi[SUPERSITE_MAX_ALTS + 1];
        for (int c = 0; c <= (int)ss.var_count; ++c) {
            sum_lo[c] = _mm256_set1_pd(0.0);
            sum_hi[c] = _mm256_set1_pd(0.0);
        }
        __m256d denom_lo = _mm256_set1_pd(0.0);
        __m256d denom_hi = _mm256_set1_pd(0.0);

        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int code = (int)ss_cond_codes[k];
            if (code > (int)ss.var_count) code = 0;

            __m256d _prob_lo = _mm256_load_pd(&prob[i]);
            __m256d _prob_hi = _mm256_load_pd(&prob[i + 4]);
            __m256d _alpha_lo = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i]);
            __m256d _alpha_hi = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i + 4]);
            __m256d term_lo = _mm256_mul_pd(_mm256_mul_pd(_alpha_lo, _alphaInv_lo), _prob_lo);
            __m256d term_hi = _mm256_mul_pd(_mm256_mul_pd(_alpha_hi, _alphaInv_hi), _prob_hi);

            sum_lo[code] = _mm256_add_pd(sum_lo[code], term_lo);
            sum_hi[code] = _mm256_add_pd(sum_hi[code], term_hi);
            denom_lo = _mm256_add_pd(denom_lo, term_lo);
            denom_hi = _mm256_add_pd(denom_hi, term_hi);
        }

        alignas(32) double sum_lo_vals[SUPERSITE_MAX_ALTS + 1][4];
        alignas(32) double sum_hi_vals[SUPERSITE_MAX_ALTS + 1][4];
        alignas(32) double denom_lo_vals[4];
        alignas(32) double denom_hi_vals[4];

        for (int c = 0; c <= (int)ss.var_count; ++c) {
            _mm256_store_pd(sum_lo_vals[c], sum_lo[c]);
            _mm256_store_pd(sum_hi_vals[c], sum_hi[c]);
        }
        _mm256_store_pd(denom_lo_vals, denom_lo);
        _mm256_store_pd(denom_hi_vals, denom_hi);

        for (int lane = 0; lane < 4; ++lane) {
            double denom = denom_lo_vals[lane];
            double numer = sum_lo_vals[target_class][lane];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + lane] = (float)((denom > 0.0) ? (numer / denom) : 0.0);
        }
        for (int lane = 0; lane < 4; ++lane) {
            double denom = denom_hi_vals[lane];
            double numer = sum_hi_vals[target_class][lane];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + 4 + lane] = (float)((denom > 0.0) ? (numer / denom) : 0.0);
        }
    } else {
        __m256d sum_lo[2];
        __m256d sum_hi[2];
        sum_lo[0] = _mm256_set1_pd(0.0);
        sum_lo[1] = _mm256_set1_pd(0.0);
        sum_hi[0] = _mm256_set1_pd(0.0);
        sum_hi[1] = _mm256_set1_pd(0.0);

        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);

            __m256d _prob_lo = _mm256_load_pd(&prob[i]);
            __m256d _prob_hi = _mm256_load_pd(&prob[i + 4]);
            __m256d _alpha_lo = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i]);
            __m256d _alpha_hi = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i + 4]);
            __m256d term_lo = _mm256_mul_pd(_mm256_mul_pd(_alpha_lo, _alphaInv_lo), _prob_lo);
            __m256d term_hi = _mm256_mul_pd(_mm256_mul_pd(_alpha_hi, _alphaInv_hi), _prob_hi);

            sum_lo[ah] = _mm256_add_pd(sum_lo[ah], term_lo);
            sum_hi[ah] = _mm256_add_pd(sum_hi[ah], term_hi);
        }

        alignas(32) double class_lo[2][4];
        alignas(32) double class_hi[2][4];
        _mm256_store_pd(class_lo[0], sum_lo[0]);
        _mm256_store_pd(class_lo[1], sum_lo[1]);
        _mm256_store_pd(class_hi[0], sum_hi[0]);
        _mm256_store_pd(class_hi[1], sum_hi[1]);

        for (int lane = 0; lane < 4; ++lane) {
            double denom = class_lo[0][lane] + class_lo[1][lane];
            double numer = class_lo[1][lane];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + lane] = (float)((denom > 0.0) ? (numer / denom) : 0.0);
        }
        for (int lane = 0; lane < 4; ++lane) {
            double denom = class_hi[0][lane] + class_hi[1][lane];
            double numer = class_hi[1][lane];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + 4 + lane] = (float)((denom > 0.0) ? (numer / denom) : 0.0);
        }
    }
}

// Phase 3: Multivariant imputation for supersites
// Computes P(class_c | Alpha, Beta) for each class c in {0=REF, 1=ALT1, ..., n_alts}
// Writes to SC buffer at ss.class_prob_offset
inline
void haplotype_segment_double::IMPUTE_SUPERSITE_MULTIVARIATE(std::vector < float > & SC, const SuperSite& ss, int ss_idx, int rel_missing_index, const std::vector<uint32_t>* supersite_sc_offset) {
    // Return early if panel codes not available (no supersites built)
    if (panel_codes == nullptr) {
        return;
    }

    ss_load_cond_codes(ss, ss_idx);
    
    const int C = ss.n_classes;  // 1 + n_alts
    
    // Calculate offset - use thread-local offset if available, otherwise fallback to simple calculation
    uint32_t offset;
    if (supersite_sc_offset != nullptr) {
        offset = (*supersite_sc_offset)[ss_idx];  // Use thread-local offset to prevent race conditions
    } else {
        // Fallback for tests: assume simple layout with supersite index * (lanes * classes)
        offset = static_cast<uint32_t>(ss_idx) * HAP_NUMBER * C;
    }
    
    const size_t required = static_cast<size_t>(offset) + static_cast<size_t>(HAP_NUMBER) * static_cast<size_t>(C);
    assert(required <= SC.size());
    
    // Compute 1 / AlphaSum for normalization
    __m256d _alphaSum0 = _mm256_load_pd(&AlphaSumMissing[rel_missing_index][0]);
    __m256d _alphaSum1 = _mm256_load_pd(&AlphaSumMissing[rel_missing_index][4]);
    __m256d _ones = _mm256_set1_pd(1.0);
    _alphaSum0 = _mm256_div_pd(_ones, _alphaSum0);
    _alphaSum1 = _mm256_div_pd(_ones, _alphaSum1);
    
    // Accumulators: per-class, per-lane (lo=lanes 0-3, hi=lanes 4-7)
    __m256d sum_lo[SUPERSITE_MAX_ALTS + 1];
    __m256d sum_hi[SUPERSITE_MAX_ALTS + 1];
    for (int c = 0; c < C; ++c) {
        sum_lo[c] = _mm256_set1_pd(0.0);
        sum_hi[c] = _mm256_set1_pd(0.0);
    }
    
    __m256d denom_lo = _mm256_set1_pd(0.0);
    __m256d denom_hi = _mm256_set1_pd(0.0);
    
    // Accumulate: for each conditioning haplotype k, add its contribution to class bucket
    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        int code = (int)ss_cond_codes[k];  // 0=REF, 1..n_alts
        if (code >= C) code = 0;  // Safety: invalid code -> REF
        
        // Load Alpha (forward) and Beta (backward prob)
        __m256d _alpha0 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i]);
        __m256d _alpha1 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+4]);
        __m256d _beta0  = _mm256_load_pd(&prob[i]);
        __m256d _beta1  = _mm256_load_pd(&prob[i+4]);
        
        // Posterior contribution: Alpha * Beta / AlphaSum
        __m256d _term0 = _mm256_mul_pd(_mm256_mul_pd(_alpha0, _alphaSum0), _beta0);
        __m256d _term1 = _mm256_mul_pd(_mm256_mul_pd(_alpha1, _alphaSum1), _beta1);
        
        // Add to class bucket
        sum_lo[code] = _mm256_add_pd(sum_lo[code], _term0);
        sum_hi[code] = _mm256_add_pd(sum_hi[code], _term1);
        
        // Accumulate denominator (total across all classes)
        denom_lo = _mm256_add_pd(denom_lo, _term0);
        denom_hi = _mm256_add_pd(denom_hi, _term1);
    }
    
    // Normalize and write to SC buffer
    // Layout: SC[offset + h*C + c] = P(class c | lane h)
    alignas(32) double lo_buf[4], hi_buf[4], dlo[4], dhi[4];
    _mm256_store_pd(dlo, denom_lo);
    _mm256_store_pd(dhi, denom_hi);
    
    for (int c = 0; c < C; ++c) {
        _mm256_store_pd(lo_buf, sum_lo[c]);
        _mm256_store_pd(hi_buf, sum_hi[c]);
        
        // Lanes 0-3 (lo)
        for (int h = 0; h < 4; ++h) {
            float p = (dlo[h] > 0.0) ? (float)(lo_buf[h] / dlo[h]) : (1.0f / C);  // Fallback to uniform
            SC[offset + h * C + c] = p;
        }
        
        // Lanes 4-7 (hi)
        for (int h = 4; h < HAP_NUMBER; ++h) {
            int hidx = h - 4;
            float p = (dhi[hidx] > 0.0) ? (float)(hi_buf[hidx] / dhi[hidx]) : (1.0f / C);
            SC[offset + h * C + c] = p;
        }
    }

}

#endif

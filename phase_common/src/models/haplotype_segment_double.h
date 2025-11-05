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
#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>
#include <models/site_emission_types.h>
#include <models/super_site_macros.h>

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
	int curr_segment_index;
	int curr_segment_locus;
	int curr_abs_locus;
	int prev_abs_locus;
	int curr_rel_locus;
	int curr_rel_locus_offset;
	int curr_abs_ambiguous;
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
	const std::vector<bool>* is_super_site;
	const std::vector<int>* locus_to_super_idx;
	const uint8_t* panel_codes;
	const std::vector<int>* super_site_var_index;
	const std::vector<unsigned int>* cond_idx;
	aligned_vector32<uint8_t> ss_cond_codes;
	aligned_vector32<double> ss_emissions;
	aligned_vector32<double> ss_emissions_h1;
	std::vector<bool> ss_cached;  // Phase 3: cache flags per supersite

	// Anchor MIS mapping: record rel-missing index per locus in window
	std::vector<int> missing_index_by_locus;

	// EMISSION HELPERS
	MatchMask init_match_mask;

	//INLINED AND UNROLLED ROUTINES
	void INIT_HOM();
	void INIT_AMB();
	void INIT_MIS();
	void INIT_FROM_MASK(const MatchMask& mask, double mismatch_penalty);
	bool RUN_HOM(char);
	void RUN_AMB();
	void RUN_MIS();
	void COLLAPSE_HOM();
	void COLLAPSE_AMB();
	void COLLAPSE_MIS();
	
	// Supersite-specific helper functions
	void SS_INIT_HOM();
	void SS_INIT_AMB();
	void SS_INIT_MIS();
	bool SS_RUN_HOM();
	void SS_RUN_AMB();
	void SS_RUN_MIS();
	void SS_COLLAPSE_HOM();
	void SS_COLLAPSE_AMB();
	void SS_COLLAPSE_MIS();
	
	// Caching helper
	void ss_load_cond_codes(const SuperSite& ss, int ss_idx);
	
	void SUMK();
	void IMPUTE(std::vector < float > & );
	void IMPUTE_SUPERSITE_MULTIVARIATE(std::vector < float > & SC, const SuperSite& ss, int ss_idx, int rel_missing_index);  // Phase 3
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(std::vector < double > & );
	int SET_OTHER_TRANS(std::vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment_double(genotype *, bitmatrix &, std::vector < unsigned int > &, window &, hmm_parameters &,
		const std::vector<SuperSite>* _super_sites = nullptr,
		const std::vector<bool>* _is_super_site = nullptr,
		const std::vector<int>* _locus_to_super_idx = nullptr,
		const uint8_t* _panel_codes = nullptr,
		const std::vector<int>* _super_site_var_index = nullptr);
	~haplotype_segment_double();

	//void fetch();
	void forward();
	int backward(std::vector < double > &, std::vector < float > &, 
	            std::vector < float > * SC = nullptr,
	            const std::vector < bool > * anchor_has_missing = nullptr);  // Phase 3: optional supersite posteriors
	
	// Supersite cache management (Phase 3)
	// Cache is per-segment and automatically reset when new segment created.
	// Segments are local variables in phaser_algorithm.cpp phaseWindow(),
	// so cache resets automatically between windows.
	void clear_supersite_cache() {
		std::fill(ss_cached.begin(), ss_cached.end(), false);
	}
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
    // Return if already cached
    if (ss_cached[ss_idx]) return;
    
    // Unpack and cache conditioning haplotype codes
    ss_cond_codes.resize(n_cond_haps);
    for (int k = 0; k < (int)n_cond_haps; ++k) {
        unsigned int gh = (*cond_idx)[k];
        ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
    }
    
    // Mark as cached
    ss_cached[ss_idx] = true;
}

// Note: pack_expected_codes_pd() helper removed - was causing half-lane split antipattern
// Now using proper per-lane emission logic throughout SS_*_AMB() functions

inline
void haplotype_segment_double::SS_INIT_HOM() {
    // Assumes: ss_idx >= 0, at anchor, classified as HOM or AMB (c0==c1)
    // Uses member: super_sites, locus_to_super_idx, curr_abs_locus, super_site_var_index
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
    uint8_t sample_code = c0;  // For HOM, c0==c1
    
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
}

inline
void haplotype_segment_double::SS_INIT_AMB() {
    // Assumes: ss_idx >= 0, at anchor, classified as AMB (c0 != c1)
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);

    // Load conditioning haplotype codes (cached after first call)
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

    __m256d match_d = _mm256_set1_pd(1.0);
    __m256d mis_d   = _mm256_set1_pd(M.ed/M.ee);

    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);

    // Anchor ALT code for this locus
    int anchor_code = 0; 
    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
        if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { anchor_code = (int)ai + 1; break; }
    }

    if (M.ss_anchor_split_emissions) {
        // Split semantics: per-lane expected ALT at anchor vs donor carries anchor ALT
        alignas(32) int exp_is_alt[HAP_NUMBER];
        for (int h = 0; h < HAP_NUMBER; ++h) exp_is_alt[h] = (expected_class[h] == anchor_code) ? 1 : 0;
        
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int donor_is_alt = ((int)ss_cond_codes[k] == anchor_code) ? 1 : 0;
            
            // Create per-lane emissions (process as two 4-lane halves for double precision)
            alignas(32) double E8[HAP_NUMBER];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                E8[h] = (donor_is_alt == exp_is_alt[h]) ? 1.0 : (M.ed/M.ee);
            }
            
            // Load as two 256-bit double vectors
            __m256d emit_lo = _mm256_load_pd(&E8[0]);
            __m256d emit_hi = _mm256_load_pd(&E8[4]);
            
            _sum0 = _mm256_add_pd(_sum0, emit_lo);
            _sum1 = _mm256_add_pd(_sum1, emit_hi);
            _mm256_store_pd(&prob[i+0], emit_lo);
            _mm256_store_pd(&prob[i+4], emit_hi);
        }
    } else {
        // Strict 4-bit equality semantics: per-lane comparison
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int dc = (int)ss_cond_codes[k];
            
            // Create per-lane emissions (process as two 4-lane halves for double precision)
            alignas(32) double E8[HAP_NUMBER];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                E8[h] = (dc == (int)expected_class[h]) ? 1.0 : (M.ed/M.ee);
            }
            
            // Load as two 256-bit double vectors
            __m256d emit_lo = _mm256_load_pd(&E8[0]);
            __m256d emit_hi = _mm256_load_pd(&E8[4]);
            
            _sum0 = _mm256_add_pd(_sum0, emit_lo);
            _sum1 = _mm256_add_pd(_sum1, emit_hi);
            _mm256_store_pd(&prob[i+0], emit_lo);
            _mm256_store_pd(&prob[i+4], emit_hi);
        }
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::SS_INIT_MIS() {
    // Missing data initialization - same as biallelic
    fill(prob.begin(), prob.end(), 1.0/(HAP_NUMBER * n_cond_haps));
    fill(probSumH.begin(), probSumH.end(), 1.0/HAP_NUMBER);
    probSumT = 1.0;
}

inline
bool haplotype_segment_double::SS_RUN_HOM() {
    // Assumes: ss_idx >= 0, at anchor, classified as HOM or AMB (c0==c1)
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
    uint8_t sample_code = c0;
    
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
        _prob0 = _mm256_mul_pd(_prob0, _emit);
        _prob1 = _mm256_mul_pd(_prob1, _emit);
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

inline
void haplotype_segment_double::SS_RUN_AMB() {
    // Assumes: ss_idx >= 0, at anchor, classified as AMB (c0 != c1)
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);

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

    __m256d match_d = _mm256_set1_pd(1.0);
    __m256d mis_d   = _mm256_set1_pd(M.ed/M.ee);

    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
    __m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
    __m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
    _tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
    _tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);

    // Anchor ALT code for this locus
    int anchor_code = 0; 
    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
        if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { anchor_code = (int)ai + 1; break; }
    }

    if (M.ss_anchor_split_emissions) {
        // Split semantics: per-lane expected ALT at anchor vs donor carries anchor ALT
        alignas(32) int exp_is_alt[HAP_NUMBER];
        for (int h = 0; h < HAP_NUMBER; ++h) exp_is_alt[h] = (expected_class[h] == anchor_code) ? 1 : 0;
        
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256d _prob0 = _mm256_load_pd(&prob[i]);
            __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
            _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
            _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
            
            int donor_is_alt = ((int)ss_cond_codes[k] == anchor_code) ? 1 : 0;
            
            // Create per-lane emissions (process as two 4-lane halves for double precision)
            alignas(32) double E8[HAP_NUMBER];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                E8[h] = (donor_is_alt == exp_is_alt[h]) ? 1.0 : (M.ed/M.ee);
            }
            
            // Load as two 256-bit double vectors
            __m256d emit_lo = _mm256_load_pd(&E8[0]);
            __m256d emit_hi = _mm256_load_pd(&E8[4]);
            
            _prob0 = _mm256_mul_pd(_prob0, emit_lo);
            _prob1 = _mm256_mul_pd(_prob1, emit_hi);
            _sum0 = _mm256_add_pd(_sum0, _prob0);
            _sum1 = _mm256_add_pd(_sum1, _prob1);
            _mm256_store_pd(&prob[i], _prob0);
            _mm256_store_pd(&prob[i+4], _prob1);
        }
    } else {
        // Strict 4-bit equality semantics: per-lane comparison
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256d _prob0 = _mm256_load_pd(&prob[i]);
            __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
            _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
            _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);

            int dc = (int)ss_cond_codes[k];
            
            // Create per-lane emissions (process as two 4-lane halves for double precision)
            alignas(32) double E8[HAP_NUMBER];
            for (int h = 0; h < HAP_NUMBER; ++h) {
                E8[h] = (dc == (int)expected_class[h]) ? 1.0 : (M.ed/M.ee);
            }
            
            // Load as two 256-bit double vectors
            __m256d emit_lo = _mm256_load_pd(&E8[0]);
            __m256d emit_hi = _mm256_load_pd(&E8[4]);

            _prob0 = _mm256_mul_pd(_prob0, emit_lo);
            _prob1 = _mm256_mul_pd(_prob1, emit_hi);
            _sum0  = _mm256_add_pd(_sum0, _prob0);
            _sum1  = _mm256_add_pd(_sum1, _prob1);
            _mm256_store_pd(&prob[i],   _prob0);
            _mm256_store_pd(&prob[i+4], _prob1);
        }
    }
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::SS_RUN_MIS() {
    // Missing data run - same as biallelic
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
void haplotype_segment_double::SS_COLLAPSE_HOM() {
    // Assumes: ss_idx >= 0, at anchor, classified as HOM or AMB (c0==c1)
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
    uint8_t sample_code = c0;
    
    // Load conditioning haplotype codes (cached after first call)
    ss_load_cond_codes(ss, ss_idx);
    
    precomputeSuperSiteEmissions_AVX2(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0, M.ed / M.ee, ss_emissions);
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        __m256d _emit = _mm256_set1_pd(ss_emissions[k]);
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
        _prob0 = _mm256_mul_pd(_prob0, _emit);
        _prob1 = _mm256_mul_pd(_prob1, _emit);
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
void haplotype_segment_double::SS_COLLAPSE_AMB() {
    // Assumes: ss_idx >= 0, at anchor, classified as AMB (c0 != c1)
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    const SuperSite& ss = (*super_sites)[ss_idx];
    uint8_t c0, c1;
    SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);

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

    // Precompute expected vector based on mode
    __m256d match_d = _mm256_set1_pd(1.0);
    __m256d mis_d   = _mm256_set1_pd(M.ed/M.ee);
    alignas(32) int exp_arr[HAP_NUMBER];
    
    if (M.ss_anchor_split_emissions) {
        // Split semantics: compare binary (is anchor ALT?) for both donor and expected
        int anchor_code = 0; 
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { 
                anchor_code = (int)ai + 1; 
                break; 
            }
        }
        for (int h = 0; h < HAP_NUMBER; ++h) 
            exp_arr[h] = (expected_class[h] == anchor_code) ? 1 : 0;
    } else {
        // Strict 4-bit class equality semantics
        for (int h = 0; h < HAP_NUMBER; ++h) 
            exp_arr[h] = (int)expected_class[h];
    }

    // Unified loop: collapse from probSumK with per-lane emissions
    // BUG FIX #5: Single code path with parameterized emission computation
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);

    for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
        // Transition: collapse from previous segment boundary
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
        
        // Emission: compute donor code based on mode, then compare to expected per lane
        // BUG #6 DOCUMENTED: Supersite uses per-lane array computation (required for double precision)
        // vs. biallelic inline conditional (optimized for binary alleles)
        // Both implement: emit[h] = (donor_matches_expected[h]) ? 1.0 : (ed/ee)
        int donor_code;
        if (M.ss_anchor_split_emissions) {
            int anchor_code = 0; 
            for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
                if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { 
                    anchor_code = (int)ai + 1; 
                    break; 
                }
            }
            donor_code = ((int)ss_cond_codes[k] == anchor_code) ? 1 : 0;
        } else {
            donor_code = (int)ss_cond_codes[k];
        }
        
        // Build per-lane emissions (process as two 4-lane halves for double precision)
        alignas(32) double E8[HAP_NUMBER];
        for (int h = 0; h < HAP_NUMBER; ++h) {
            E8[h] = (donor_code == exp_arr[h]) ? 1.0 : (M.ed/M.ee);
        }
        
        // Load as two 256-bit double vectors and apply
        __m256d emit_lo = _mm256_load_pd(&E8[0]);
        __m256d emit_hi = _mm256_load_pd(&E8[4]);
        _prob0 = _mm256_mul_pd(_prob0, emit_lo);
        _prob1 = _mm256_mul_pd(_prob1, emit_hi);
        
        // Accumulate and store
        _sum0 = _mm256_add_pd(_sum0, _prob0);
        _sum1 = _mm256_add_pd(_sum1, _prob1);
        _mm256_store_pd(&prob[i],   _prob0);
        _mm256_store_pd(&prob[i+4], _prob1);
    }
    
    _mm256_store_pd(&probSumH[0], _sum0);
    _mm256_store_pd(&probSumH[4], _sum1);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_double::SS_COLLAPSE_MIS() {
    // Missing data collapse - same as biallelic
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    __m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
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
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment_double::INIT_FROM_MASK(const MatchMask& mask, double mismatch_penalty) {
	const __m256d match_vec = _mm256_set1_pd(1.0);
	const __m256d mismatch_vec = _mm256_set1_pd(mismatch_penalty);
	__m256d sum0 = _mm256_set1_pd(0.0);
	__m256d sum1 = _mm256_set1_pd(0.0);
	const __m256i zero = _mm256_setzero_si256();
	const uint8_t* mask_data = mask.by_donor_lane.data();

	for (unsigned int k = 0, idx = 0; k < n_cond_haps; ++k, idx += HAP_NUMBER) {
		__m128i mask_u8 = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(mask_data + idx));
		__m256i mask_epi32 = _mm256_cvtepu8_epi32(mask_u8);
		__m128i mask_lo_128 = _mm256_castsi256_si128(mask_epi32);
		__m128i mask_hi_128 = _mm256_extracti128_si256(mask_epi32, 1);
		__m256i mask_lo_epi64 = _mm256_cvtepi32_epi64(mask_lo_128);
		__m256i mask_hi_epi64 = _mm256_cvtepi32_epi64(mask_hi_128);
		__m256i sign_lo = _mm256_cmpgt_epi64(mask_lo_epi64, zero);
		__m256i sign_hi = _mm256_cmpgt_epi64(mask_hi_epi64, zero);
		__m256d emit0 = _mm256_blendv_pd(mismatch_vec, match_vec, _mm256_castsi256_pd(sign_lo));
		__m256d emit1 = _mm256_blendv_pd(mismatch_vec, match_vec, _mm256_castsi256_pd(sign_hi));
		sum0 = _mm256_add_pd(sum0, emit0);
		sum1 = _mm256_add_pd(sum1, emit1);
		_mm256_store_pd(&prob[idx], emit0);
		_mm256_store_pd(&prob[idx + 4], emit1);
	}

	_mm256_store_pd(&probSumH[0], sum0);
	_mm256_store_pd(&probSumH[4], sum1);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}


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
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_INIT_HOM(); return;
            case SSClass::AMB: 
                // For HOM entry point with c0==c1, route to HOM helper
                if (c0 == c1) { SS_INIT_HOM(); return; }
                SS_INIT_AMB(); return;
        }
    }
    // Default biallelic path
    const bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    init_match_mask.resize(static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER);
    uint8_t* mask_ptr = init_match_mask.by_donor_lane.data();
    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const bool ah = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
        const bool is_match = (ag == ah);
        const uint8_t value = is_match ? MatchMask::kMatch : MatchMask::kMismatch;
        std::fill_n(mask_ptr + static_cast<std::size_t>(k) * HAP_NUMBER, HAP_NUMBER, value);
        if (is_match) {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                init_match_mask.any_match_lane[h] = true;
            }
        }
    }
    INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
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
            // Sibling: treat as uninformative locus (BUG FIX #2)
            RUN_MIS();
            return true;
        }
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return true;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: return SS_RUN_HOM();
            case SSClass::AMB: 
                if (c0 == c1) return SS_RUN_HOM();
                SS_RUN_AMB(); return true;
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
		__m256d _mismatch = _mm256_set1_pd(M.ed/M.ee);
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
            // Sibling: treat as uninformative locus (BUG FIX #2)
            COLLAPSE_MIS();
            return;
        }
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_COLLAPSE_HOM(); return;
            case SSClass::AMB: 
                if (c0 == c1) { SS_COLLAPSE_HOM(); return; }
                SS_COLLAPSE_AMB(); return;
        }
    }
    // Biallelic path
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);
    // BUG #4 EXPERIMENTAL: Test normalization behavior
    static const char* norm_env = std::getenv("SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION");
    static const bool use_normalization = (norm_env && norm_env[0] != '\0' && norm_env[0] != '0');
    __m256d _tFreq = use_normalization
        ? _mm256_set1_pd((yt * probSumT) / n_cond_haps)
        : _mm256_set1_pd(yt / n_cond_haps);
    __m256d _nt = _mm256_set1_pd(nt / probSumT);
    __m256d _mismatch = _mm256_set1_pd(M.ed/M.ee);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256d _prob0 = _mm256_set1_pd(probSumK[k]);
        __m256d _prob1 = _mm256_set1_pd(probSumK[k]);
        _prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
        _prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
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
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: INIT_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_INIT_HOM(); return;
            case SSClass::AMB: SS_INIT_AMB(); return;
        }
    }
    // Biallelic path
    unsigned char amb_code = (curr_abs_ambiguous >= ambiguous_first && curr_abs_ambiguous <= ambiguous_last)
                             ? G->Ambiguous[curr_abs_ambiguous] : 0u;
    init_match_mask.resize(static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER);
    uint8_t* mask_ptr = init_match_mask.by_donor_lane.data();
    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const bool donor_alt = Hvar.get(curr_rel_locus + curr_rel_locus_offset, k);
        const uint8_t donor_code = donor_alt ? 1u : 0u;
        const std::size_t base = static_cast<std::size_t>(k) * HAP_NUMBER;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool lane_wants_alt = ((amb_code >> h) & 1U) != 0;
            const uint8_t expected = lane_wants_alt ? 1u : 0u;
            const bool match = (donor_code == expected);
            mask_ptr[base + h] = match ? MatchMask::kMatch : MatchMask::kMismatch;
            init_match_mask.any_match_lane[h] = init_match_mask.any_match_lane[h] || match;
        }
    }
    INIT_FROM_MASK(init_match_mask, M.ed/M.ee);
}

inline
void haplotype_segment_double::RUN_AMB() {
    // Supersite dispatcher
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: treat as uninformative locus (BUG FIX #2)
            RUN_MIS();
            return;
        }
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: RUN_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_RUN_HOM(); return;
            case SSClass::AMB: SS_RUN_AMB(); return;
        }
    }
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0;
        g1[h] = HAP_GET(amb_code,h)?1.0:M.ed/M.ee;
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
            // Sibling: treat as uninformative locus (BUG FIX #2)
            COLLAPSE_MIS();
            return;
        }
        
        // Classify and dispatch to SS_* helpers
        uint8_t c0, c1;
        SSClass cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        switch (cls) {
            case SSClass::MIS: COLLAPSE_MIS(); return;  // BUG FIX #1: Use biallelic MIS
            case SSClass::HOM: SS_COLLAPSE_HOM(); return;
            case SSClass::AMB: SS_COLLAPSE_AMB(); return;
        }
    }
    // Biallelic path
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0;
        g1[h] = HAP_GET(amb_code,h)?1.0:M.ed/M.ee;
    }
	__m256d _sum0 = _mm256_set1_pd(0.0f);
	__m256d _sum1 = _mm256_set1_pd(0.0f);
	__m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps);
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
		_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
		_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
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
	// BUG #4 EXPERIMENTAL: Test normalization behavior
	static const char* norm_env = std::getenv("SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION");
	static const bool use_normalization = (norm_env && norm_env[0] != '\0' && norm_env[0] != '0');
	__m256d _tFreq = use_normalization
	    ? _mm256_set1_pd((yt * probSumT) / n_cond_haps)
	    : _mm256_set1_pd(yt / n_cond_haps);
	__m256d _nt = _mm256_set1_pd(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256d _prob0 = _mm256_set1_pd(probSumK[k]);
		__m256d _prob1 = _mm256_set1_pd(probSumK[k]);
		_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq);
		_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq);
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
	double fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		__m256d _sum0 = _mm256_set1_pd(0.0);
		__m256d _sum1 = _mm256_set1_pd(0.0);
		double fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * yt / n_cond_haps;
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
void haplotype_segment_double::IMPUTE(std::vector < float > & missing_probabilities) {
    // Supersite-aware imputation: classify cond hap codes as REF(0) vs ALT(>0)
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];

    __m256d _sum0 = _mm256_set1_pd(0.0);
    __m256d _sum1 = _mm256_set1_pd(0.0);

    __m256d _sumA0[2], _sumA1[2];
    _sumA0[0] = _mm256_set1_pd(0.0);
    _sumA0[1] = _mm256_set1_pd(0.0);
    _sumA1[0] = _mm256_set1_pd(0.0);
    _sumA1[1] = _mm256_set1_pd(0.0);

    __m256d _alphaSum0 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][0]);
    __m256d _alphaSum1 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][4]);
    __m256d _ones = _mm256_set1_pd(1.0);
    _alphaSum0 = _mm256_div_pd(_ones, _alphaSum0);
    _alphaSum1 = _mm256_div_pd(_ones, _alphaSum1);

    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Identify target class for current missing variant within supersite
        int target_class = 0; // 0 = REF, 1..n_alts = ALT index
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) {
                target_class = (int)ai + 1;
                break;
            }
        }
        // Unpack conditioning hap codes
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        // Initialize per-class accumulators (up to SUPERSITE_MAX_ALTS+1)
        __m256d sum_lo_classes[SUPERSITE_MAX_ALTS + 1];
        __m256d sum_hi_classes[SUPERSITE_MAX_ALTS + 1];
        for (int c = 0; c <= (int)ss.var_count; ++c) {
            sum_lo_classes[c] = _mm256_set1_pd(0.0);
            sum_hi_classes[c] = _mm256_set1_pd(0.0);
        }
        __m256d denom_lo = _mm256_set1_pd(0.0);
        __m256d denom_hi = _mm256_set1_pd(0.0);
        // Accumulate per conditioning haplotype into class buckets
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int code = (int)ss_cond_codes[k]; // 0..n_alts
            if (code > (int)ss.var_count) code = 0; // safety
            __m256d _prob0 = _mm256_load_pd(&prob[i]);
            __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
            __m256d _alpha0 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+0]);
            __m256d _alpha1 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+4]);
            _sum0 = _mm256_mul_pd(_mm256_mul_pd(_alpha0, _alphaSum0), _prob0);
            _sum1 = _mm256_mul_pd(_mm256_mul_pd(_alpha1, _alphaSum1), _prob1);
            sum_lo_classes[code] = _mm256_add_pd(sum_lo_classes[code], _sum0);
            sum_hi_classes[code] = _mm256_add_pd(sum_hi_classes[code], _sum1);
            denom_lo = _mm256_add_pd(denom_lo, _sum0);
            denom_hi = _mm256_add_pd(denom_hi, _sum1);
        }
        // Compute per-lane posterior for target class = ALT for this specific split variant
        __m256d p_lo = sum_lo_classes[target_class];
        __m256d p_hi = sum_hi_classes[target_class];
        // Store results
        alignas(32) double plo[4], phi[4], dlo[4], dhi[4];
        _mm256_store_pd(plo, p_lo);
        _mm256_store_pd(phi, p_hi);
        _mm256_store_pd(dlo, denom_lo);
        _mm256_store_pd(dhi, denom_hi);
        for (int h = 0; h < 4; ++h) {
            double denom = dlo[h];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = (float)((denom > 0.0) ? (plo[h] / denom) : 0.0);
        }
        for (int h = 0; h < 4; ++h) {
            double denom = dhi[h];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + (4 + h)] = (float)((denom > 0.0) ? (phi[h] / denom) : 0.0);
        }
        return;
    } else {
        // Biallelic path
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256d _prob0 = _mm256_load_pd(&prob[i]);
            __m256d _prob1 = _mm256_load_pd(&prob[i+4]);
            __m256d _alpha0 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+0]);
            __m256d _alpha1 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+4]);
            _sum0 = _mm256_mul_pd(_mm256_mul_pd(_alpha0, _alphaSum0), _prob0);
            _sum1 = _mm256_mul_pd(_mm256_mul_pd(_alpha1, _alphaSum1), _prob1);
            _sumA0[ah] = _mm256_add_pd(_sumA0[ah], _sum0);
            _sumA1[ah] = _mm256_add_pd(_sumA1[ah], _sum1);
        }
    }

    double * prob0 = (double*)&_sumA0[0];
    double * prob1 = (double*)&_sumA0[1];
    for (int h = 0 ; h < 4 ; h ++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h]+prob1[h]);
    prob0 = (double*)&_sumA1[0];
    prob1 = (double*)&_sumA1[1];
    for (int h = 4 ; h < HAP_NUMBER ; h ++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h]+prob1[h]);
}

// Phase 3: Multivariant imputation for supersites
// Computes P(class_c | Alpha, Beta) for each class c in {0=REF, 1=ALT1, ..., n_alts}
// Writes to SC buffer at ss.class_prob_offset
inline
void haplotype_segment_double::IMPUTE_SUPERSITE_MULTIVARIATE(std::vector < float > & SC, const SuperSite& ss, int ss_idx, int rel_missing_index) {
    // Unpack conditioning haplotype codes once
    for (int k = 0; k < (int)n_cond_haps; ++k) {
        unsigned int gh = (*cond_idx)[k];
        ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
    }
    
    const int C = ss.n_classes;  // 1 + n_alts
    const uint32_t offset = ss.class_prob_offset;
    
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

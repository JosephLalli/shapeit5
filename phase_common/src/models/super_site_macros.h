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

#ifndef _SUPER_SITE_MACROS_H
#define _SUPER_SITE_MACROS_H

#include <models/super_site_accessor.h>
#include <models/super_site_emissions.h>

// Note: These macros are placeholders that will be integrated into the phasing loop
// They assume the existence of certain variables from the HMM context (prob, probSumH, etc.)

// ============================================================================
// INIT_HOM_SUPERSITE - Initialize at first position in window
// ============================================================================

#define INIT_HOM_SUPERSITE() \
do { \
	const SuperSite& ss = super_sites[curr_abs_locus]; \
	uint8_t sample_code = getSampleSuperSiteAlleleCode(G->SuperSiteGenotypes, curr_abs_locus, 0); \
	\
	if (sample_code != SUPERSITE_CODE_MISSING) { \
		/* Unpack all conditioning haplotype codes for this super-site */ \
		aligned_vector32<uint8_t> cond_codes(n_cond_haps); \
		unpackSuperSiteCodesBatch( \
			G->SuperSitePanelCodes, \
			ss.panel_offset, \
			0, \
			n_cond_haps, \
			cond_codes.data()); \
		\
		/* Precompute emissions using AVX2 */ \
		aligned_vector32<double> emissions(n_cond_haps); \
		double match_prob = 1.0; \
		double mismatch_prob = M.ed / M.ee; \
		precomputeSuperSiteEmissions_AVX2( \
			cond_codes.data(), \
			n_cond_haps, \
			sample_code, \
			match_prob, \
			mismatch_prob, \
			emissions); \
		\
		/* Standard initialization FMA loop with precomputed emissions */ \
		__m256d _sum0 = _mm256_set1_pd(0.0); \
		__m256d _sum1 = _mm256_set1_pd(0.0); \
		__m256d _f0 = _mm256_set1_pd(1.0 / n_cond_haps); \
		\
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) { \
			__m256d _emit = _mm256_set1_pd(emissions[k]); \
			__m256d _prob0 = _mm256_mul_pd(_emit, _f0); \
			__m256d _prob1 = _mm256_mul_pd(_emit, _f0); \
			_sum0 = _mm256_add_pd(_sum0, _prob0); \
			_sum1 = _mm256_add_pd(_sum1, _prob1); \
			_mm256_store_pd(&prob[i], _prob0); \
			_mm256_store_pd(&prob[i+4], _prob1); \
		} \
		\
		_mm256_store_pd(&probSumH[0], _sum0); \
		_mm256_store_pd(&probSumH[4], _sum1); \
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7]; \
	} else { \
		/* Missing: use standard INIT_MIS */ \
		INIT_MIS(); \
	} \
} while(0)

// ============================================================================
// RUN_HOM_SUPERSITE - Standard update during segment
// ============================================================================

#define RUN_HOM_SUPERSITE() \
do { \
	const SuperSite& ss = super_sites[curr_abs_locus]; \
	uint8_t sample_code = getSampleSuperSiteAlleleCode(G->SuperSiteGenotypes, curr_abs_locus, 0); \
	\
	if (sample_code != SUPERSITE_CODE_MISSING) { \
		/* Unpack conditioning haplotype codes */ \
		aligned_vector32<uint8_t> cond_codes(n_cond_haps); \
		unpackSuperSiteCodesBatch( \
			G->SuperSitePanelCodes, \
			ss.panel_offset, \
			0, \
			n_cond_haps, \
			cond_codes.data()); \
		\
		/* Precompute emissions */ \
		aligned_vector32<double> emissions(n_cond_haps); \
		double match_prob = 1.0; \
		double mismatch_prob = M.ed / M.ee; \
		precomputeSuperSiteEmissions_AVX2( \
			cond_codes.data(), \
			n_cond_haps, \
			sample_code, \
			match_prob, \
			mismatch_prob, \
			emissions); \
		\
		/* Standard HMM update with emissions */ \
		__m256d _sum0 = _mm256_set1_pd(0.0); \
		__m256d _sum1 = _mm256_set1_pd(0.0); \
		__m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT)); \
		__m256d _tFreq0 = _mm256_load_pd(&probSumH[0]); \
		__m256d _tFreq1 = _mm256_load_pd(&probSumH[4]); \
		_tFreq0 = _mm256_mul_pd(_tFreq0, _factor); \
		_tFreq1 = _mm256_mul_pd(_tFreq1, _factor); \
		__m256d _nt = _mm256_set1_pd(nt / probSumT); \
		\
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) { \
			__m256d _emit = _mm256_set1_pd(emissions[k]); \
			__m256d _prob0 = _mm256_load_pd(&prob[i]); \
			__m256d _prob1 = _mm256_load_pd(&prob[i+4]); \
			_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0); \
			_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1); \
			_prob0 = _mm256_mul_pd(_prob0, _emit); \
			_prob1 = _mm256_mul_pd(_prob1, _emit); \
			_sum0 = _mm256_add_pd(_sum0, _prob0); \
			_sum1 = _mm256_add_pd(_sum1, _prob1); \
			_mm256_store_pd(&prob[i], _prob0); \
			_mm256_store_pd(&prob[i+4], _prob1); \
		} \
		\
		_mm256_store_pd(&probSumH[0], _sum0); \
		_mm256_store_pd(&probSumH[4], _sum1); \
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7]; \
	} else { \
		RUN_MIS(); \
	} \
} while(0)

// ============================================================================
// COLLAPSE_HOM_SUPERSITE - Segment boundary update
// ============================================================================

#define COLLAPSE_HOM_SUPERSITE() \
do { \
	const SuperSite& ss = super_sites[curr_abs_locus]; \
	uint8_t sample_code = getSampleSuperSiteAlleleCode(G->SuperSiteGenotypes, curr_abs_locus, 0); \
	\
	if (sample_code != SUPERSITE_CODE_MISSING) { \
		/* Unpack conditioning haplotype codes */ \
		aligned_vector32<uint8_t> cond_codes(n_cond_haps); \
		unpackSuperSiteCodesBatch( \
			G->SuperSitePanelCodes, \
			ss.panel_offset, \
			0, \
			n_cond_haps, \
			cond_codes.data()); \
		\
		/* Precompute emissions */ \
		aligned_vector32<double> emissions(n_cond_haps); \
		double match_prob = 1.0; \
		double mismatch_prob = M.ed / M.ee; \
		precomputeSuperSiteEmissions_AVX2( \
			cond_codes.data(), \
			n_cond_haps, \
			sample_code, \
			match_prob, \
			mismatch_prob, \
			emissions); \
		\
		/* Collapse update using per-haplotype sums */ \
		__m256d _sum0 = _mm256_set1_pd(0.0); \
		__m256d _sum1 = _mm256_set1_pd(0.0); \
		__m256d _tFreq = _mm256_set1_pd(yt / n_cond_haps); \
		__m256d _nt = _mm256_set1_pd(nt / probSumT); \
		\
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) { \
			__m256d _emit = _mm256_set1_pd(emissions[k]); \
			__m256d _prob0 = _mm256_set1_pd(probSumK[k]); \
			__m256d _prob1 = _mm256_set1_pd(probSumK[k]); \
			_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq); \
			_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq); \
			_prob0 = _mm256_mul_pd(_prob0, _emit); \
			_prob1 = _mm256_mul_pd(_prob1, _emit); \
			_sum0 = _mm256_add_pd(_sum0, _prob0); \
			_sum1 = _mm256_add_pd(_sum1, _prob1); \
			_mm256_store_pd(&prob[i], _prob0); \
			_mm256_store_pd(&prob[i+4], _prob1); \
		} \
		\
		_mm256_store_pd(&probSumH[0], _sum0); \
		_mm256_store_pd(&probSumH[4], _sum1); \
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7]; \
	} else { \
		COLLAPSE_MIS(); \
	} \
} while(0)

#endif // _SUPER_SITE_MACROS_H

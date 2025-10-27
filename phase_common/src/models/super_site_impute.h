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

#ifndef _SUPER_SITE_IMPUTE_H
#define _SUPER_SITE_IMPUTE_H

#include <models/super_site_accessor.h>

// ============================================================================
// IMPUTE_SUPERSITE - Impute missing genotype at super-site
// ============================================================================

#define IMPUTE_SUPERSITE(...) \
do { \
	const SuperSite& ss = super_sites[curr_abs_locus]; \
	\
	/* Unpack conditioning codes */ \
	aligned_vector32<uint8_t> cond_codes(n_cond_haps); \
	unpackSuperSiteCodesBatch( \
		G->SuperSitePanelCodes, \
		ss.panel_offset, \
		0, \
		n_cond_haps, \
		cond_codes.data()); \
	\
	/* For each conditioning haplotype, accumulate posterior */ \
	__m256d _alphaSum0 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][0]); \
	__m256d _alphaSum1 = _mm256_load_pd(&AlphaSumMissing[curr_rel_missing][4]); \
	__m256d _ones = _mm256_set1_pd(1.0f); \
	_alphaSum0 = _mm256_div_pd(_ones, _alphaSum0); \
	_alphaSum1 = _mm256_div_pd(_ones, _alphaSum1); \
	\
	__m256d _sumA0[2] = {_mm256_set1_pd(0.0f), _mm256_set1_pd(0.0f)}; \
	__m256d _sumA1[2] = {_mm256_set1_pd(0.0f), _mm256_set1_pd(0.0f)}; \
	\
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) { \
		uint8_t code = cond_codes[k]; \
		int allele = (code > 0) ? 1 : 0;  /* 0=REF, 1=ALT */ \
		\
		__m256d _prob0 = _mm256_load_pd(&prob[i]); \
		__m256d _prob1 = _mm256_load_pd(&prob[i+4]); \
		__m256d _alpha0 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i]); \
		__m256d _alpha1 = _mm256_load_pd(&AlphaMissing[curr_rel_missing][i+4]); \
		\
		__m256d _prod0 = _mm256_mul_pd(_mm256_mul_pd(_alpha0, _alphaSum0), _prob0); \
		__m256d _prod1 = _mm256_mul_pd(_mm256_mul_pd(_alpha1, _alphaSum1), _prob1); \
		\
		_sumA0[allele] = _mm256_add_pd(_sumA0[allele], _prod0); \
		_sumA1[allele] = _mm256_add_pd(_sumA1[allele], _prod1); \
	} \
	\
	/* Store imputation probabilities */ \
	double* p0 = (double*)&_sumA0[0]; \
	double* p1 = (double*)&_sumA0[1]; \
	for (int h = 0; h < 4; h++) { \
		missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = \
			p1[h] / (p0[h] + p1[h]); \
	} \
	p0 = (double*)&_sumA1[0]; \
	p1 = (double*)&_sumA1[1]; \
	for (int h = 4; h < HAP_NUMBER; h++) { \
		missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = \
			p1[h] / (p0[h] + p1[h]); \
	} \
} while(0)

#endif // _SUPER_SITE_IMPUTE_H

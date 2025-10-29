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

#ifndef _SUPER_SITE_EMISSIONS_H
#define _SUPER_SITE_EMISSIONS_H

#include <immintrin.h>
#include <cstdint>
#include <algorithm>
#include <models/super_site_accessor.h>

// ============================================================================
// AVX2 Vectorized Emission Precomputation
// ============================================================================

// Precompute emissions for all conditioning haplotypes using AVX2
// This is the critical performance function
inline void precomputeSuperSiteEmissions_AVX2(
	const uint8_t* cond_codes,    // Array of conditioning haplotype codes
	uint32_t n_cond_haps,
	uint8_t sample_code,
	double match_prob,            // Typically 1.0
	double mismatch_prob,         // Typically M.ed/M.ee
	aligned_vector32<double>& emissions_out) {

    // Process in batches of 8 for AVX2 (4 doubles per 256-bit register)
    for (uint32_t k = 0; k < n_cond_haps; k += 8) {
        uint32_t batch_size = std::min(8u, n_cond_haps - k);
        
        // Load up to 8 conditioning codes
        uint8_t codes_batch[16] = {0};
        for (uint32_t j = 0; j < batch_size; j++) {
            codes_batch[j] = cond_codes[k + j];
        }
        
        // Load 8 bytes into 128-bit register (will expand to 256-bit below)
        __m128i _codes_128 = _mm_loadu_si128((__m128i*)codes_batch);
        
        // Expand uint8 to uint32 (so we can use cmpeq_epi32)
        // Lower 8 codes: from bytes [0..7] of _codes_128
        __m256i _codes_lo = _mm256_cvtepu8_epi32(_codes_128);                // lanes 0..7 from bytes 0..7
        __m256i _codes_hi = _mm256_cvtepu8_epi32(_mm_srli_si128(_codes_128, 4)); // lanes 0..7 from bytes 4..11 (bytes 8..11 are zero)
        
        // Broadcast sample code to 8 lanes (as uint32)
        __m256i _sample_code_vec = _mm256_set1_epi32((uint32_t)sample_code);
        
        // Compare: produces 0xFFFFFFFF for match, 0x00000000 for mismatch
        __m256i _cmp_lo = _mm256_cmpeq_epi32(_codes_lo, _sample_code_vec);
        __m256i _cmp_hi = _mm256_cmpeq_epi32(_codes_hi, _sample_code_vec);
        
        // Create per-lane double masks using widening to epi64 then cast to pd
        // Take low 4 lanes from each compare result
        __m128i _cmp_lo_128 = _mm256_castsi256_si128(_cmp_lo);             // lanes 0..3
        __m128i _cmp_hi_128 = _mm256_castsi256_si128(_cmp_hi);             // lanes 0..3 correspond to codes 4..7
        __m256i _mask_lo_i64 = _mm256_cvtepi32_epi64(_cmp_lo_128);         // 4 x 64-bit: all-ones or zero
        __m256i _mask_hi_i64 = _mm256_cvtepi32_epi64(_cmp_hi_128);
        __m256d _mask_lo = _mm256_castsi256_pd(_mask_lo_i64);
        __m256d _mask_hi = _mm256_castsi256_pd(_mask_hi_i64);
        
        // Create probability vectors
        __m256d _match_vec = _mm256_set1_pd(match_prob);
        __m256d _mismatch_vec = _mm256_set1_pd(mismatch_prob);
        
        // Blend using bitwise masks: (mask & match) | (~mask & mismatch)
        __m256d _emit_lo = _mm256_or_pd(_mm256_and_pd(_mask_lo, _match_vec), _mm256_andnot_pd(_mask_lo, _mismatch_vec));
        __m256d _emit_hi = _mm256_or_pd(_mm256_and_pd(_mask_hi, _match_vec), _mm256_andnot_pd(_mask_hi, _mismatch_vec));

        // Store emissions without overrunning the buffer
        alignas(32) double emit_lo_buf[4];
        alignas(32) double emit_hi_buf[4];
        _mm256_storeu_pd(emit_lo_buf, _emit_lo);
        for (uint32_t j = 0; j < std::min<uint32_t>(batch_size, 4u); ++j) {
            emissions_out[k + j] = emit_lo_buf[j];
        }
        if (batch_size > 4) {
            _mm256_storeu_pd(emit_hi_buf, _emit_hi);
            uint32_t hi_count = batch_size - 4;
            for (uint32_t j = 0; j < hi_count; ++j) {
                emissions_out[k + 4 + j] = emit_hi_buf[j];
            }
        }
	}
}

// Fallback scalar version (for when AVX2 isn't worth the overhead)
inline void precomputeSuperSiteEmissions_Scalar(
    const uint8_t* cond_codes,
    uint32_t n_cond_haps,
    uint8_t sample_code,
    double match_prob,
    double mismatch_prob,
    aligned_vector32<double>& emissions_out) {

	for (uint32_t k = 0; k < n_cond_haps; k++) {
		emissions_out[k] = (cond_codes[k] == sample_code) ? match_prob : mismatch_prob;
	}
}

// Float version (scalar) for single-precision path
inline void precomputeSuperSiteEmissions_FloatScalar(
    const uint8_t* cond_codes,
    uint32_t n_cond_haps,
    uint8_t sample_code,
    float match_prob,
    float mismatch_prob,
    aligned_vector32<float>& emissions_out) {

    for (uint32_t k = 0; k < n_cond_haps; k++) {
        emissions_out[k] = (cond_codes[k] == sample_code) ? match_prob : mismatch_prob;
    }
}

#endif // _SUPER_SITE_EMISSIONS_H

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

#ifndef _SUPER_SITE_ACCESSOR_H
#define _SUPER_SITE_ACCESSOR_H

#include <cstdint>
#include <vector>
#include <boost/align/aligned_allocator.hpp>

template <typename T>
using aligned_vector32 = std::vector<T, boost::alignment::aligned_allocator<T, 32>>;

// ============================================================================
// Constants
// ============================================================================

#define SUPERSITE_MAX_ALTS 15
#define SUPERSITE_CODE_REF 0
#define SUPERSITE_CODE_MISSING 255

// ============================================================================
// Data Structures
// ============================================================================

struct SuperSite {
	uint32_t global_site_id;   // Index in overall variant array
	uint32_t chr;
	uint32_t bp;
	uint8_t n_alts;            // Number of ALTs (1-15, since code 0 is REF)
	uint32_t panel_offset;     // Byte offset in packed_allele_codes buffer
	// Note: bitwidth is always 4 bits (fixed), so codes are 4 bits each
};

// ============================================================================
// Unpacking Functions (Scalar)
// ============================================================================

// Unpack a single allele code for a conditioning haplotype at a super-site
// With fixed 4-bit codes: 2 codes per byte
inline uint8_t unpackSuperSiteCode(
	const uint8_t* packed_buffer,
	uint32_t panel_offset,
	uint32_t hap_idx) {

	// Byte index: 2 haplotypes per byte with 4-bit codes
	uint32_t byte_idx = panel_offset + hap_idx / 2;
	uint8_t byte_val = packed_buffer[byte_idx];

	// Bit shift: even haplotypes use lower 4 bits, odd use upper 4 bits
	uint8_t shift = (hap_idx % 2) * 4;
	return (byte_val >> shift) & 0x0F;
}

// Unpack multiple codes into a local buffer (for vectorization)
// Process in batches of 8, 16, or 32 for cache efficiency
inline void unpackSuperSiteCodesBatch(
	const uint8_t* packed_buffer,
	uint32_t panel_offset,
	uint32_t start_hap,
	uint32_t n_haps,
	uint8_t* codes_out) {

	for (uint32_t k = 0; k < n_haps; k++) {
		codes_out[k] = unpackSuperSiteCode(packed_buffer, panel_offset, start_hap + k);
	}
}

// Vectorized unpacking using PEXT (BMI2, Intel Haswell+)
// For maximum performance: requires -mbmi2 compiler flag
inline void unpackSuperSiteCodesVectorized_PEXT(
	const uint8_t* packed_buffer,
	uint32_t panel_offset,
	uint32_t n_haps,
	uint8_t* codes_out) {

	// This is complex and system-dependent; for now, use scalar
	// Can be optimized later if profiling shows it's a bottleneck
	unpackSuperSiteCodesBatch(packed_buffer, panel_offset, 0, n_haps, codes_out);
}

#endif

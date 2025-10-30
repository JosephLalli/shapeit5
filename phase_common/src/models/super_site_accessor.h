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
#include <utils/otools.h>
#include <objects/genotype/genotype_header.h>

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
	// Member variants span in super_site_var_index (start offset and count)
	uint32_t var_start;        // Start index into flattened super_site_var_index
	uint16_t var_count;        // Number of constituent variant indices
	// Note: bitwidth is always 4 bits (fixed), so codes are 4 bits each
	
	// Phase 3: Multi-class posterior indexing (per-sample, set during window setup)
	uint32_t class_prob_offset; // Offset into CurrentSuperClassPosteriors (HAP_NUMBER * C floats)
	uint8_t n_classes;          // C = 1 + n_alts (REF + ALT1..ALTn), cached for convenience
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

// ============================================================================
// Sample code accessor
// ============================================================================

// Returns 0 (REF) if no ALTs carried across member variants, ALT index+1 for first ALT seen,
// or SUPERSITE_CODE_MISSING if all member records are missing for this sample at this locus.
inline uint8_t getSampleSuperSiteAlleleCode(
    const genotype* G,
    const SuperSite& ss,
    const std::vector<int>& super_site_var_index,
    int hap // 0 or 1
) {
    bool any_non_missing = false;
    for (uint32_t i = 0; i < ss.var_count; ++i) {
        int v_idx = super_site_var_index[ss.var_start + i];
        unsigned char v = G->Variants[DIV2(v_idx)];
        if (VAR_GET_MIS(MOD2(v_idx), v)) continue; // missing at this component
        any_non_missing = true;
        bool carries = (hap == 0) ? VAR_GET_HAP0(MOD2(v_idx), v) : VAR_GET_HAP1(MOD2(v_idx), v);
        if (carries) return static_cast<uint8_t>(i + 1);
    }
    if (!any_non_missing) return SUPERSITE_CODE_MISSING;
    return SUPERSITE_CODE_REF;
}

// ============================================================================
// Supersite Classification
// ============================================================================

enum class SSClass { MIS, HOM, AMB };

// Classify a supersite genotype based on both haplotype codes
inline SSClass classify_supersite(
    const genotype* G,
    const SuperSite& ss,
    const std::vector<int>& super_site_var_index,
    uint8_t& c0,
    uint8_t& c1
) {
    c0 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 0);
    c1 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 1);
    
    if (c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING) 
        return SSClass::MIS;
    if (c0 == c1) 
        return SSClass::HOM;
    return SSClass::AMB;
}

#endif

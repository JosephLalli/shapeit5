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
#include <cstdlib>
#include <cinttypes>
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
#define SUPERSITE_CODE_CONFLICT 254

// Enforce a deterministic ordering for class pairs used in emissions.
// Keeps the lowest non-missing code in c0 so biallelic parity does not
// depend on hap0/hap1 orientation.
inline void canonicalize_class_pair(uint8_t& c0, uint8_t& c1) {
	if (c0 != SUPERSITE_CODE_MISSING && c1 != SUPERSITE_CODE_MISSING && c0 > c1) {
		uint8_t tmp = c0;
		c0 = c1;
		c1 = tmp;
	}
}

// ============================================================================
// Data Structures
// ============================================================================

struct SuperSite {
	uint32_t global_site_id;   // Index in overall variant array
	uint32_t chr;
	uint32_t bp;
	uint8_t n_alts;            // Number of ALTs (1-15, since code 0 is REF)
	uint32_t panel_offset;     // Byte offset in packed_allele_codes buffer
	uint32_t panel_span_bytes; // Number of bytes reserved in packed buffer for this supersite
	// Member variants span in super_site_var_index (start offset and count)
	uint32_t var_start;        // Start index into flattened super_site_var_index
	uint16_t var_count;        // Number of constituent variant indices
	// Note: bitwidth is always 4 bits (fixed), so codes are 4 bits each
	
	// Phase 3: Multi-class posterior indexing (n_classes cached for convenience)
	// NOTE: class_prob_offset moved to thread-local storage to fix race condition
	uint8_t n_classes;          // C = 1 + n_alts (REF + ALT1..ALTn), cached for convenience

	// Rare code gating (mirrors biallelic rare_allele[])
	// Bit i = 1 means code i is rare (freq < RARE_VARIANT_FREQ in panel)
	// Used by RUN_HOM to skip emission updates for common alleles
	uint16_t rare_code_mask;
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

	// Handle null buffer case (no supersites built)
	if (packed_buffer == nullptr) {
		return 0; // Return REF code as safe fallback
	}

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

// Classify observed genotype from immutable c0/c1 values
// Input: c0, c1 from getSupersiteObservedGt() (can be 0-15 or MISSING=255)
// Output: SSClass (MIS, HOM, or AMB)
inline SSClass classifyObservedGt(uint8_t c0, uint8_t c1) {
	if (c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING)
		return SSClass::MIS;
	if (c0 == c1)
		return SSClass::HOM;
	return SSClass::AMB;
}

// Classify a supersite genotype based on both haplotype codes
// NOTE: This function reads from MUTABLE Variants[] and should only be used
// for dynamic classification. For emission-based dispatch, use classifyObservedGt()
// with values from getSupersiteObservedGt() instead.
inline SSClass classify_supersite(
    const genotype* G,
    const SuperSite& ss,
    const std::vector<int>& super_site_var_index,
    uint8_t& c0,
    uint8_t& c1
) {
    c0 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 0);
    c1 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 1);
    canonicalize_class_pair(c0, c1);
    
    if (c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING) 
        return SSClass::MIS;
    if (c0 == c1) 
        return SSClass::HOM;
    return SSClass::AMB;
}

// Check if a supersite is heterozygous (different alleles on each haplotype)
inline bool isSuperSiteHeterozygous(
    const genotype* G,
    const SuperSite& ss,
    const std::vector<int>& super_site_var_index
) {
    uint8_t c0 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 0);
    uint8_t c1 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, 1);
    return c0 != c1;
}

#endif

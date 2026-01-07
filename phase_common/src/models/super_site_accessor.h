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
#include <array>
#include <vector>
#include <boost/align/aligned_allocator.hpp>
#include <utils/otools.h>
#include <objects/genotype/genotype_header.h>

// ============================================================================
// Constants
// ============================================================================

#define SUPERSITE_MAX_ALTS 255
#define SUPERSITE_CODE_REF 0
// Diagnostic-only sentinel codes used by test helpers (not serialized).
// These may overlap with valid allele codes if n_alts is very large.
#define SUPERSITE_CODE_CONFLICT 254
#define SUPERSITE_CODE_MISSING 255

// Enforce a deterministic ordering for class pairs used in emissions.
// Keeps the lowest code in c0 so emissions do not depend on hap0/hap1 orientation.
inline void canonicalize_class_pair(uint8_t& c0, uint8_t& c1) {
	if (c0 > c1) {
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
	uint16_t n_alts;           // Number of ALTs (1-255, since code 0 is REF)
	uint32_t panel_offset;     // Byte offset in packed_allele_codes buffer
	uint32_t panel_span_bytes; // Number of bytes reserved in packed buffer for this supersite
	// Member variants span in super_site_var_index (start offset and count)
	uint32_t var_start;        // Start index into flattened super_site_var_index
	uint16_t var_count;        // Number of constituent variant indices
	// Note: codes are 1 byte each (0=REF, 1..n_alts=ALT1..ALTn)
	
	// Phase 3: Multi-class posterior indexing (n_classes cached for convenience)
	// NOTE: class_prob_offset moved to thread-local storage to fix race condition
	uint16_t n_classes;          // C = 1 + n_alts (REF + ALT1..ALTn), cached for convenience

	// Rare code gating (mirrors biallelic rare_allele[])
	// Bit i = 1 means code i is rare (freq < RARE_VARIANT_FREQ in panel)
	// Used by RUN_HOM to skip emission updates for common alleles
	std::array<uint64_t, 4> rare_code_mask;
};

inline bool supersite_code_is_rare(const SuperSite& ss, uint8_t code) {
	if (code > ss.n_alts) return false;
	const uint16_t idx = static_cast<uint16_t>(code) >> 6;
	const uint16_t bit = static_cast<uint16_t>(code) & 63u;
	return (ss.rare_code_mask[idx] & (1ULL << bit)) != 0u;
}

inline bool supersite_has_rare_mask(const SuperSite& ss) {
	return (ss.rare_code_mask[0] | ss.rare_code_mask[1] |
	        ss.rare_code_mask[2] | ss.rare_code_mask[3]) != 0u;
}

// ============================================================================
// Unpacking Functions (Scalar)
// ============================================================================

// Unpack a single allele code for a conditioning haplotype at a supersite
// One byte per haplotype
inline uint8_t unpackSuperSiteCode(
	const uint8_t* packed_buffer,
	uint32_t panel_offset,
	uint32_t hap_idx) {

	// Handle null buffer case (no supersites built)
	if (packed_buffer == nullptr) {
		return 0; // Return REF code as safe fallback
	}

	// One byte per haplotype
	return packed_buffer[panel_offset + hap_idx];
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

// Infer the per-haplotype allele code for a supersite in a sample.
// Prefers stored supersite codes when available; falls back to split-variant bits.
inline uint8_t getSampleSuperSiteAlleleCode(
	const genotype* g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	int hap) {
	if (!g || (hap != 0 && hap != 1)) return SUPERSITE_CODE_MISSING;

	int ss_idx = -1;
	if (g->locus_to_super_idx &&
	    ss.global_site_id < g->locus_to_super_idx->size()) {
		ss_idx = (*g->locus_to_super_idx)[ss.global_site_id];
	}
	if (ss_idx >= 0) {
		if (g->supersiteIsMissing(ss_idx)) return SUPERSITE_CODE_MISSING;
		const size_t offset = static_cast<size_t>(ss_idx) * 2u;
		if (g->ss_phased_gts.size() >= offset + 2u) {
			return (hap == 0) ? g->ss_phased_gts[offset]
			                  : g->ss_phased_gts[offset + 1];
		}
		if (g->ss_observed_gts.size() >= offset + 2u) {
			return (hap == 0) ? g->ss_observed_gts[offset]
			                  : g->ss_observed_gts[offset + 1];
		}
	}

	if (super_site_var_index.empty() || ss.var_count == 0) {
		return SUPERSITE_CODE_REF;
	}
	bool missing = false;
	int alt_index = -1;
	for (uint32_t i = 0; i < ss.var_count; ++i) {
		const size_t offset = ss.var_start + i;
		if (offset >= super_site_var_index.size()) break;
		const int v_idx = super_site_var_index[offset];
		if (v_idx < 0) continue;
		const unsigned char byte = g->Variants[DIV2(v_idx)];
		if (VAR_GET_MIS(MOD2(v_idx), byte)) {
			missing = true;
			break;
		}
		const bool carries = (hap == 0)
			? VAR_GET_HAP0(MOD2(v_idx), byte)
			: VAR_GET_HAP1(MOD2(v_idx), byte);
		if (carries) {
			if (alt_index >= 0) return SUPERSITE_CODE_CONFLICT;
			alt_index = static_cast<int>(i);
		}
	}
	if (missing) return SUPERSITE_CODE_MISSING;
	if (alt_index < 0) return SUPERSITE_CODE_REF;
	return static_cast<uint8_t>(alt_index + 1);
}

inline bool isSuperSiteHeterozygous(
	const genotype* g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index) {
	const uint8_t c0 = getSampleSuperSiteAlleleCode(g, ss, super_site_var_index, 0);
	const uint8_t c1 = getSampleSuperSiteAlleleCode(g, ss, super_site_var_index, 1);
	if (c0 == SUPERSITE_CODE_MISSING || c1 == SUPERSITE_CODE_MISSING) return false;
	if (c0 == SUPERSITE_CODE_CONFLICT || c1 == SUPERSITE_CODE_CONFLICT) return true;
	return c0 != c1;
}

// ============================================================================
// Supersite Classification
// ============================================================================

enum class SSClass { MIS, HOM, AMB };

// Classify observed genotype from immutable c0/c1 values and missing mask.
// Input: c0, c1 from getSupersiteObservedGt(); is_missing from supersiteIsMissing().
// Output: SSClass (MIS, HOM, or AMB)
inline SSClass classifyObservedGt(uint8_t c0, uint8_t c1, bool is_missing) {
	if (is_missing) return SSClass::MIS;
	if (c0 == c1) return SSClass::HOM;
	return SSClass::AMB;
}

inline SSClass classify_supersite(
	const genotype* g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	uint8_t& c0_out,
	uint8_t& c1_out) {
	c0_out = getSampleSuperSiteAlleleCode(g, ss, super_site_var_index, 0);
	c1_out = getSampleSuperSiteAlleleCode(g, ss, super_site_var_index, 1);
	if (c0_out == SUPERSITE_CODE_MISSING || c1_out == SUPERSITE_CODE_MISSING) {
		return SSClass::MIS;
	}
	if (c0_out == SUPERSITE_CODE_CONFLICT || c1_out == SUPERSITE_CODE_CONFLICT) {
		return SSClass::AMB;
	}
	canonicalize_class_pair(c0_out, c1_out);
	return (c0_out == c1_out) ? SSClass::HOM : SSClass::AMB;
}

#endif

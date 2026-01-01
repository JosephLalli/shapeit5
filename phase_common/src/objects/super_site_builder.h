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

#ifndef _SUPER_SITE_BUILDER_H
#define _SUPER_SITE_BUILDER_H

#include <vector>
#include <cstdint>
#include <models/super_site_accessor.h>

// Forward declarations (these types will be defined in the actual codebase)
class variant_map;
class conditioning_set;
class genotype_set;
class genotype;

// ============================================================================
// Supersite Construction (1:1 with multiallelic variants)
// ============================================================================

void buildSuperSites(
    variant_map& V,
    conditioning_set& H,
    std::vector<SuperSite>& super_sites_out,
    std::vector<bool>& is_super_site_out,
    std::vector<uint8_t>& packed_allele_codes_out,
    std::vector<int>& locus_to_super_idx_out,
    std::vector<int>& super_site_var_index_out);

// Build a lookup table mapping each locus to its supersite anchor (or -1 if none).
std::vector<int> buildSupersiteAnchorMap(const std::vector<SuperSite>& super_sites,
                                         const std::vector<int>& super_site_var_index,
                                         size_t n_loci);

// Update anchor variant encoding to reflect supersite genotype status
// Must be called after setSuperSiteContext() and before genotype::build()
void updateSuperSiteAnchorEncoding(genotype_set& G,
                                   const std::vector<SuperSite>& super_sites,
                                   const std::vector<int>& super_site_var_index);

// Resolve per-sample supersite allele classes from stored observed codes.
void resolveSupersiteClasses(
	genotype& g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	uint8_t& c0_out,
	uint8_t& c1_out);

#endif // _SUPER_SITE_BUILDER_H

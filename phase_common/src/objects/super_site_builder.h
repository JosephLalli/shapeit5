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

// ============================================================================
// Super-Site Detection and Construction
// ============================================================================

void buildSuperSites(
	variant_map& V,
	conditioning_set& H,
	std::vector<SuperSite>& super_sites_out,
	std::vector<bool>& is_super_site_out,
	std::vector<uint8_t>& packed_allele_codes_out,
	std::vector<uint8_t>& sample_supersite_genotypes_out);

#endif // _SUPER_SITE_BUILDER_H

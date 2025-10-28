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

#include <objects/super_site_builder.h>
#include <models/super_site_accessor.h>
#include <map>
#include <algorithm>
#include <utility>

// Note: This implementation is a skeleton that will need to be integrated
// with the actual variant_map and conditioning_set classes when they are available

/*
// ============================================================================
// Super-Site Detection and Construction
// ============================================================================

void buildSuperSites(
	variant_map& V,
	conditioning_set& H,
	std::vector<SuperSite>& super_sites_out,
	std::vector<bool>& is_super_site_out,
	std::vector<uint8_t>& packed_allele_codes_out,
	std::vector<uint8_t>& sample_supersite_genotypes_out) {

	super_sites_out.clear();
	is_super_site_out.assign(V.size(), false);
	packed_allele_codes_out.clear();

	// Group variants by (chr, bp)
	std::map<std::pair<int,int>, std::vector<int>> sites_by_pos;
	for (int v = 0; v < V.size(); v++) {
		auto key = std::make_pair(V.vec_pos[v]->chr, V.vec_pos[v]->bp);
		sites_by_pos[key].push_back(v);
	}

	// Process each group
	uint32_t current_panel_offset = 0;

	for (auto& [pos_key, variant_indices] : sites_by_pos) {
		if (variant_indices.size() == 1) {
			// Normal biallelic site, skip
			continue;
		}
		
		// Multiple variants at same position → super-site
		uint8_t n_alts = std::min((uint8_t)variant_indices.size(), (uint8_t)SUPERSITE_MAX_ALTS);
		
		if (n_alts > SUPERSITE_MAX_ALTS) {
			// TODO: Handle >16 ALTs by splitting into multiple super-sites
			// vrb.warning("Super-site at " + std::to_string(pos_key.first) + ":" + std::to_string(pos_key.second) + 
			//            " has " + std::to_string(n_alts) + " ALTs; truncating to " + std::to_string(SUPERSITE_MAX_ALTS));
			n_alts = SUPERSITE_MAX_ALTS;
		}
		
		SuperSite ss;
		ss.global_site_id = variant_indices[0];
		ss.chr = pos_key.first;
		ss.bp = pos_key.second;
		ss.n_alts = n_alts;
		ss.panel_offset = current_panel_offset;
		
		// For each reference haplotype, determine allele code
		std::vector<uint8_t> haplotype_codes(H.n_haplotypes);
		
		for (int h = 0; h < H.n_haplotypes; h++) {
			uint8_t code = SUPERSITE_CODE_REF;  // Default to REF
			int n_alts_seen = 0;
			
			for (int alt_idx = 0; alt_idx < n_alts; alt_idx++) {
				int v = variant_indices[alt_idx];
				
				// Check if haplotype h carries ALT at variant v
				// Using existing bitmatrix accessor
				if (H.Hhap.get(v, h)) {
					if (code == SUPERSITE_CODE_REF) {
						// First ALT this haplotype carries
						code = alt_idx + 1;  // Code 1, 2, 3, ...
						n_alts_seen = 1;
					} else {
						// Haplotype carries multiple ALTs → violation
						// Canonicalization: keep first ALT, or apply deterministic policy
						// For STRs, pick longer allele (TODO: implement)
						n_alts_seen++;
					}
				}
			}
			
			if (n_alts_seen > 1) {
				// Warn but don't fail; use first ALT
				// vrb.warning("Reference haplotype " + std::to_string(h) + " carries " + std::to_string(n_alts_seen) +
				//            " ALTs at position " + std::to_string(ss.bp) + "; using first");
			}
			
			haplotype_codes[h] = code;
		}
		
		// Pack codes: 2 codes per byte (4 bits each)
		uint32_t n_bytes = (H.n_haplotypes + 1) / 2;
		for (uint32_t byte_idx = 0; byte_idx < n_bytes; byte_idx++) {
			uint8_t packed_byte = 0;
			for (int code_idx = 0; code_idx < 2; code_idx++) {
				int hap_idx = byte_idx * 2 + code_idx;
				if (hap_idx < H.n_haplotypes) {
					uint8_t code = haplotype_codes[hap_idx];
					packed_byte |= (code << (code_idx * 4));
				}
			}
			packed_allele_codes_out.push_back(packed_byte);
		}
		
		current_panel_offset += n_bytes;
		
		// Mark all variants in this group as part of super-site
		super_sites_out.push_back(ss);
		for (int v : variant_indices) {
			is_super_site_out[v] = true;
		}
	}

	// vrb.bullet("Built " + std::to_string(super_sites_out.size()) + " super-sites");
}
*/

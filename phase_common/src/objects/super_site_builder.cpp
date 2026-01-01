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
#include <containers/variant_map.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/genotype_set.h>
#include <objects/genotype/genotype_header.h>
#include <cassert>
#include <stdexcept>

namespace {
inline void rare_mask_set(std::array<uint64_t, 4>& mask, uint16_t code) {
	const uint16_t idx = code >> 6;
	const uint16_t bit = code & 63u;
	mask[idx] |= (1ULL << bit);
}
}

// Build supersites directly from multiallelic records.
// One supersite per multiallelic variant; packed codes are 1 byte per haplotype.
void buildSuperSites(
    variant_map& V,
    conditioning_set& H,
    std::vector<SuperSite>& super_sites_out,
    std::vector<bool>& is_super_site_out,
    std::vector<uint8_t>& packed_allele_codes_out,
    std::vector<int>& locus_to_super_idx_out,
    std::vector<int>& super_site_var_index_out)
{
    super_sites_out.clear();
    packed_allele_codes_out.clear();
    super_site_var_index_out.clear();
    is_super_site_out.assign(V.size(), false);
    locus_to_super_idx_out.assign(V.size(), -1);

    const size_t expected_codes = static_cast<size_t>(H.n_supersites) * H.n_hap;
    if (H.n_supersites > 0 && H.H_supersite_codes.size() < expected_codes) {
        throw std::runtime_error("Supersite codes missing: expected " +
                                 std::to_string(expected_codes) + " codes, got " +
                                 std::to_string(H.H_supersite_codes.size()));
    }

    uint32_t current_panel_offset = 0;
    size_t ss_idx = 0;

    for (int v = 0; v < V.size(); v++) {
        variant* vp = V.vec_pos[v];
        if (vp->n_alts <= 1) continue;
        if (vp->n_alts > SUPERSITE_MAX_ALTS) {
            throw std::runtime_error("Supersite n_alts exceeds SUPERSITE_MAX_ALTS at " +
                                     vp->chr + ":" + std::to_string(vp->bp));
        }

        SuperSite ss;
        ss.global_site_id = static_cast<uint32_t>(v);
        ss.chr = 0; // chr is not used downstream; keep 0 to avoid string->int issues
        ss.bp = static_cast<uint32_t>(vp->bp);
        ss.n_alts = vp->n_alts;
        ss.n_classes = static_cast<uint16_t>(vp->n_alts + 1u);
        ss.panel_offset = current_panel_offset;
        ss.panel_span_bytes = static_cast<uint32_t>(H.n_hap);
        ss.var_start = static_cast<uint32_t>(super_site_var_index_out.size());
        ss.var_count = 1;
        ss.rare_code_mask = {0u, 0u, 0u, 0u};

        super_site_var_index_out.push_back(v);
        is_super_site_out[v] = true;
        locus_to_super_idx_out[v] = static_cast<int>(super_sites_out.size());

        const size_t base = ss_idx * static_cast<size_t>(H.n_hap);
        if (base + H.n_hap > H.H_supersite_codes.size()) {
            throw std::runtime_error("Supersite code offset out of range for supersite index " +
                                     std::to_string(ss_idx));
        }

        std::vector<uint32_t> code_counts(static_cast<size_t>(ss.n_classes), 0u);
        for (unsigned long h = 0; h < H.n_hap; ++h) {
            uint8_t code = H.H_supersite_codes[base + h];
            if (code >= ss.n_classes) code = 0;
            code_counts[code]++;
        }

        for (uint16_t c = 0; c < ss.n_classes; ++c) {
            float freq = static_cast<float>(code_counts[c]) / static_cast<float>(H.n_hap);
            if (freq < RARE_VARIANT_FREQ) {
                rare_mask_set(ss.rare_code_mask, c);
            }
        }

        for (unsigned long h = 0; h < H.n_hap; ++h) {
            uint8_t code = H.H_supersite_codes[base + h];
            if (code >= ss.n_classes) code = 0;
            packed_allele_codes_out.push_back(code);
        }

        current_panel_offset += static_cast<uint32_t>(H.n_hap);
        super_sites_out.push_back(ss);
        ++ss_idx;
    }
}

// Update anchor variant encoding to reflect supersite genotype status
// Must be called after setSuperSiteContext() and before genotype::build()
void updateSuperSiteAnchorEncoding(genotype_set& G,
                                   const std::vector<SuperSite>& super_sites,
                                   const std::vector<int>& super_site_var_index)
{
	(void)super_site_var_index;
	// Update anchor encoding for all samples
    for (unsigned int i = 0; i < G.n_ind; i++) {
        genotype* g = G.vecG[i];
        
        // Supersite context must be set before calling this function
        assert(g->super_sites && g->locus_to_super_idx && g->super_site_var_index);
        
        // Update anchor encoding based on observed supersite genotype
        for (size_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
            const SuperSite& ss = super_sites[ss_idx];
            unsigned char& anchor_byte = g->Variants[DIV2(ss.global_site_id)];
            if (g->supersiteIsMissing(static_cast<int>(ss_idx))) {
                VAR_SET_MIS(MOD2(ss.global_site_id), anchor_byte);
                VAR_CLR_HAP0(MOD2(ss.global_site_id), anchor_byte);
                VAR_CLR_HAP1(MOD2(ss.global_site_id), anchor_byte);
                continue;
            }
			if (VAR_GET_SCA(MOD2(ss.global_site_id), anchor_byte)) {
				continue;
			}
            uint8_t c0 = 0u;
            uint8_t c1 = 0u;
            g->getSupersiteObservedGt(static_cast<int>(ss_idx), c0, c1);
            if (c0 == c1) {
                VAR_SET_HOM(MOD2(ss.global_site_id), anchor_byte);
            } else {
                VAR_SET_HET(MOD2(ss.global_site_id), anchor_byte);
            }
		}
	}
}

void resolveSupersiteClasses(
	genotype& g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	uint8_t& c0_out,
	uint8_t& c1_out)
{
	(void)super_site_var_index;
	const int locus = static_cast<int>(ss.global_site_id);
	int ss_idx = -1;
	if (g.locus_to_super_idx && locus >= 0 &&
	    locus < static_cast<int>(g.locus_to_super_idx->size())) {
		ss_idx = (*g.locus_to_super_idx)[locus];
	}
	c0_out = 0u;
	c1_out = 0u;
	if (ss_idx >= 0) {
		g.getSupersiteObservedGt(ss_idx, c0_out, c1_out);
	}
}

std::vector<int> buildSupersiteAnchorMap(const std::vector<SuperSite>& super_sites,
                                         const std::vector<int>& super_site_var_index,
                                         size_t n_loci) {
	std::vector<int> redirect(n_loci, -1);
	for (const auto& ss : super_sites) {
		const int anchor = static_cast<int>(ss.global_site_id);
		if (anchor >= 0 && static_cast<size_t>(anchor) < n_loci) redirect[anchor] = anchor;
		for (uint16_t ai = 0; ai < ss.var_count; ++ai) {
			const size_t offset = ss.var_start + ai;
			if (offset >= super_site_var_index.size()) continue;
			const int member = super_site_var_index[offset];
			if (member >= 0 && static_cast<size_t>(member) < n_loci) {
				redirect[member] = anchor;
			}
		}
	}
	return redirect;
}

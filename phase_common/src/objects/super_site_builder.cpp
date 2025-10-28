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
#include <map>
#include <algorithm>
#include <utility>

// Build super-sites by collapsing split biallelic records at identical (chr,bp)
// positions into multi-allelic "super-sites" and packing per-haplotype 4-bit codes.
void buildSuperSites(
    variant_map& V,
    conditioning_set& H,
    std::vector<SuperSite>& super_sites_out,
    std::vector<bool>& is_super_site_out,
    std::vector<uint8_t>& packed_allele_codes_out,
    std::vector<int>& locus_to_super_idx_out,
    std::vector<int>& super_site_var_index_out,
    std::vector<uint8_t>& /*sample_supersite_genotypes_out*/)
{
    super_sites_out.clear();
    packed_allele_codes_out.clear();
    super_site_var_index_out.clear();
    is_super_site_out.assign(V.size(), false);
    locus_to_super_idx_out.assign(V.size(), -1);

    // Group variants by (chr,bp)
    using Key = std::pair<std::string,int>;
    std::map<Key, std::vector<int>> sites_by_pos;
    for (int v = 0; v < V.size(); v++) {
        variant* vp = V.vec_pos[v];
        sites_by_pos[{vp->chr, vp->bp}].push_back(v);
    }

    uint32_t current_panel_offset = 0;

    for (auto& kv : sites_by_pos) {
        const std::vector<int>& variant_indices = kv.second;
        if (variant_indices.size() <= 1) continue; // not multi-allelic

        // Split groups larger than SUPERSITE_MAX_ALTS into deterministic chunks
        for (size_t chunk_start = 0; chunk_start < variant_indices.size(); chunk_start += SUPERSITE_MAX_ALTS) {
            uint8_t n_alts = static_cast<uint8_t>(std::min<size_t>(SUPERSITE_MAX_ALTS, variant_indices.size() - chunk_start));

            SuperSite ss;
            ss.global_site_id = variant_indices[chunk_start];
            ss.chr = 0; // chr is not used downstream; keep 0 to avoid string->int issues
            ss.bp = static_cast<uint32_t>(V.vec_pos[variant_indices[chunk_start]]->bp);
            ss.n_alts = n_alts;
            ss.panel_offset = current_panel_offset;
            ss.var_start = static_cast<uint32_t>(super_site_var_index_out.size());
            ss.var_count = n_alts;

            // Record member variant indices and mark mappings
            for (uint8_t ai = 0; ai < n_alts; ++ai) {
                int v_idx = variant_indices[chunk_start + ai];
                super_site_var_index_out.push_back(v_idx);
                is_super_site_out[v_idx] = true;
                locus_to_super_idx_out[v_idx] = static_cast<int>(super_sites_out.size());
            }

            // Build per-haplotype codes (0=REF, 1..n_alts=first ALT seen)
            std::vector<uint8_t> hap_codes;
            hap_codes.resize(static_cast<size_t>(H.n_hap), 0);

            for (unsigned long h = 0; h < H.n_hap; ++h) {
                uint8_t code = SUPERSITE_CODE_REF;
                for (uint8_t ai = 0; ai < n_alts; ++ai) {
                    int v_idx = variant_indices[chunk_start + ai];
                    if (H.H_opt_var.get(static_cast<unsigned int>(v_idx), static_cast<unsigned int>(h))) {
                        code = static_cast<uint8_t>(ai + 1);
                        break; // take first ALT encountered
                    }
                }
                hap_codes[static_cast<size_t>(h)] = code;
            }

            // Pack codes: 2 per byte (4 bits each)
            uint32_t n_bytes = static_cast<uint32_t>((H.n_hap + 1) / 2);
            for (uint32_t byte_idx = 0; byte_idx < n_bytes; ++byte_idx) {
                uint8_t packed = 0;
                uint32_t hap0 = byte_idx * 2;
                uint32_t hap1 = hap0 + 1;
                if (hap0 < H.n_hap) packed |= static_cast<uint8_t>((hap_codes[hap0] & 0x0F) << 0);
                if (hap1 < H.n_hap) packed |= static_cast<uint8_t>((hap_codes[hap1] & 0x0F) << 4);
                packed_allele_codes_out.push_back(packed);
            }
            current_panel_offset += n_bytes;

            super_sites_out.push_back(ss);
        }
    }
}

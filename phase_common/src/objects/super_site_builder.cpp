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
#include <map>
#include <algorithm>
#include <utility>
#include <cassert>

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
    std::vector<uint8_t>& /*sample_supersite_genotypes_out*/,
    int mac_threshold)
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

        std::vector<int> kept_indices;
        kept_indices.reserve(variant_indices.size());
        for (int v_idx : variant_indices) {
            variant* vp = V.vec_pos[v_idx];
            if (mac_threshold > 0 && vp->getMAC() < static_cast<unsigned int>(mac_threshold)) {
                continue; // treat low-MAC sibling as biallelic
            }
            kept_indices.push_back(v_idx);
        }

        if (kept_indices.size() <= 1) continue; // nothing left to form a supersite

        // Split groups larger than SUPERSITE_MAX_ALTS into deterministic chunks
        for (size_t chunk_start = 0; chunk_start < kept_indices.size(); chunk_start += SUPERSITE_MAX_ALTS) {
            uint8_t n_alts = static_cast<uint8_t>(std::min<size_t>(SUPERSITE_MAX_ALTS, kept_indices.size() - chunk_start));
            if (n_alts <= 1) continue; // chunk degenerates to biallelic after filtering

            SuperSite ss;
            ss.global_site_id = kept_indices[chunk_start];
            ss.chr = 0; // chr is not used downstream; keep 0 to avoid string->int issues
            ss.bp = static_cast<uint32_t>(V.vec_pos[kept_indices[chunk_start]]->bp);
            ss.n_alts = n_alts;
            ss.panel_offset = current_panel_offset;
            ss.var_start = static_cast<uint32_t>(super_site_var_index_out.size());
            ss.var_count = n_alts;
            ss.n_classes = static_cast<uint8_t>(1 + n_alts);

            // Record member variant indices and mark mappings
            for (uint8_t ai = 0; ai < n_alts; ++ai) {
                int v_idx = kept_indices[chunk_start + ai];
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
                    int v_idx = kept_indices[chunk_start + ai];
                    if (H.H_opt_var.get(static_cast<unsigned int>(v_idx), static_cast<unsigned int>(h))) {
                        code = static_cast<uint8_t>(ai + 1);
                        break; // take first ALT encountered
                    }
                }
                hap_codes[static_cast<size_t>(h)] = code;
            }

            // Pack codes: 2 per byte (4 bits each)
            uint32_t n_bytes = static_cast<uint32_t>((H.n_hap + 1) / 2);
            ss.panel_span_bytes = n_bytes;
            
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

    // Diagnostics: summarize build when tracing is enabled
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    if (tr && tr[0] != '\0' && tr[0] != '0') {
        std::fprintf(stdout, "buildSuperSites: n_sites=%zu n_hap=%lu n_supersites=%zu packed_bytes=%zu\n",
                     (size_t)V.size(), H.n_hap, super_sites_out.size(), packed_allele_codes_out.size());
        for (size_t i = 0; i < super_sites_out.size(); ++i) {
            const SuperSite& s = super_sites_out[i];
            std::fprintf(stdout, "  ss[%zu]: anchor=%u var_count=%u panel_off=%u span_bytes=%u\n",
                         i, s.global_site_id, (unsigned)s.var_count, (unsigned)s.panel_offset, (unsigned)s.panel_span_bytes);
        }
    }
}

// Update anchor variant encoding to reflect supersite genotype status
// Must be called after setSuperSiteContext() and before genotype::build()
void updateSuperSiteAnchorEncoding(genotype_set& G,
                                   const std::vector<SuperSite>& super_sites,
                                   const std::vector<int>& super_site_var_index)
{
    // Update anchor encoding for all samples
    for (unsigned int i = 0; i < G.n_ind; i++) {
        genotype* g = G.vecG[i];
        
        // Supersite context must be set before calling this function
        assert(g->super_sites && g->locus_to_super_idx && g->super_site_var_index);
        
        // Update anchor encoding for heterozygous multiallelic sites only
        // (HOM and MIS anchors are already correct from VCF reading)
        for (size_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
            const SuperSite& ss = super_sites[ss_idx];
            
            // Only heterozygous multiallelic sites need correction
            // (e.g., ALT2|ALT3 appears as 0/0 at anchor but is actually HET)
            if (isSuperSiteHeterozygous(g, ss, super_site_var_index)) {
                VAR_SET_HET(MOD2(ss.global_site_id), g->Variants[DIV2(ss.global_site_id)]);
            }
        }
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

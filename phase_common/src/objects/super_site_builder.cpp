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
#include <map>
#include <algorithm>
#include <utility>
#include <cassert>
#include <random>
#include <stdexcept>
#include <unordered_set>
#include <sstream>

namespace {
// Parse comma-separated bp targets from SHAPEIT5_TRACE_BP
const std::unordered_set<int>& bp_trace_targets() {
	static std::unordered_set<int> targets;
	static bool initialized = false;
	if (!initialized) {
		initialized = true;
		const char* env = std::getenv("SHAPEIT5_TRACE_BP");
		if (env && env[0]) {
			std::stringstream ss(env);
			std::string tok;
			while (std::getline(ss, tok, ',')) {
				try {
					int bp = std::stoi(tok);
					targets.insert(bp);
				} catch (const std::exception&) {
					// ignore parse errors
				}
			}
		}
	}
	return targets;
}

bool should_trace_bp(int bp) {
	const auto& targets = bp_trace_targets();
	return !targets.empty() && targets.find(bp) != targets.end();
}
} // namespace

// Build super-sites by collapsing split biallelic records at identical (chr,bp)
// positions into multi-allelic "super-sites" and packing per-haplotype 4-bit codes.
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
        const int site_bp = kv.first.second;
        const bool trace_bp = should_trace_bp(site_bp);
        if (trace_bp) {
            std::fprintf(stdout,
                         "[BP_SUPERSITE_TRACE] bp=%d variants=%zu indices=",
                         site_bp,
                         variant_indices.size());
            for (size_t i = 0; i < variant_indices.size(); ++i) {
                if (i) std::fprintf(stdout, ",");
                std::fprintf(stdout, "%d", variant_indices[i]);
            }
            std::fprintf(stdout, "\n");
        }
        if (variant_indices.size() <= 1) continue; // not multi-allelic

        const std::vector<int>& kept_indices = variant_indices;

        // Split groups larger than SUPERSITE_MAX_ALTS into deterministic chunks
        for (size_t chunk_start = 0; chunk_start < kept_indices.size(); chunk_start += SUPERSITE_MAX_ALTS) {
            uint8_t n_alts = static_cast<uint8_t>(std::min<size_t>(SUPERSITE_MAX_ALTS, kept_indices.size() - chunk_start));
            if (n_alts <= 1) continue; // chunk degenerates to a single split

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

            // Compute rare_code_mask: count code frequencies and mark rare ones
            // Bit i = 1 means code i is rare (freq < RARE_VARIANT_FREQ)
            std::vector<uint32_t> code_counts(static_cast<size_t>(n_alts + 1), 0);
            for (unsigned long h = 0; h < H.n_hap; ++h) {
                uint8_t c = hap_codes[h];
                if (c <= n_alts) code_counts[c]++;
            }
            ss.rare_code_mask = 0;
            for (uint8_t c = 0; c <= n_alts; ++c) {
                float freq = static_cast<float>(code_counts[c]) / static_cast<float>(H.n_hap);
                if (freq < RARE_VARIANT_FREQ) {
                    ss.rare_code_mask |= static_cast<uint16_t>(1u << c);
                }
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

            // Guard: verify packed codes match hap_codes when enabled
            if (supersite_debug::guard_checks_enabled()) {
                const uint8_t* buf = packed_allele_codes_out.data();
                for (unsigned int h = 0; h < H.n_hap; ++h) {
                    const uint8_t packed_code = unpackSuperSiteCode(buf, ss.panel_offset, h);
                    if (packed_code != hap_codes[h]) {
                        std::fprintf(stderr,
                                     "[supersite-guard] packed code mismatch ss_idx=%zu hap=%u expected=%u got=%u\n",
                                     super_sites_out.size(),
                                     h,
                                     static_cast<unsigned>(hap_codes[h]),
                                     static_cast<unsigned>(packed_code));
                        // Keep going; this is diagnostic only.
                    }
                }
            }

            super_sites_out.push_back(ss);

            if (trace_bp) {
                const size_t ss_idx = super_sites_out.size() - 1;
                std::fprintf(stdout,
                             "[BP_SUPERSITE_TRACE] bp=%d ss_idx=%zu anchor_locus=%u n_alts=%u members=",
                             site_bp,
                             ss_idx,
                             ss.global_site_id,
                             static_cast<unsigned>(ss.n_alts));
                for (uint16_t ai = 0; ai < ss.n_alts; ++ai) {
                    const size_t offset = ss.var_start + ai;
                    if (offset >= super_site_var_index_out.size()) break;
                    if (ai) std::fprintf(stdout, ",");
                    std::fprintf(stdout, "%d", super_site_var_index_out[offset]);
                }
                std::fprintf(stdout, " panel_offset=%u span_bytes=%u\n",
                             ss.panel_offset,
                             ss.panel_span_bytes);
            }
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

void resolveSupersiteClasses(
	genotype& g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	uint8_t& c0_out,
	uint8_t& c1_out)
{
	// First pass: check for any missing variants
	bool has_missing = false;
	bool has_called = false;
	for (uint16_t ai = 0; ai < ss.var_count; ++ai) {
		const int v_idx = super_site_var_index[ss.var_start + ai];
		unsigned char v_byte = g.Variants[DIV2(v_idx)];
		if (VAR_GET_MIS(MOD2(v_idx), v_byte)) {
			has_missing = true;
		} else {
			has_called = true;
		}
	}

	// If mixed state detected, mark all as missing
	if (has_missing && has_called) {
		for (uint16_t ai = 0; ai < ss.var_count; ++ai) {
			const int v_idx = super_site_var_index[ss.var_start + ai];
			unsigned char v_byte = g.Variants[DIV2(v_idx)];
			VAR_SET_MIS(MOD2(v_idx), v_byte);
			g.Variants[DIV2(v_idx)] = v_byte;
		}
		// Mark anchor as missing
		unsigned char anchor_byte = g.Variants[DIV2(ss.global_site_id)];
		VAR_SET_MIS(MOD2(ss.global_site_id), anchor_byte);
		g.Variants[DIV2(ss.global_site_id)] = anchor_byte;
		c0_out = 0;
		c1_out = 0;
		return;
	}

	uint8_t c0 = 0, c1 = 0;           // 0 = REF / not yet assigned
	int missing_state = -1;           // -1 unknown, 0 = non-missing, 1 = missing

	for (uint16_t ai = 0; ai < ss.var_count; ++ai) {
		const int v_idx = super_site_var_index[ss.var_start + ai];
		unsigned char v_byte = g.Variants[DIV2(v_idx)];
		const bool is_mis = VAR_GET_MIS(MOD2(v_idx), v_byte);

		// Track missing state (now guaranteed to be consistent)
		if (missing_state == -1) {
			missing_state = is_mis ? 1 : 0;
		}
		if (is_mis) continue;

		const uint8_t alt_code = static_cast<uint8_t>(ai + 1);
		const bool h0 = VAR_GET_HAP0(MOD2(v_idx), v_byte);
		const bool h1 = VAR_GET_HAP1(MOD2(v_idx), v_byte);
		if (!h0 && !h1) continue;  // both REF at this split

		// ALT on both haps at this split: treat as single class on both haps
		if (h0 && h1) {
			if (c0 == 0 && c1 == 0) {
				c0 = c1 = alt_code;
				continue;
			}
			if (c0 == alt_code && c1 == 0) { c1 = alt_code; continue; }
			if (c1 == alt_code && c0 == 0) { c0 = alt_code; continue; }
			if (c0 == alt_code && c1 == alt_code) continue;
			// Introducing a new ALT class on both haps when classes already assigned -> conflict
			throw std::runtime_error(
				"Supersite conflict: sample " + g.name +
				" carries >2 ALT classes at supersite bp=" + std::to_string(ss.bp));
		}

		// ALT on a single hap at this split
		if (c0 == 0) {
			bool flipped = false;
			c0 = alt_code;
			if (h1 && !h0) { // move ALT to hap0
				VAR_CLR_HAP1(MOD2(v_idx), v_byte);
				VAR_SET_HAP0(MOD2(v_idx), v_byte);
				flipped = true;
			}
			if (flipped) g.Variants[DIV2(v_idx)] = v_byte;
			continue;
		}

		if (c1 == 0) {
			bool flipped = false;
			c1 = alt_code;
			if (h0 && !h1) { // move ALT to hap1
				VAR_CLR_HAP0(MOD2(v_idx), v_byte);
				VAR_SET_HAP1(MOD2(v_idx), v_byte);
				flipped = true;
			}
			if (flipped) g.Variants[DIV2(v_idx)] = v_byte;
			continue;
		}

		// Third distinct ALT class -> conflict
		if (alt_code != c0 && alt_code != c1) {
			throw std::runtime_error(
				"Supersite conflict: sample " + g.name +
				" carries >2 ALT classes at supersite bp=" + std::to_string(ss.bp));
		}
	}

	// Return c0/c1 (unordered; caller may canonicalize)
	c0_out = c0;
	c1_out = c1;

	// Anchor genotype flag: set HET if distinct classes, HOM otherwise
	unsigned char anchor_byte = g.Variants[DIV2(ss.global_site_id)];
	if (missing_state == 1) {
		VAR_SET_MIS(MOD2(ss.global_site_id), anchor_byte);
	} else if (c0 != c1) {
		VAR_SET_HET(MOD2(ss.global_site_id), anchor_byte);
	} else {
		VAR_SET_HOM(MOD2(ss.global_site_id), anchor_byte);
	}
	g.Variants[DIV2(ss.global_site_id)] = anchor_byte;
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

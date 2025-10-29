#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
    std::cout << "Testing super-site builder...\n";

    // Build a variant_map with 3 split biallelic records at the same position (1:100)
    variant_map V;
    V.push(make_var("1", 100, "rs100_A_T", "A", "T", 0));
    V.push(make_var("1", 100, "rs100_A_C", "A", "C", 1));
    V.push(make_var("1", 100, "rs100_A_G", "A", "G", 2));

    // Make a small conditioning panel with 4 haplotypes
    conditioning_set H;
    H.allocate(/*n_main_samples*/0, /*n_ref_samples*/2, /*n_variants*/3); // n_hap = 4

    // Set haplotype alleles at each split record per desired codes:
    // hap0 -> ALT at v=1 (code 2), hap1 -> ALT at v=0 (code 1), hap2 -> REF (code 0), hap3 -> ALT at v=2 (code 3)
    H.H_opt_var.set(/*row=v*/1, /*col=h*/0, /*bit*/1);
    H.H_opt_var.set(/*row=v*/0, /*col=h*/1, /*bit*/1);
    // hap2 stays REF (all zeros)
    H.H_opt_var.set(/*row=v*/2, /*col=h*/3, /*bit*/1);

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> dummy;

    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, locus_to_super_idx, super_site_var_index, dummy);

    // Expect one supersite
    assert(super_sites.size() == 1);
    const SuperSite& ss = super_sites[0];
    assert(ss.n_alts == 3);
    assert(ss.var_count == 3);
    assert(ss.var_start == 0);
    assert(ss.panel_offset == 0);

    // locus_to_super_idx maps all three variants to 0
    assert(locus_to_super_idx.size() == (size_t)V.size());
    assert(locus_to_super_idx[0] == 0 && locus_to_super_idx[1] == 0 && locus_to_super_idx[2] == 0);

    // super_site_var_index contains {0,1,2}
    assert(super_site_var_index.size() == 3);
    assert(super_site_var_index[0] == 0 && super_site_var_index[1] == 1 && super_site_var_index[2] == 2);

    // Packed codes for 4 haps -> 2 bytes: [hap1:1, hap0:2] => 0x12, [hap3:3, hap2:0] => 0x30
    assert(packed_codes.size() == 2);
    assert(packed_codes[0] == (uint8_t)((1u<<4)|2u));
    assert(packed_codes[1] == (uint8_t)((3u<<4)|0u));

    // Markers
    assert(is_super_site.size() == (size_t)V.size());
    assert(is_super_site[0] && is_super_site[1] && is_super_site[2]);

    std::cout << "  OK\n";
    return 0;
}


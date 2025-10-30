/*******************************************************************************
 * Test supersite builder - grouping and panel encoding
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/variant.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
    std::cout << "Testing supersite builder..." << std::endl;
    
    // Test 1: Build a simple 2-split supersite
    variant_map V;
    V.push(make_var("1", 1000, "split1", "A", "C", 0));
    V.push(make_var("1", 1000, "split2", "A", "G", 1));
    V.push(make_var("1", 2000, "single", "T", "G", 2));
    
    // Conditioning panel: 2 reference samples => 4 haplotypes
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/2, /*n_variants*/V.size());
    
    // Set haplotype alleles for the multiallelic site
    // Hap0: ALT at split1 (carries ALT1 = code 1)
    H.H_opt_var.set(0, 0, 1);
    H.H_opt_hap.set(0, 0, 1);
    
    // Hap1: ALT at split2 (carries ALT2 = code 2)
    H.H_opt_var.set(1, 1, 1);
    H.H_opt_hap.set(1, 1, 1);
    
    // Hap2: REF at both (code 0)
    // Hap3: ALT at split1 (code 1)
    H.H_opt_var.set(0, 3, 1);
    H.H_opt_hap.set(3, 0, 1);
    
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;
    
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);
    
    // Verify results
    assert(super_sites.size() == 1);  // One supersite
    assert(super_sites[0].var_count == 2);  // Two split records
    assert(super_sites[0].n_alts == 2);  // Two alternates
    // Note: n_classes is computed as 1 + n_alts in downstream code, not set by builder
    assert(super_sites[0].global_site_id == 0);  // First variant is anchor
    
    std::cout << "  Basic supersite grouping: OK" << std::endl;
    
    // Test 2: Verify locus_to_super_idx mapping
    assert(locus_to_super_idx.size() == V.size());
    assert(locus_to_super_idx[0] == 0);  // Both splits map to supersite 0
    assert(locus_to_super_idx[1] == 0);
    assert(locus_to_super_idx[2] == -1); // Single variant not in supersite
    
    std::cout << "  Locus to supersite mapping: OK" << std::endl;
    
    // Test 3: Verify var_index array
    assert(super_site_var_index.size() == 2);  // Two member variants
    assert(super_site_var_index[0] == 0);
    assert(super_site_var_index[1] == 1);
    
    std::cout << "  Variant index array: OK" << std::endl;
    
    // Test 4: Verify panel code packing
    assert(!packed_codes.empty());
    const SuperSite& ss = super_sites[0];
    
    // Unpack and verify conditioning haplotype codes
    uint8_t code0 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 0);
    uint8_t code1 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 1);
    uint8_t code2 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 2);
    uint8_t code3 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 3);
    
    assert(code0 == 1);  // Hap0: ALT1
    assert(code1 == 2);  // Hap1: ALT2
    assert(code2 == 0);  // Hap2: REF
    assert(code3 == 1);  // Hap3: ALT1
    
    std::cout << "  Panel code packing: OK" << std::endl;
    
    // Test 5: Multiple supersites at different positions
    variant_map V2;
    V2.push(make_var("1", 1000, "ss1_a", "A", "C", 0));
    V2.push(make_var("1", 1000, "ss1_b", "A", "G", 1));
    V2.push(make_var("1", 2000, "ss2_a", "T", "A", 2));
    V2.push(make_var("1", 2000, "ss2_b", "T", "C", 3));
    V2.push(make_var("1", 2000, "ss2_c", "T", "G", 4));
    
    conditioning_set H2;
    H2.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V2.size());
    
    std::vector<SuperSite> super_sites2;
    std::vector<bool> is_super_site2;
    std::vector<uint8_t> packed_codes2;
    std::vector<int> locus_to_super_idx2;
    std::vector<int> super_site_var_index2;
    std::vector<uint8_t> unused_sample_codes2;
    
    buildSuperSites(V2, H2, super_sites2, is_super_site2, packed_codes2,
                    locus_to_super_idx2, super_site_var_index2, unused_sample_codes2);
    
    assert(super_sites2.size() == 2);  // Two supersites
    assert(super_sites2[0].var_count == 2);  // First has 2 variants
    assert(super_sites2[1].var_count == 3);  // Second has 3 variants
    assert(super_sites2[0].n_alts == 2);
    assert(super_sites2[1].n_alts == 3);
    
    std::cout << "  Multiple supersites: OK" << std::endl;
    
    // Test 6: No supersites (all single variants)
    variant_map V3;
    V3.push(make_var("1", 1000, "v1", "A", "C", 0));
    V3.push(make_var("1", 2000, "v2", "T", "G", 1));
    V3.push(make_var("1", 3000, "v3", "C", "A", 2));
    
    conditioning_set H3;
    H3.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V3.size());
    
    std::vector<SuperSite> super_sites3;
    std::vector<bool> is_super_site3;
    std::vector<uint8_t> packed_codes3;
    std::vector<int> locus_to_super_idx3;
    std::vector<int> super_site_var_index3;
    std::vector<uint8_t> unused_sample_codes3;
    
    buildSuperSites(V3, H3, super_sites3, is_super_site3, packed_codes3,
                    locus_to_super_idx3, super_site_var_index3, unused_sample_codes3);
    
    assert(super_sites3.empty());  // No supersites
    assert(locus_to_super_idx3[0] == -1);
    assert(locus_to_super_idx3[1] == -1);
    assert(locus_to_super_idx3[2] == -1);
    
    std::cout << "  No supersites case: OK" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

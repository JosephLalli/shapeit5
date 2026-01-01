/*******************************************************************************
 * Test supersite builder - grouping and panel encoding
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>
#include <string>

#include "test_common.h"
#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/variant.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

int main() {
    TEST_INIT("test_supersite_builder");

    // Test 1: Build a simple 2-split supersite
    TEST_START("basic_supersite_grouping", "Build 2-split supersite");
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
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index);

    bool test1_pass = (super_sites.size() == 1 &&
                       super_sites[0].var_count == 2 &&
                       super_sites[0].n_alts == 2 &&
                       super_sites[0].global_site_id == 0);
    TEST_CHECK(test1_pass, "basic_supersite_grouping",
               test1_pass ? "" : "Unexpected supersite structure");

    // Test 2: Verify locus_to_super_idx mapping
    TEST_START("locus_to_supersite_mapping", "Verify locus-to-supersite index mapping");
    bool test2_pass = (locus_to_super_idx.size() == V.size() &&
                       locus_to_super_idx[0] == 0 &&
                       locus_to_super_idx[1] == 0 &&
                       locus_to_super_idx[2] == -1);
    TEST_CHECK(test2_pass, "locus_to_supersite_mapping",
               test2_pass ? "" : "Incorrect locus mapping");

    // Test 3: Verify var_index array
    TEST_START("variant_index_array", "Verify variant index array");
    bool test3_pass = (super_site_var_index.size() == 2 &&
                       super_site_var_index[0] == 0 &&
                       super_site_var_index[1] == 1);
    TEST_CHECK(test3_pass, "variant_index_array",
               test3_pass ? "" : "Incorrect variant indices");

    // Test 4: Verify panel code packing
    TEST_START("panel_code_packing", "Verify panel haplotype code packing");
    const SuperSite& ss = super_sites[0];
    uint8_t code0 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 0);
    uint8_t code1 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 1);
    uint8_t code2 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 2);
    uint8_t code3 = unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 3);

    bool test4_pass = (!packed_codes.empty() &&
                       code0 == 1 && code1 == 2 && code2 == 0 && code3 == 1);
    TEST_CHECK(test4_pass, "panel_code_packing",
               test4_pass ? "" : "Incorrect panel codes");
    
    // Test 5: Multiple supersites at different positions
    TEST_START("multiple_supersites", "Build multiple supersites at different positions");
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
    buildSuperSites(V2, H2, super_sites2, is_super_site2, packed_codes2,
                    locus_to_super_idx2, super_site_var_index2);

    bool test5_pass = (super_sites2.size() == 2 &&
                       super_sites2[0].var_count == 2 &&
                       super_sites2[1].var_count == 3 &&
                       super_sites2[0].n_alts == 2 &&
                       super_sites2[1].n_alts == 3);
    TEST_CHECK(test5_pass, "multiple_supersites",
               test5_pass ? "" : "Incorrect multiple supersite structure");

    // Test 6: No supersites (all single variants)
    TEST_START("no_supersites_case", "Handle case with no multiallelic sites");
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
    buildSuperSites(V3, H3, super_sites3, is_super_site3, packed_codes3,
                    locus_to_super_idx3, super_site_var_index3);

    bool test6_pass = (super_sites3.empty() &&
                       locus_to_super_idx3[0] == -1 &&
                       locus_to_super_idx3[1] == -1 &&
                       locus_to_super_idx3[2] == -1);
    TEST_CHECK(test6_pass, "no_supersites_case",
               test6_pass ? "" : "Should have no supersites for biallelic-only input");

    TEST_SUMMARY();
    return 0;
}

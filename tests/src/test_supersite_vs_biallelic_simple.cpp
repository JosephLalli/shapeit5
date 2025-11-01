/*******************************************************************************
 * Test that buildSuperSites correctly distinguishes biallelic vs multiallelic
 * Ensures supersite representation doesn't affect biallelic sites
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/objects/variant.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/haplotype_set.h"

int main() {
    std::cout << "Testing buildSuperSites behavior..." << std::endl;
    
    // Test scenario: 10-variant context
    // - 8 biallelic variants at unique positions (should NOT become supersites)
    // - 1 multiallelic site (2 ALTs = 2 split records at same position)
    
    variant_map V;
    haplotype_set H;
    
    // Create 10 variants total
    // Variants 0-7: biallelic at unique positions
    // Variants 8-9: multiallelic (both at position 9000)
    
    std::vector<variant*> variants;
    
    // Biallelic variants at unique positions 1000-8000
    for (int i = 0; i < 8; i++) {
        variant* v = new variant();
        v->chr = 1;
        v->bp = 1000 + i * 1000;  // positions 1000, 2000, ..., 8000
        v->cm = 0.01 * (i + 1);
        v->ref = "A";
        v->alt = "T";
        v->idx = i;
        variants.push_back(v);
    }
    
    // Multiallelic site: 2 split records at position 9000
    // Split 0: REF=A, ALT=T (represents first ALT)
    variant* v8 = new variant();
    v8->chr = 1;
    v8->bp = 9000;
    v8->cm = 0.09;
    v8->ref = "A";
    v8->alt = "T";
    v8->idx = 8;
    variants.push_back(v8);
    
    // Split 1: REF=A, ALT=G (represents second ALT)
    variant* v9 = new variant();
    v9->chr = 1;
    v9->bp = 9000;  // Same position!
    v9->cm = 0.09;
    v9->ref = "A";
    v9->alt = "G";
    v9->idx = 9;
    variants.push_back(v9);
    
    // Initialize variant_map
    V.vec_pos.resize(10);
    for (int i = 0; i < 10; i++) {
        V.vec_pos[i] = variants[i];
    }
    
    // Initialize haplotype_set with dummy data (4 haplotypes)
    H.n_site = 10;
    H.n_haps = 4;
    
    // Create dummy haplotype data (not critical for this test)
    // Just need to have valid bitvector structures
    H.H_opt_hap.resize(10 * 4);  // 4 haplotypes × 10 sites
    for (size_t i = 0; i < H.H_opt_hap.size(); i++) {
        H.H_opt_hap[i] = (i % 2) ? true : false;
    }
    
    // Build supersites
    std::cout << "  Building supersites from 10-variant context..." << std::endl;
    std::vector<SuperSite> super_sites;
    std::vector<int> locus_to_super_idx;
    std::vector<bool> is_super_site;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> packed_allele_codes;
    
    buildSuperSites(V, H, super_sites, locus_to_super_idx, 
                    is_super_site, super_site_var_index, packed_allele_codes);
    
    std::cout << "  Verifying results..." << std::endl;
    
    // Verify: locus_to_super_idx should have size 10
    assert(locus_to_super_idx.size() == 10);
    std::cout << "    ✓ locus_to_super_idx has correct size (10)" << std::endl;
    
    // Verify: Variants 0-7 should NOT be in a supersite (index = -1)
    for (int i = 0; i < 8; i++) {
        assert(locus_to_super_idx[i] == -1);
    }
    std::cout << "    ✓ Biallelic variants (0-7) not in supersites" << std::endl;
    
    // Verify: Variants 8-9 should be in the SAME supersite
    assert(locus_to_super_idx[8] == locus_to_super_idx[9]);
    assert(locus_to_super_idx[8] >= 0);
    std::cout << "    ✓ Multiallelic variants (8-9) mapped to same supersite" << std::endl;
    
    // Verify: Exactly 1 supersite created
    assert(super_sites.size() == 1);
    std::cout << "    ✓ Exactly 1 supersite created" << std::endl;
    
    // Verify: Supersite has correct properties
    const SuperSite& ss = super_sites[0];
    assert(ss.chr == 1);
    assert(ss.bp == 9000);
    assert(ss.n_alts == 2);
    assert(ss.var_count == 2);
    assert(ss.global_site_id == 8);  // Anchor is first split
    std::cout << "    ✓ Supersite has correct properties (chr=1, bp=9000, n_alts=2)" << std::endl;
    
    // Verify: super_site_var_index contains indices 8 and 9
    assert(super_site_var_index.size() >= 2);
    assert(super_site_var_index[ss.var_start] == 8);
    assert(super_site_var_index[ss.var_start + 1] == 9);
    std::cout << "    ✓ Variant index array contains correct indices (8, 9)" << std::endl;
    
    // Cleanup
    for (auto* v : variants) {
        delete v;
    }
    
    std::cout << "✓ SUCCESS: buildSuperSites correctly distinguishes biallelic vs multiallelic" << std::endl;
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

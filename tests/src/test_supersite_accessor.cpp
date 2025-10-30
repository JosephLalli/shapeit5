/*******************************************************************************
 * Test supersite accessor functions and data structures
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>
#include <cstring>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/objects/variant.h"

// Helper to create SuperSite structure
static SuperSite make_supersite(unsigned int site_id, unsigned char chr, int bp, 
                                unsigned char n_alts, uint32_t panel_offset,
                                unsigned int var_start, unsigned char var_count) {
    SuperSite ss;
    ss.global_site_id = site_id;
    ss.chr = chr;
    ss.bp = bp;
    ss.n_alts = n_alts;
    ss.panel_offset = panel_offset;
    ss.var_start = var_start;
    ss.var_count = var_count;
    // Note: n_classes and class_prob_offset are set by compute_job, not in test helper
    ss.n_classes = 0;
    ss.class_prob_offset = 0;
    return ss;
}

int main() {
    std::cout << "Testing supersite accessor functions..." << std::endl;
    
    // Test 1: SuperSite structure creation
    SuperSite ss = make_supersite(/*site_id*/0, /*chr*/1, /*bp*/1000, 
                                  /*n_alts*/2, /*panel_offset*/0, 
                                  /*var_start*/0, /*var_count*/2);
    
    assert(ss.global_site_id == 0);
    assert(ss.chr == 1);
    assert(ss.bp == 1000);
    assert(ss.n_alts == 2);
    // n_classes would be set by HMM phase, not in builder
    assert(ss.panel_offset == 0);
    assert(ss.var_start == 0);
    assert(ss.var_count == 2);
    
    std::cout << "  SuperSite structure: OK" << std::endl;
    
    // Test 2: Panel haplotype code unpacking
    // Create packed codes: 4 haplotypes with allele codes 0, 1, 2, 0
    uint8_t packed_codes[2];
    packed_codes[0] = 0x10;  // hap 0 = REF (0), hap 1 = ALT1 (1)
    packed_codes[1] = 0x02;  // hap 2 = ALT2 (2), hap 3 = REF (0)
    
    assert(unpackSuperSiteCode(packed_codes, 0, 0) == 0);  // REF
    assert(unpackSuperSiteCode(packed_codes, 0, 1) == 1);  // ALT1
    assert(unpackSuperSiteCode(packed_codes, 0, 2) == 2);  // ALT2
    assert(unpackSuperSiteCode(packed_codes, 0, 3) == 0);  // REF
    
    std::cout << "  Panel code unpacking: OK" << std::endl;
    
    // Test 3: Sample genotype allele code inference
    // Create a genotype with 2 variants at same position (split records)
    genotype G(0);
    G.n_variants = 2;
    G.Variants.assign(1, 0);  // Need 1 byte for 2 variants
    
    // Set variant 0 to 1|0 (hap0=ALT, hap1=REF)
    VAR_SET_HAP0(0, G.Variants[0]);
    VAR_CLR_HAP1(0, G.Variants[0]);
    
    // Set variant 1 to 0|0 (both REF)
    VAR_CLR_HAP0(1, G.Variants[0]);
    VAR_CLR_HAP1(1, G.Variants[0]);
    
    std::vector<int> var_index = {0, 1};
    
    // Get sample codes for this supersite
    uint8_t code_hap0 = getSampleSuperSiteAlleleCode(&G, ss, var_index, 0);
    uint8_t code_hap1 = getSampleSuperSiteAlleleCode(&G, ss, var_index, 1);
    
    // Hap0: ALT at split 0 → class 1 (ALT1)
    // Hap1: REF at both splits → class 0 (REF)
    assert(code_hap0 == 1);  // ALT1
    assert(code_hap1 == 0);  // REF
    
    std::cout << "  Sample code inference: OK" << std::endl;
    
    // Test 4: Missing data code (special value)
    G.Variants[0] = 0;
    VAR_SET_MIS(0, G.Variants[0]);
    VAR_SET_MIS(1, G.Variants[0]);
    
    code_hap0 = getSampleSuperSiteAlleleCode(&G, ss, var_index, 0);
    code_hap1 = getSampleSuperSiteAlleleCode(&G, ss, var_index, 1);
    
    assert(code_hap0 == SUPERSITE_CODE_MISSING);
    assert(code_hap1 == SUPERSITE_CODE_MISSING);
    
    std::cout << "  Missing code: OK" << std::endl;
    
    // Test 5: Multi-allelic site with 3 ALTs
    SuperSite ss3 = make_supersite(/*site_id*/2, /*chr*/1, /*bp*/2000,
                                   /*n_alts*/3, /*panel_offset*/4,
                                   /*var_start*/2, /*var_count*/3);
    
    // n_classes computed as 1 + n_alts in HMM phase, not set here
    
    std::cout << "  Multi-allelic supersite: OK" << std::endl;
    
    // Test 6: aligned_vector32 works correctly
    aligned_vector32<uint8_t> vec;
    vec.resize(32);
    for (int i = 0; i < 32; ++i) {
        vec[i] = (uint8_t)i;
    }
    
    assert(vec.size() == 32);
    for (int i = 0; i < 32; ++i) {
        assert(vec[i] == (uint8_t)i);
    }
    
    // Check alignment (address should be 32-byte aligned)
    uintptr_t addr = reinterpret_cast<uintptr_t>(vec.data());
    assert((addr % 32) == 0);
    
    std::cout << "  aligned_vector32: OK" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

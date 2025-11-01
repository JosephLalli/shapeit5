/*******************************************************************************
 * Test REAL supersite emission computation with actual supersites
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/models/super_site_emissions.h"

int main() {
    std::cout << "Testing REAL supersite emission computation..." << std::endl;
    
    // Test 1: Create a 2-ALT supersite with 4 conditioning haplotypes
    // Panel codes: [0, 1, 2, 0] representing REF, ALT1, ALT2, REF
    uint8_t panel_codes[4] = {0, 1, 2, 0};
    const uint32_t n_haps = 4;
    
    aligned_vector32<double> emissions;
    emissions.resize(n_haps);
    
    // Test 1a: Sample has REF (code 0) - should match haps 0 and 3
    std::cout << "  Test 1a: Sample=REF, expect matches at haps 0,3..." << std::endl;
    precomputeSuperSiteEmissions_AVX2(panel_codes, n_haps, 0, 1.0, 0.01, emissions);
    
    assert(std::fabs(emissions[0] - 1.0) < 1e-9);   // REF matches REF
    assert(std::fabs(emissions[1] - 0.01) < 1e-9);  // ALT1 mismatches REF
    assert(std::fabs(emissions[2] - 0.01) < 1e-9);  // ALT2 mismatches REF
    assert(std::fabs(emissions[3] - 1.0) < 1e-9);   // REF matches REF
    std::cout << "    OK: emissions = [1.0, 0.01, 0.01, 1.0]" << std::endl;
    
    // Test 1b: Sample has ALT1 (code 1) - should match hap 1
    std::cout << "  Test 1b: Sample=ALT1, expect match at hap 1..." << std::endl;
    precomputeSuperSiteEmissions_AVX2(panel_codes, n_haps, 1, 1.0, 0.01, emissions);
    
    assert(std::fabs(emissions[0] - 0.01) < 1e-9);  // REF mismatches ALT1
    assert(std::fabs(emissions[1] - 1.0) < 1e-9);   // ALT1 matches ALT1
    assert(std::fabs(emissions[2] - 0.01) < 1e-9);  // ALT2 mismatches ALT1
    assert(std::fabs(emissions[3] - 0.01) < 1e-9);  // REF mismatches ALT1
    std::cout << "    OK: emissions = [0.01, 1.0, 0.01, 0.01]" << std::endl;
    
    // Test 1c: Sample has ALT2 (code 2) - should match hap 2
    std::cout << "  Test 1c: Sample=ALT2, expect match at hap 2..." << std::endl;
    precomputeSuperSiteEmissions_AVX2(panel_codes, n_haps, 2, 1.0, 0.01, emissions);
    
    assert(std::fabs(emissions[0] - 0.01) < 1e-9);  // REF mismatches ALT2
    assert(std::fabs(emissions[1] - 0.01) < 1e-9);  // ALT1 mismatches ALT2
    assert(std::fabs(emissions[2] - 1.0) < 1e-9);   // ALT2 matches ALT2
    assert(std::fabs(emissions[3] - 0.01) < 1e-9);  // REF mismatches ALT2
    std::cout << "    OK: emissions = [0.01, 0.01, 1.0, 0.01]" << std::endl;
    
    // Test 2: Different error rate
    std::cout << "  Test 2: Different mismatch probability (0.1)..." << std::endl;
    precomputeSuperSiteEmissions_AVX2(panel_codes, n_haps, 0, 1.0, 0.1, emissions);
    
    assert(std::fabs(emissions[0] - 1.0) < 1e-9);   // Match
    assert(std::fabs(emissions[1] - 0.1) < 1e-9);   // Mismatch with higher error rate
    assert(std::fabs(emissions[2] - 0.1) < 1e-9);   // Mismatch with higher error rate
    assert(std::fabs(emissions[3] - 1.0) < 1e-9);   // Match
    std::cout << "    OK: emissions = [1.0, 0.1, 0.1, 1.0]" << std::endl;
    
    // Test 3: Larger panel (8 haplotypes) - tests AVX2 vectorization fully
    std::cout << "  Test 3: 8 conditioning haplotypes (full AVX2 vector)..." << std::endl;
    uint8_t panel_codes_8[8] = {0, 1, 2, 0, 1, 2, 0, 1};
    aligned_vector32<double> emissions_8;
    emissions_8.resize(8);
    
    precomputeSuperSiteEmissions_AVX2(panel_codes_8, 8, 1, 1.0, 0.01, emissions_8);
    
    // Sample=ALT1 should match haps 1, 4, 7
    assert(std::fabs(emissions_8[0] - 0.01) < 1e-9);  // REF
    assert(std::fabs(emissions_8[1] - 1.0) < 1e-9);   // ALT1 - match
    assert(std::fabs(emissions_8[2] - 0.01) < 1e-9);  // ALT2
    assert(std::fabs(emissions_8[3] - 0.01) < 1e-9);  // REF
    assert(std::fabs(emissions_8[4] - 1.0) < 1e-9);   // ALT1 - match
    assert(std::fabs(emissions_8[5] - 0.01) < 1e-9);  // ALT2
    assert(std::fabs(emissions_8[6] - 0.01) < 1e-9);  // REF
    assert(std::fabs(emissions_8[7] - 1.0) < 1e-9);   // ALT1 - match
    std::cout << "    OK: 8-haplotype emissions computed correctly" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

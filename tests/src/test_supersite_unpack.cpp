/*******************************************************************************
 * Test supersite code unpacking (4-bit allele codes)
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <cstring>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/models/super_site_accessor.h"

int main() {
    std::cout << "Testing supersite code unpacking..." << std::endl;
    
    // Test 1: Single code unpacking
    // Pack format: lower 4 bits = even haplotype, upper 4 bits = odd haplotype
    uint8_t packed_buffer[4];
    packed_buffer[0] = 0x10;  // hap 0 = 0x0, hap 1 = 0x1
    packed_buffer[1] = 0x32;  // hap 2 = 0x2, hap 3 = 0x3
    packed_buffer[2] = 0x54;  // hap 4 = 0x4, hap 5 = 0x5
    packed_buffer[3] = 0x76;  // hap 6 = 0x6, hap 7 = 0x7
    
    assert(unpackSuperSiteCode(packed_buffer, 0, 0) == 0x0);
    assert(unpackSuperSiteCode(packed_buffer, 0, 1) == 0x1);
    assert(unpackSuperSiteCode(packed_buffer, 0, 2) == 0x2);
    assert(unpackSuperSiteCode(packed_buffer, 0, 3) == 0x3);
    assert(unpackSuperSiteCode(packed_buffer, 0, 4) == 0x4);
    assert(unpackSuperSiteCode(packed_buffer, 0, 5) == 0x5);
    assert(unpackSuperSiteCode(packed_buffer, 0, 6) == 0x6);
    assert(unpackSuperSiteCode(packed_buffer, 0, 7) == 0x7);
    
    std::cout << "  Single code unpacking: OK" << std::endl;
    
    // Test 2: Batch unpacking
    uint8_t codes_out[8];
    unpackSuperSiteCodesBatch(packed_buffer, 0, 0, 8, codes_out);
    
    for (int i = 0; i < 8; ++i) {
        assert(codes_out[i] == i);
    }
    std::cout << "  Batch code unpacking: OK" << std::endl;
    
    // Test 3: Vectorized unpacking
    uint8_t codes_vec[8];
    unpackSuperSiteCodesVectorized_PEXT(packed_buffer, 0, 8, codes_vec);
    
    for (int i = 0; i < 8; ++i) {
        assert(codes_vec[i] == i);
    }
    std::cout << "  Vectorized code unpacking: OK" << std::endl;
    
    // Test 4: Test with offset
    uint8_t packed_buffer2[8];
    for (int i = 0; i < 8; ++i) {
        packed_buffer2[i] = (uint8_t)((i * 2 + 1) << 4 | (i * 2));
    }
    
    // Offset by 4 bytes (8 haplotypes)
    assert(unpackSuperSiteCode(packed_buffer2, 4, 0) == 8);
    assert(unpackSuperSiteCode(packed_buffer2, 4, 1) == 9);
    
    std::cout << "  Offset unpacking: OK" << std::endl;
    
    // Test 5: REF code (0) vs ALT codes (1-15)
    uint8_t ref_alt_buffer[2];
    ref_alt_buffer[0] = 0x00;  // Both REF
    ref_alt_buffer[1] = 0xFF;  // Both max ALT (15)
    
    assert(unpackSuperSiteCode(ref_alt_buffer, 0, 0) == 0);  // REF
    assert(unpackSuperSiteCode(ref_alt_buffer, 0, 1) == 0);  // REF
    assert(unpackSuperSiteCode(ref_alt_buffer, 0, 2) == 15); // ALT15
    assert(unpackSuperSiteCode(ref_alt_buffer, 0, 3) == 15); // ALT15
    
    std::cout << "  REF vs ALT code unpacking: OK" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

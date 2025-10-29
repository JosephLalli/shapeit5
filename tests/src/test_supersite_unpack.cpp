#include <cassert>
#include <iostream>
#include <vector>

#include "../../phase_common/src/models/super_site_accessor.h"

int main() {
    std::cout << "Testing supersite code unpack...\n";

    // codes for 4 haps: [2,1,0,3] => bytes: 0x12, 0x30
    uint8_t packed[2];
    packed[0] = (uint8_t)((1u << 4) | 2u);
    packed[1] = (uint8_t)((3u << 4) | 0u);

    // Single
    assert(unpackSuperSiteCode(packed, 0, 0) == 2);
    assert(unpackSuperSiteCode(packed, 0, 1) == 1);
    assert(unpackSuperSiteCode(packed, 0, 2) == 0);
    assert(unpackSuperSiteCode(packed, 0, 3) == 3);

    // Batch
    uint8_t out[4] = {0};
    unpackSuperSiteCodesBatch(packed, 0, 0, 4, out);
    assert(out[0] == 2 && out[1] == 1 && out[2] == 0 && out[3] == 3);

    std::cout << "  OK\n";
    return 0;
}


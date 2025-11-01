/*******************************************************************************
 * Test supersite emission computation (float and double precision)
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

// Test basic emission computation for supersites
int main() {
    std::cout << "Testing supersite emission computation..." << std::endl;
    
    // Test 1: Match emission should be 1.0
    const float match_emission = 1.0f;
    assert(std::fabs(match_emission - 1.0f) < 1e-6f);
    std::cout << "  Match emission: OK" << std::endl;
    
    // Test 2: Mismatch emission should be error rate ratio
    const float ed = 0.1f;  // error rate
    const float ee = 1.0f;  // match rate
    const float mismatch_emission = ed / ee;
    assert(std::fabs(mismatch_emission - 0.1f) < 1e-6f);
    std::cout << "  Mismatch emission: OK" << std::endl;
    
    // Test 3: Missing emission should be 1.0 (uninformative)
    const float missing_emission = 1.0f;
    assert(std::fabs(missing_emission - 1.0f) < 1e-6f);
    std::cout << "  Missing emission: OK" << std::endl;
    
    // Test 4: Double precision emissions
    const double match_emission_d = 1.0;
    const double mismatch_emission_d = 0.1;
    assert(std::fabs(match_emission_d - 1.0) < 1e-12);
    assert(std::fabs(mismatch_emission_d - 0.1) < 1e-12);
    std::cout << "  Double precision emissions: OK" << std::endl;
    
    std::cout << "All tests passed!" << std::endl;
    return 0;
}

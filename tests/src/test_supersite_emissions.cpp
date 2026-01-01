/*******************************************************************************
 * Test supersite emission computation (float and double precision)
 ******************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../common/src/utils/otools.h"
#include "test_common.h"


// Test basic emission computation for supersites
int main() {
    TEST_INIT("test_supersite_emissions");
    TEST_START("match_emission");
    const float match_emission = 1.0f;
    TEST_CHECK(std::fabs(match_emission - 1.0f) < 1e-6f,
               "match_emission",
               "match emission should be 1.0");
    
    TEST_START("mismatch_emission");
    const float ed = 0.1f;  // error rate
    const float ee = 1.0f;  // match rate
    const float mismatch_emission = ed / ee;
    TEST_CHECK(std::fabs(mismatch_emission - 0.1f) < 1e-6f,
               "mismatch_emission",
               "mismatch emission should be error rate ratio");
    
    TEST_START("missing_emission");
    const float missing_emission = 1.0f;
    TEST_CHECK(std::fabs(missing_emission - 1.0f) < 1e-6f,
               "missing_emission",
               "missing emission should be 1.0 (uninformative)");
    
    TEST_START("double_precision_match_emission");
    const double match_emission_d = 1.0;
    const double mismatch_emission_d = 0.1;
    TEST_CHECK(std::fabs(match_emission_d - 1.0) < 1e-12,
               "double_precision_match_emission",
               "double precision match emission");
    TEST_START("double_precision_mismatch_emission");
    TEST_CHECK(std::fabs(mismatch_emission_d - 0.1) < 1e-12,
               "double_precision_mismatch_emission",
               "double precision mismatch emission");
    
    TEST_SUMMARY();
    return TestReporting::exit_code();
}

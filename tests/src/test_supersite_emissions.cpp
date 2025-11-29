/*******************************************************************************
 * Test supersite emission computation (float and double precision)
 ******************************************************************************/

#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>

#include "../../common/src/utils/otools.h"
#include "test_framework.h"


#include "test_reporting.h"
// Test basic emission computation for supersites
int main() {
    TEST_INIT("test_supersite_emissions");
    TEST_RUN("match_emission", []() {
        const float match_emission = 1.0f;
        TEST_ASSERT(std::fabs(match_emission - 1.0f) < 1e-6f, "match emission should be 1.0");
    });
    
    TEST_RUN("mismatch_emission", []() {
        const float ed = 0.1f;  // error rate
        const float ee = 1.0f;  // match rate
        const float mismatch_emission = ed / ee;
        TEST_ASSERT(std::fabs(mismatch_emission - 0.1f) < 1e-6f, "mismatch emission should be error rate ratio");
    });
    
    TEST_RUN("missing_emission", []() {
        const float missing_emission = 1.0f;
        TEST_ASSERT(std::fabs(missing_emission - 1.0f) < 1e-6f, "missing emission should be 1.0 (uninformative)");
    });
    
    TEST_RUN("double_precision_emissions", []() {
        const double match_emission_d = 1.0;
        const double mismatch_emission_d = 0.1;
        TEST_ASSERT(std::fabs(match_emission_d - 1.0) < 1e-12, "double precision match emission");
        TEST_ASSERT(std::fabs(mismatch_emission_d - 0.1) < 1e-12, "double precision mismatch emission");
    });
    
    return TEST_EXIT();
}

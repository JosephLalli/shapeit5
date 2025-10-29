#include <cassert>
#include <iostream>
#include <vector>

#include "../../phase_common/src/models/super_site_emissions.h"

int main() {
    std::cout << "Testing supersite emissions (float/double)...\n";

    // cond codes: 0,1,2,0,2,1,0,1
    uint8_t cond_codes[8] = {0,1,2,0,2,1,0,1};
    uint8_t sample_code = 2;
    double match_d = 1.0, mism_d = 0.1;
    float  match_f = 1.0f, mism_f = 0.1f;

    aligned_vector32<double> ed(8, 0.0);
    aligned_vector32<float>  ef(8, 0.0f);

    precomputeSuperSiteEmissions_AVX2(cond_codes, 8, sample_code, match_d, mism_d, ed);
    precomputeSuperSiteEmissions_FloatScalar(cond_codes, 8, sample_code, match_f, mism_f, ef);

    for (int i = 0; i < 8; ++i) {
        double expect = (cond_codes[i] == sample_code) ? match_d : mism_d;
        float  expectf = (cond_codes[i] == sample_code) ? match_f : mism_f;
        assert(std::abs(ed[i] - expect) < 1e-12);
        assert(std::abs(ef[i] - expectf) < 1e-6f);
    }

    std::cout << "  Emission precompute passed\n";
    return 0;
}


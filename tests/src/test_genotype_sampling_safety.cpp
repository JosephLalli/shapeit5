
#include <vector>
#include <iostream>
#include <cassert>
#include <limits>
#include "utils/otools.h"
#include "objects/genotype/genotype_header.h"
#include "test_reporting.h"

// Declare the global tools required by the SHAPEIT5 codebase
random_number_generator rng;
string_utils stb;
basic_algos alg;
verbose vrb;
timer tac;

void test_sample_forward_safety() {
    TEST_START("sample_forward_safety", "Testing sampleForward vector sizing safety");

    // 1. Instantiate a genotype object with 1 segment
    // Note: Constructor argument is 'index', not 'n_segments'.
    unsigned int n_segments = 1;
    genotype g(0); 
    g.name = "TEST_SAMPLE_SAFETY";
    g.n_segments = n_segments;
    
    // 2. Manually setup internal vectors
    // The constructor initializes these to size 0. We must resize them to n_segments.
    g.Diplotypes.resize(n_segments);
    g.Lengths.resize(n_segments, 0); // Length 0 implies no variants in segment
    g.Lengths_bio.resize(n_segments, 0);
    
    // We enable only diplotypes 0 and 1 (binary 11 = 3)
    // valid diplotypes = 2.
    g.Diplotypes[0] = 3UL; 

    // 3. Setup Transitions
    // We need probs for the 2 diplotypes.
    // We set them such that sum is 1.0.
    std::vector<double> trans_probs = {0.1, 0.9};
    std::vector<float> missing_probs; 

    // 4. Stress test
    // The fix resizes currProbs to 2 inside sampleForward.
    // rng.sample iterates up to 1.
    // It returns max 1. Dipcode 1 is valid.
    
    for (int i = 0; i < 100; ++i) {
        g.sampleForward(trans_probs, missing_probs);
    }

    TEST_PASS("sample_forward_safety");
}

int main() {
    TEST_INIT("test_genotype_sampling_safety");
    try {
        test_sample_forward_safety();
    } catch (const std::exception& e) {
        TEST_FAIL("sample_forward_safety", std::string("Exception: ") + e.what());
        return 1;
    } catch (...) {
        TEST_FAIL("sample_forward_safety", "Unknown exception");
        return 1;
    }
    TEST_SUMMARY();
    return 0;
}

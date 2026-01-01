// Test if we can create a genotype with Lengths[0] = 0

#include <iostream>
#include <vector>
#include <objects/genotype/genotype_header.h>


#include "test_common.h"
int main() {
    TEST_INIT("test_empty_segment");
    genotype G(0);

    // Test 1: One HOM variant (should create 1 segment with length 1)
    TEST_START("one_hom_variant", "Single homozygous variant creates one segment");
    G.n_variants = 1;
    G.Variants.assign(1, 0);
    VAR_SET_HOM(0, G.Variants[0]);
    G.build();

    bool test1_pass = (G.n_segments == 1 && G.Lengths[0] == 1);
    TEST_CHECK(test1_pass, "one_hom_variant",
               test1_pass ? "" : "Expected 1 segment with length 1");

    // Test 2: Four HET+SCA variants (should create segments)
    TEST_START("four_het_sca_variants", "Four HET+SCA variants trigger segmentation");
    G.n_variants = 4;
    G.Variants.assign(2, 0);
    for (int v = 0; v < 4; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);  // 0|1 phasing
    }
    G.build();

    bool test2_pass = (G.n_segments >= 1);
    for (unsigned int s = 0; s < G.n_segments; s++) {
        test2_pass = test2_pass && (G.Lengths[s] > 0);
    }
    TEST_CHECK(test2_pass, "four_het_sca_variants",
               test2_pass ? "" : "Expected all segment lengths > 0");

    TEST_SUMMARY();
    return 0;
}

// Test to catch empty segment underflow bug
// This test creates a scenario where:
// 1. A segment boundary creates a zero-length segment
// 2. Forward pass doesn't populate AlphaSumSum for that segment
// 3. Backward pass encounters zero AlphaSumSum and triggers underflow

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <objects/genotype/genotype_header.h>
#include <models/haplotype_segment_single.h>
#include <models/haplotype_segment_double.h>
#include <containers/bitmatrix.h>
#include <objects/hmm_parameters.h>


#include "test_reporting.h"
// Helper to create a minimal genetic map
void create_minimal_map(hmm_parameters &M, int n_loci, double cm_distance = 0.0001) {
    M.t.assign(n_loci, 0.0f);
    M.nt.assign(n_loci, 1.0f);
    M.cm.assign(n_loci, 0.0f);
    M.rare_allele.assign(n_loci, 0);

    for (int i = 0; i < n_loci; i++) {
        float cm_val = static_cast<float>(i * cm_distance);
        M.cm[i] = cm_val;
        M.t[i] = static_cast<float>(1.0 - std::exp(-2.0 * cm_val));
        M.nt[i] = 1.0f - M.t[i];
    }

    M.ed = 0.0001;
    M.ee = 0.9999;
    M.Neff = 100;
    M.Nhap = 100;
}

// Helper to create conditioning haplotypes
void create_conditioning_haps(bitmatrix &H, int n_samples, int n_loci) {
    H.allocate(n_samples, n_loci);
    // Fill with alternating pattern for diversity
    for (int s = 0; s < n_samples; s++) {
        for (int l = 0; l < n_loci; l++) {
            H.set(s, l, (s + l) % 2);
        }
    }
}

// Test case 1: Empty segment at the end
bool test_empty_final_segment() {
    std::cout << "\n=== Test 1: Empty Final Segment ===" << std::endl;

    genotype G(0);
    G.name = "TEST_SAMPLE";

    // Create scenario: 4 HET+SCA variants that trigger boundary, then nothing
    // This forces: segment[0] with 4 variants, then segment[1] with 0 variants
    G.n_variants = 4;
    G.Variants.assign(2, 0);

    for (int v = 0; v < 4; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);  // 0|1 phasing
    }

    G.build();

    std::cout << "  n_segments = " << G.n_segments << std::endl;
    std::cout << "  n_variants = " << G.n_variants << std::endl;

    bool has_empty_segment = false;
    for (unsigned int s = 0; s < G.n_segments; s++) {
        std::cout << "  Lengths[" << s << "] = " << G.Lengths[s] << std::endl;
        if (G.Lengths[s] == 0) {
            has_empty_segment = true;
            std::cout << "  ⚠️  FOUND EMPTY SEGMENT at index " << s << std::endl;
        }
    }

    if (has_empty_segment) {
        std::cout << "  ❌ FAIL: Empty segment created during build" << std::endl;
        return false;
    }

    std::cout << "  ✓ PASS: No empty segments" << std::endl;
    return true;
}

// Test case 2: Empty segment in middle
bool test_empty_middle_segment() {
    std::cout << "\n=== Test 2: Empty Middle Segment ===" << std::endl;

    genotype G(0);
    G.name = "TEST_SAMPLE";

    // Create: [4 HET+SCA] -> boundary -> [4 HET+SCA] -> boundary
    // Goal: force empty segment between boundaries
    G.n_variants = 8;
    G.Variants.assign(4, 0);

    // First group of 4 HET+SCA
    for (int v = 0; v < 4; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);
    }

    // Second group of 4 HET+SCA
    for (int v = 4; v < 8; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);  // Different phase
    }

    G.build();

    std::cout << "  n_segments = " << G.n_segments << std::endl;

    bool has_empty_segment = false;
    for (unsigned int s = 0; s < G.n_segments; s++) {
        std::cout << "  Lengths[" << s << "] = " << G.Lengths[s] << std::endl;
        if (G.Lengths[s] == 0) {
            has_empty_segment = true;
            std::cout << "  ⚠️  FOUND EMPTY SEGMENT at index " << s << std::endl;
        }
    }

    if (has_empty_segment) {
        std::cout << "  ❌ FAIL: Empty segment created" << std::endl;
        return false;
    }

    std::cout << "  ✓ PASS: No empty segments" << std::endl;
    return true;
}

// Test case 3: Full forward/backward pass with potential empty segment
bool test_forward_backward_with_empty_segment() {
    std::cout << "\n=== Test 3: Forward/Backward Pass with Segment Validation ===" << std::endl;

    genotype G(0);
    G.name = "TEST_SAMPLE_FB";

    // Create a realistic scenario with multiple segments
    G.n_variants = 12;
    G.Variants.assign(6, 0);

    // Pattern: 4 HET (boundary), 2 HOM, 4 HET (boundary), 2 HOM
    for (int v = 0; v < 4; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
    }
    for (int v = 4; v < 6; v++) {
        VAR_SET_HOM(MOD2(v), G.Variants[DIV2(v)]);
    }
    for (int v = 6; v < 10; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
    }
    for (int v = 10; v < 12; v++) {
        VAR_SET_HOM(MOD2(v), G.Variants[DIV2(v)]);
    }

    G.build();

    std::cout << "  n_segments = " << G.n_segments << std::endl;
    std::cout << "  n_variants = " << G.n_variants << std::endl;

    // Check for empty segments
    bool has_empty_segment = false;
    for (unsigned int s = 0; s < G.n_segments; s++) {
        std::cout << "  Lengths[" << s << "] = " << G.Lengths[s] << std::endl;
        if (G.Lengths[s] == 0) {
            has_empty_segment = true;
            std::cout << "  ⚠️  FOUND EMPTY SEGMENT at index " << s << std::endl;
        }
    }

    if (has_empty_segment) {
        std::cout << "  ❌ FAIL: Empty segment would cause underflow" << std::endl;
        return false;
    }

    // Now run actual forward/backward pass
    try {
        int n_cond_haps = 10;
        bitmatrix H;
        create_conditioning_haps(H, n_cond_haps, G.n_variants);

        std::vector<unsigned int> idxH(n_cond_haps);
        for (int i = 0; i < n_cond_haps; i++) idxH[i] = i;

        hmm_parameters M;
        create_minimal_map(M, G.n_variants);

        window W;
        W.start_segment = 0;
        W.stop_segment = G.n_segments - 1;
        W.start_locus = 0;
        W.stop_locus = G.n_variants - 1;
        W.start_ambiguous = 0;
        W.stop_ambiguous = G.n_ambiguous - 1;
        W.start_missing = 0;
        W.stop_missing = G.n_missing - 1;
        W.start_transition = 0;
        W.stop_transition = 100;  // Sufficient for test

        // Try single precision first
        haplotype_segment_single HS(&G, H, idxH, W, M);
        HS.forward();

        std::vector<double> trans_probs(200, 0.0);
        std::vector<float> miss_probs(G.n_missing, 0.0f);

        int outcome = HS.backward(trans_probs, miss_probs, nullptr, nullptr, nullptr);

        if (outcome == -1) {
            std::cout << "  ❌ FAIL: Haploid underflow (single precision)" << std::endl;

            // Try double precision
            haplotype_segment_double HD(&G, H, idxH, W, M);
            HD.forward();
            outcome = HD.backward(trans_probs, miss_probs, nullptr, nullptr, nullptr);

            if (outcome == -1) {
                std::cout << "  ❌ FAIL: Haploid underflow still occurs in double precision" << std::endl;
                return false;
            } else if (outcome == -2) {
                std::cout << "  ❌ FAIL: Diploid underflow in double precision" << std::endl;
                return false;
            } else {
                std::cout << "  ⚠️  Recovered with double precision (outcome=" << outcome << ")" << std::endl;
            }
        } else if (outcome == -2) {
            std::cout << "  ❌ FAIL: Diploid underflow" << std::endl;
            return false;
        } else {
            std::cout << "  ✓ Forward/backward completed successfully (outcome=" << outcome << ")" << std::endl;
        }

    } catch (const std::exception &e) {
        std::cout << "  ❌ EXCEPTION: " << e.what() << std::endl;
        return false;
    }

    std::cout << "  ✓ PASS: Forward/backward with segment validation" << std::endl;
    return true;
}

// Test case 4: Check for zero-length segments specifically
bool test_zero_length_segment_detection() {
    std::cout << "\n=== Test 4: Zero-Length Segment Detection ===" << std::endl;

    genotype G(0);
    G.name = "TEST_ZERO_LENGTH";

    // Manually construct a pathological case
    G.n_variants = 0;
    G.Variants.clear();
    G.build();

    std::cout << "  Zero variants case:" << std::endl;
    std::cout << "    n_segments = " << G.n_segments << std::endl;

    bool found_issue = false;
    for (unsigned int s = 0; s < G.n_segments; s++) {
        if (s < G.Lengths.size()) {
            std::cout << "    Lengths[" << s << "] = " << G.Lengths[s] << std::endl;
            if (G.Lengths[s] == 0 && G.n_segments > 1) {
                found_issue = true;
            }
        }
    }

    // Having one segment with length 0 when n_variants=0 is acceptable
    // Having multiple segments with any having length 0 is NOT acceptable
    if (found_issue) {
        std::cout << "  ❌ FAIL: Multiple segments with at least one empty" << std::endl;
        return false;
    }

    std::cout << "  ✓ PASS: Zero-length segment handling correct" << std::endl;
    return true;
}

int main() {
    TEST_INIT("test_empty_segment_underflow");
    std::cout << "╔════════════════════════════════════════════════╗" << std::endl;
    std::cout << "║  Empty Segment Underflow Detection Test       ║" << std::endl;
    std::cout << "╚════════════════════════════════════════════════╝" << std::endl;

    bool all_passed = true;

    all_passed &= test_empty_final_segment();
    all_passed &= test_empty_middle_segment();
    all_passed &= test_forward_backward_with_empty_segment();
    all_passed &= test_zero_length_segment_detection();

    std::cout << "\n" << std::string(50, '=') << std::endl;
    if (all_passed) {
        std::cout << "✓ ALL TESTS PASSED" << std::endl;
        std::cout << "\nNOTE: These tests verify basic empty segment prevention." << std::endl;
        std::cout << "The empty segment underflow bug occurs specifically with" << std::endl;
        std::cout << "SUPERSITES enabled when:" << std::endl;
        std::cout << "  1. Segment boundary occurs after processing anchor" << std::endl;
        std::cout << "  2. Sibling variants follow, which don't trigger boundaries" << std::endl;
        std::cout << "  3. Loop terminates, creating Lengths[seg] = 0" << std::endl;
        std::cout << "  4. Forward pass never updates AlphaSumSum[empty_seg]" << std::endl;
        std::cout << "  5. Backward pass reads AlphaSumSum[seg-1] = 0.0" << std::endl;
        std::cout << "  6. TRANS_HAP() detects zero and returns underflow" << std::endl;
        std::cout << "\nTo trigger this bug, run: test/scripts/phase.chr22.wgs.sh" << std::endl;
        std::cout << "with --enable-supersites on real data with complex variant" << std::endl;
        std::cout << "patterns that create empty segments." << std::endl;
        TEST_SUMMARY();
        return 0;
    } else {
        std::cout << "❌ SOME TESTS FAILED" << std::endl;
        std::cout << "\nThis indicates the empty segment bug is present." << std::endl;
        std::cout << "Empty segments cause AlphaSumSum[i] to remain 0.0," << std::endl;
        std::cout << "triggering underflow in backward pass transitions." << std::endl;
        return 1;
    }
}

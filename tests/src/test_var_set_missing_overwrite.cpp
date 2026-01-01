#include "test_common.h"
#include <objects/genotype/genotype_header.h>
#include <cstdio>
#include <iostream>

// Regression guard: ensure VAR_SET_MIS overwrites any pre-existing code bits.
// Simulate MOD2/DIV2 for e=0 (first variant in byte)
#define E 0

int main(int argc, char ** argv) {
    // 1. Set to HET (Unphased Heterozygous)
    // Binary: 0000 0010 (assuming e=0)
    unsigned char v = 0;
    VAR_SET_HET(E, v);
    
    std::cout << "Initial State (HET):" << std::endl;
    std::cout << "  Value: " << (int)v << std::endl;
    std::cout << "  Is HET? " << VAR_GET_HET(E, v) << std::endl;
    std::cout << "  Is MIS? " << VAR_GET_MIS(E, v) << std::endl;
    std::cout << "  Is SCA? " << VAR_GET_SCA(E, v) << std::endl;

    if (!VAR_GET_HET(E, v)) {
        std::cout << "test_var_set_missing_overwrite: FAIL (Setup failed)" << std::endl;
        return 1;
    }

    // 2. Apply VAR_SET_MIS (Simulate bug in resolveSupersiteClasses)
    // Goal: Change to Missing.
    // Expected: 0000 0001
    // Actual:   0000 0010 | 0000 0001 = 0000 0011 (SCA)
    VAR_SET_MIS(E, v);

    std::cout << "\nAfter VAR_SET_MIS:" << std::endl;
    std::cout << "  Value: " << (int)v << std::endl;
    std::cout << "  Is HET? " << VAR_GET_HET(E, v) << std::endl;
    std::cout << "  Is MIS? " << VAR_GET_MIS(E, v) << std::endl;
    std::cout << "  Is SCA? " << VAR_GET_SCA(E, v) << std::endl;

    if (VAR_GET_MIS(E, v) && !VAR_GET_HET(E, v) && !VAR_GET_SCA(E, v)) {
        std::cout << "\nResult Correct: Variant is MIS." << std::endl;
        std::cout << "test_var_set_missing_overwrite: PASS" << std::endl;
        return 0;
    }

    std::cout << "\nBUG REPRODUCED: Variant failed to become MIS (state unexpected)" << std::endl;
    std::cout << "test_var_set_missing_overwrite: FAIL" << std::endl;
    return 1;
}

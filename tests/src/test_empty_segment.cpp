// Test if we can create a genotype with Lengths[0] = 0

#include <iostream>
#include <vector>
#include <objects/genotype/genotype_header.h>

int main() {
    genotype G(0);

    // Test 2: One HOM variant (should create 1 segment with length 1)
    G.n_variants = 1;
    G.Variants.assign(1, 0);
    VAR_SET_HOM(0, G.Variants[0]);
    G.build();
    
    std::cout << "\nTest 2 (one HOM variant):" << std::endl;
    std::cout << "  n_segments = " << G.n_segments << std::endl;
    if (G.n_segments > 0) {
        std::cout << "  Lengths[0] = " << G.Lengths[0] << std::endl;
    }
    
    // Test 3: Four HET+SCA variants (should trigger boundary)
    G.n_variants = 4;
    G.Variants.assign(2, 0);
    for (int v = 0; v < 4; v++) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_SCA(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);  // 0|1 phasing
    }
    G.build();
    
    std::cout << "\nTest 3 (four HET+SCA variants):" << std::endl;
    std::cout << "  n_segments = " << G.n_segments << std::endl;
    for (unsigned int s = 0; s < G.n_segments; s++) {
        std::cout << "  Lengths[" << s << "] = " << G.Lengths[s] << std::endl;
    }
    
    return 0;
}

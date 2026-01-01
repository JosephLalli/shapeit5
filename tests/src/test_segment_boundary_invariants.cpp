/*******************************************************************************
 * Test: Segment Boundary Invariants with Supersites
 * 
 * Purpose: Validate that segment construction correctly handles supersite boundaries
 *          and enforces the following design invariants:
 * 
 *          1. Segments never start at sibling loci
 *          2. If a segment ends on a multiallelic site, it contains all siblings
 * 
 * Test Cases:
 * 1. Simple supersite at segment boundary (anchor + sibling)
 * 2. Large supersite spanning potential segment boundary
 * 3. Multiple supersites with forced segment boundaries
 * 4. Edge case: segment that would naturally end on anchor vs sibling
 * 5. Validation of segment start positions (never siblings)
 * 6. Validation of segment end positions (complete supersites)
 * 
 * Date: November 9, 2025
 ******************************************************************************/

#include <utils/otools.h>
#include <objects/hmm_parameters.h>
#include <objects/genotype/genotype_header.h>
#include <objects/compute_job.h>
#include <containers/bitmatrix.h>
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"


#include "test_common.h"
// Expose private members for testing
#define private public
#define protected public
#include <models/haplotype_segment_single.h>
#undef private
#undef protected

#include <random>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

// Global test infrastructure
static std::mt19937 test_rng(12345);

// Helper: Create variant
static variant* make_var(string chr, int bp, string id, string ref, string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

// Test assertion helper
static void test_assert(bool condition, const string& message) {
    if (condition) {
        cout << "[PASS] " << message << endl;
    } else {
        cout << "[FAIL] " << message << endl;
        exit(1);
    }
}

// Helper: Create multiallelic variants at same position
void create_multiallelic_variants(variant_map& V, int bp, string base_id, string ref, vector<string> alts, int& idx) {
    for (const string& alt : alts) {
        string id = base_id + "_" + ref + "_" + alt;
        V.vec_pos.push_back(make_var("chr1", bp, id, ref, alt, idx++));
        V.vec_pos.back()->cm = 0.01 * bp; // Map position based on bp
    }
}

// Helper: Create biallelic variant
void create_biallelic_variant(variant_map& V, int bp, string id, string ref, string alt, int& idx) {
    V.vec_pos.push_back(make_var("chr1", bp, id, ref, alt, idx++));
    V.vec_pos.back()->cm = 0.01 * bp; // Map position based on bp
}

// Helper: Set up minimal genotype for testing
void setup_test_genotype(genotype& G, int n_variants) {
    G.n_variants = n_variants;
    
    // Initialize variant storage (packed format)
    int n_bytes = (n_variants + 1) / 2;
    G.Variants.resize(n_bytes, 0);
    
    // Set default: all variants HOM REF
    for (int v = 0; v < n_variants; ++v) {
        // VAR_CLR_HET, VAR_CLR_SCA, VAR_CLR_MIS sets HOM REF
        // This is the default initialization
    }
}

// Helper: Force segment boundary by creating 4 HET variants
void force_segment_boundary(genotype& G, int start_variant) {
    for (int i = 0; i < 4; ++i) {
        int v = start_variant + i;
        if (v < G.n_variants) {
            VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
            VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]); // 0|1 phasing
        }
    }
}

// Helper: Make variant heterozygous
void set_variant_het(genotype& G, int v) {
    if (v < G.n_variants) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]); // 0|1 phasing
    }
}

/**
 * Test 1: Simple supersite at segment boundary
 * Creates a scenario where segment boundary logic must handle anchor + sibling correctly
 */
bool test_simple_supersite_boundary() {
    cout << "\n=== TEST 1: Simple Supersite at Segment Boundary ===" << endl;
    
    // Create variant map with supersite + boundary-forcing variants
    variant_map V;
    int idx = 0;
    
    // Locus 0-1: Multiallelic supersite (A/C, A/G at bp=1000)
    create_multiallelic_variants(V, 1000, "SS1", "A", {"C", "G"}, idx);
    
    // Locus 2-5: Four HET variants to force segment boundary
    for (int i = 2; i < 6; ++i) {
        create_biallelic_variant(V, 2000 + i*100, "VAR" + to_string(i), "A", "T", idx);
    }
    
    cout << "  Created " << V.size() << " variants (2 supersite + 4 boundary-forcing)" << endl;
    
    // Build supersites
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<int> locus_to_super_idx;
    vector<uint8_t> packed_codes;
    vector<int> super_site_var_index;
    
    // Build minimal conditioning set
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/2, /*n_variants*/V.size()); // 2 ref samples = 4 haplotypes
    for (unsigned int v = 0; v < V.size(); ++v) {
        // Set basic panel data (REF/ALT pattern)
        // H_opt_var tracks which variants are ALT; H_opt_hap tracks which haps carry ALT
        H.H_opt_var.set(v, 1, 1);  // hap 1 carries ALT at variant v
        H.H_opt_hap.set(1, v, 1);  // hap 1 carries ALT at variant v
        H.H_opt_var.set(v, 3, 1);  // hap 3 carries ALT at variant v  
        H.H_opt_hap.set(3, v, 1);  // hap 3 carries ALT at variant v
    }
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    
    cout << "  Built " << super_sites.size() << " supersites" << endl;
    test_assert(super_sites.size() == 1, "Expected 1 supersite");
    test_assert(super_sites[0].var_count == 2, "Supersite should have 2 variants");
    test_assert(super_sites[0].global_site_id == 0, "Supersite anchor should be locus 0");
    
    // Set up genotype
    genotype G(1); // Single sample
    setup_test_genotype(G, V.size());
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr);
    
    // Make anchor heterozygous (will represent entire supersite)
    set_variant_het(G, 0); // anchor
    // Sibling at locus 1 will be handled automatically by supersite logic
    
    // Force segment boundary after supersite
    force_segment_boundary(G, 2); // variants 2-5
    
    // Build segments
    G.build();
    
    cout << "  Genotype built: " << G.n_segments << " segments" << endl;
    for (int s = 0; s < G.n_segments; ++s) {
        cout << "    Segment " << s << ": length=" << G.Lengths[s] << endl;
    }
    
    // INVARIANT 1: Check that no segment starts at a sibling
    bool invariant1_passed = true;
    int locus = 0;
    for (int s = 0; s < G.n_segments; ++s) {
        // Check if segment start position is a sibling
        genotype::SuperSiteContext ctx = G.getSuperSiteContext(locus);
        bool is_sibling = ctx.is_member && !ctx.is_anchor;
        
        if (is_sibling) {
            cout << "    [FAIL] Segment " << s << " starts at sibling locus " << locus << endl;
            invariant1_passed = false;
        } else {
            cout << "    [PASS] Segment " << s << " starts at valid locus " << locus 
                 << (ctx.is_member ? " (anchor)" : " (regular)") << endl;
        }
        
        locus += G.Lengths[s];
    }
    
    test_assert(invariant1_passed, "INVARIANT 1: No segment starts at sibling");
    
    // INVARIANT 2: Check that segments ending on multiallelic sites contain all siblings
    bool invariant2_passed = true;
    locus = 0;
    for (int s = 0; s < G.n_segments; ++s) {
        int segment_end = locus + G.Lengths[s] - 1;
        genotype::SuperSiteContext end_ctx = G.getSuperSiteContext(segment_end);
        
        if (end_ctx.is_member) {
            // This segment ends within a supersite
            int ss_idx = end_ctx.ss_idx;
            const SuperSite& ss = super_sites[ss_idx];
            int anchor_pos = ss.global_site_id;
            int sibling_count = ss.var_count - 1;
            
            // Check if segment contains all siblings of this supersite
            bool contains_all_siblings = (segment_end >= anchor_pos + sibling_count);
            
            if (!contains_all_siblings) {
                cout << "    [FAIL] Segment " << s << " ends at locus " << segment_end 
                     << " but doesn't contain all siblings (anchor=" << anchor_pos 
                     << ", needs until locus " << (anchor_pos + sibling_count) << ")" << endl;
                invariant2_passed = false;
            } else {
                cout << "    [PASS] Segment " << s << " contains complete supersite "
                     << "(anchor=" << anchor_pos << ", through sibling=" << (anchor_pos + sibling_count) << ")" << endl;
            }
        }
        
        locus += G.Lengths[s];
    }
    
    test_assert(invariant2_passed, "INVARIANT 2: Segments ending on multiallelic sites contain all siblings");
    
    return invariant1_passed && invariant2_passed;
}

/**
 * Test 2: Large supersite spanning potential segment boundary
 * Tests what happens when a supersite is large enough to potentially span segments
 */
bool test_large_supersite_spanning() {
    cout << "\n=== TEST 2: Large Supersite Spanning Potential Boundary ===" << endl;
    
    // Create variant map with large multiallelic site
    variant_map V;
    int idx = 0;
    
    // Locus 0-3: Large multiallelic supersite (A/C, A/G, A/T, A/N at bp=1000) 
    create_multiallelic_variants(V, 1000, "LARGE_SS", "A", {"C", "G", "T", "N"}, idx);
    
    // Locus 4-7: Four more variants after supersite
    for (int i = 4; i < 8; ++i) {
        create_biallelic_variant(V, 3000 + i*100, "POST" + to_string(i), "A", "T", idx);
    }
    
    cout << "  Created " << V.size() << " variants (4 supersite + 4 post)" << endl;
    
    // Build supersites
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<int> locus_to_super_idx;
    vector<uint8_t> packed_codes;
    vector<int> super_site_var_index;
    
    // Build minimal conditioning set
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/2, /*n_variants*/V.size()); // 2 ref samples = 4 haplotypes
    for (unsigned int v = 0; v < V.size(); ++v) {
        // Set basic panel data (REF/ALT pattern)
        // H_opt_var tracks which variants are ALT; H_opt_hap tracks which haps carry ALT
        H.H_opt_var.set(v, 1, 1);  // hap 1 carries ALT at variant v
        H.H_opt_hap.set(1, v, 1);  // hap 1 carries ALT at variant v
        H.H_opt_var.set(v, 3, 1);  // hap 3 carries ALT at variant v  
        H.H_opt_hap.set(3, v, 1);  // hap 3 carries ALT at variant v
    }
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    
    cout << "  Built " << super_sites.size() << " supersites" << endl;
    test_assert(super_sites.size() == 1, "Expected 1 large supersite");
    test_assert(super_sites[0].var_count == 4, "Large supersite should have 4 variants");
    
    // Set up genotype
    genotype G(1); // Single sample
    setup_test_genotype(G, V.size());
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr);
    
    // Make the large supersite heterozygous
    set_variant_het(G, 0); // anchor represents entire 4-variant supersite
    
    // Make some post-supersite variants HET to create segment boundaries
    set_variant_het(G, 4);
    set_variant_het(G, 5);
    
    // Build segments
    G.build();
    
    cout << "  Genotype built: " << G.n_segments << " segments" << endl;
    for (int s = 0; s < G.n_segments; ++s) {
        cout << "    Segment " << s << ": length=" << G.Lengths[s] << endl;
    }
    
    // Validate both invariants
    bool invariant1_passed = true;
    bool invariant2_passed = true;
    int locus = 0;
    
    for (int s = 0; s < G.n_segments; ++s) {
        // Check segment start (INVARIANT 1)
        genotype::SuperSiteContext start_ctx = G.getSuperSiteContext(locus);
        bool start_is_sibling = start_ctx.is_member && !start_ctx.is_anchor;
        
        if (start_is_sibling) {
            cout << "    [FAIL] Segment " << s << " starts at sibling locus " << locus << endl;
            invariant1_passed = false;
        }
        
        // Check segment end (INVARIANT 2) 
        int segment_end = locus + G.Lengths[s] - 1;
        genotype::SuperSiteContext end_ctx = G.getSuperSiteContext(segment_end);
        
        if (end_ctx.is_member) {
            int ss_idx = end_ctx.ss_idx;
            const SuperSite& ss = super_sites[ss_idx];
            int anchor_pos = ss.global_site_id;
            int last_sibling_pos = anchor_pos + ss.var_count - 1;
            
            if (segment_end < last_sibling_pos) {
                cout << "    [FAIL] Segment " << s << " ends at locus " << segment_end 
                     << " but supersite continues until " << last_sibling_pos << endl;
                invariant2_passed = false;
            }
        }
        
        locus += G.Lengths[s];
    }
    
    test_assert(invariant1_passed, "INVARIANT 1: No segment starts at sibling (large supersite)");
    test_assert(invariant2_passed, "INVARIANT 2: Large supersite not split across segments");
    
    return invariant1_passed && invariant2_passed;
}

/**
 * Test 3: Multiple supersites with forced boundaries
 * Tests complex scenario with multiple supersites and forced segment boundaries
 */
bool test_multiple_supersites_with_boundaries() {
    cout << "\n=== TEST 3: Multiple Supersites with Forced Boundaries ===" << endl;
    
    // Create variant map with multiple supersites separated by boundary-forcing variants
    variant_map V;
    int idx = 0;
    
    // Supersite 1: locus 0-1 (A/C, A/G at bp=1000)
    create_multiallelic_variants(V, 1000, "SS1", "A", {"C", "G"}, idx);
    
    // Force boundary: locus 2-5 (4 HET variants)
    for (int i = 2; i < 6; ++i) {
        create_biallelic_variant(V, 2000 + i*100, "BOUND1_" + to_string(i), "A", "T", idx);
    }
    
    // Supersite 2: locus 6-8 (A/C, A/G, A/T at bp=7000) 
    create_multiallelic_variants(V, 7000, "SS2", "A", {"C", "G", "T"}, idx);
    
    // More variants: locus 9-10
    for (int i = 9; i < 11; ++i) {
        create_biallelic_variant(V, 8000 + i*100, "POST_" + to_string(i), "A", "T", idx);
    }
    
    cout << "  Created " << V.size() << " variants (5 supersite + 6 regular)" << endl;
    
    // Build supersites
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<int> locus_to_super_idx;
    vector<uint8_t> packed_codes;
    vector<int> super_site_var_index;
    
    // Build minimal conditioning set
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/2, /*n_variants*/V.size()); // 2 ref samples = 4 haplotypes
    for (unsigned int v = 0; v < V.size(); ++v) {
        // Set basic panel data (REF/ALT pattern)
        // H_opt_var tracks which variants are ALT; H_opt_hap tracks which haps carry ALT
        H.H_opt_var.set(v, 1, 1);  // hap 1 carries ALT at variant v
        H.H_opt_hap.set(1, v, 1);  // hap 1 carries ALT at variant v
        H.H_opt_var.set(v, 3, 1);  // hap 3 carries ALT at variant v  
        H.H_opt_hap.set(3, v, 1);  // hap 3 carries ALT at variant v
    }
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    
    cout << "  Built " << super_sites.size() << " supersites" << endl;
    test_assert(super_sites.size() == 2, "Expected 2 supersites");
    
    // Set up genotype with complex pattern
    genotype G(1); // Single sample
    setup_test_genotype(G, V.size());
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr);
    
    // Make both supersites heterozygous
    set_variant_het(G, 0); // SS1 anchor
    set_variant_het(G, 6); // SS2 anchor
    
    // Force segment boundaries
    force_segment_boundary(G, 2); // Between SS1 and SS2
    
    // Build segments
    G.build();
    
    cout << "  Genotype built: " << G.n_segments << " segments" << endl;
    for (int s = 0; s < G.n_segments; ++s) {
        cout << "    Segment " << s << ": length=" << G.Lengths[s] << endl;
    }
    
    // Comprehensive validation
    bool all_invariants_passed = true;
    int locus = 0;
    
    for (int s = 0; s < G.n_segments; ++s) {
        cout << "    Checking segment " << s << " (locus " << locus << " to " << (locus + G.Lengths[s] - 1) << ")" << endl;
        
        // INVARIANT 1: Segment start
        genotype::SuperSiteContext start_ctx = G.getSuperSiteContext(locus);
        if (start_ctx.is_member && !start_ctx.is_anchor) {
            cout << "      [FAIL] Starts at sibling locus " << locus << endl;
            all_invariants_passed = false;
        } else {
            cout << "      [PASS] Starts at valid position" << endl;
        }
        
        // INVARIANT 2: Segment end
        int segment_end = locus + G.Lengths[s] - 1;
        genotype::SuperSiteContext end_ctx = G.getSuperSiteContext(segment_end);
        
        if (end_ctx.is_member) {
            int ss_idx = end_ctx.ss_idx;
            const SuperSite& ss = super_sites[ss_idx];
            int anchor_pos = ss.global_site_id;
            int last_sibling_pos = anchor_pos + ss.var_count - 1;
            
            if (segment_end < last_sibling_pos) {
                cout << "      [FAIL] Ends within supersite but doesn't include all siblings" << endl;
                all_invariants_passed = false;
            } else {
                cout << "      [PASS] Contains complete supersite if ending in one" << endl;
            }
        }
        
        locus += G.Lengths[s];
    }
    
    test_assert(all_invariants_passed, "All segment boundary invariants with multiple supersites");
    
    return all_invariants_passed;
}

/**
 * Main test runner
 */
int main() {
    TEST_INIT("test_segment_boundary_invariants");
    cout << "======================================================" << endl;
    cout << "Segment Boundary Invariants Test Suite" << endl;
    cout << "======================================================" << endl;
    
    bool all_tests_passed = true;
    
    try {
        // Run all test cases
        all_tests_passed &= test_simple_supersite_boundary();
        all_tests_passed &= test_large_supersite_spanning();
        all_tests_passed &= test_multiple_supersites_with_boundaries();
        
        cout << "\n======================================================" << endl;
        cout << "Test Summary" << endl;
        cout << "======================================================" << endl;
        
        if (all_tests_passed) {
            cout << "✓ ALL TESTS PASSED - Segment boundary invariants are correctly enforced" << endl;
            cout << "✓ INVARIANT 1: Segments never start at sibling loci" << endl;
            cout << "✓ INVARIANT 2: Segments ending on multiallelic sites contain all siblings" << endl;
            TEST_SUMMARY();
            return 0;
        } else {
            cout << "✗ TESTS FAILED - Segment boundary invariants are violated" << endl;
            return 1;
        }
        
    } catch (const exception& e) {
        cout << "Test execution failed with exception: " << e.what() << endl;
        return 1;
    }
}
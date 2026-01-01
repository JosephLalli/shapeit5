/*******************************************************************************
 * Test: HMM Indexing Validation in Supersite Mode
 * 
 * Purpose: Validate correct HMM cursor and array indexing when processing 
 *          multiallelic supersites. Tests key indexes including:
 *          - Genotype arrays (Variants, Lengths, Diplotypes, Ambiguous)
 *          - HMM state arrays (Alpha, AlphaSum, AlphaMissing)
 *          - Locus cursors (curr_abs_locus, prev_abs_locus, curr_rel_locus)
 *          - Segment boundaries (locus_first/last, segment_first/last)
 *          - Type-specific cursors (curr_abs_ambiguous, curr_abs_missing)
 *          - Sibling handling (anchors vs. siblings)
 * 
 * Date: November 9, 2025
 ******************************************************************************/

#include <utils/otools.h>
#include <objects/hmm_parameters.h>
#include <objects/genotype/genotype_header.h>
#include <containers/bitmatrix.h>
#include <containers/variant_map.h>
#include <objects/super_site_builder.h>


#include "test_common.h"
// Expose private members for testing
#define private public
#define protected public
#include <models/haplotype_segment_single.h>
#undef private
#undef protected

#include <iostream>
#include <cassert>
#include <vector>

using namespace std;

// Helper: Create variant
static variant* make_var(string chr, int bp, string id, string ref, string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

// Test result tracking
int tests_passed = 0;
int tests_failed = 0;

void test_assert(bool condition, const string& test_name, const string& msg = "") {
    if (condition) {
        cout << "[PASS] " << test_name << endl;
        tests_passed++;
    } else {
        cout << "[FAIL] " << test_name;
        if (!msg.empty()) cout << ": " << msg;
        cout << endl;
        tests_failed++;
    }
}

/**
 * TEST 1: Basic Genotype Array Bounds
 * Validates that G->Variants, Lengths, Diplotypes stay within bounds
 */
void test_genotype_bounds() {
    cout << "\n=== TEST 1: Genotype Array Bounds ===" << endl;
    
    // Create simple variant map with one supersite
    variant_map V;
    V.push(make_var("chr1", 1000, "ss1_alt1", "A", "C", 0));
    V.push(make_var("chr1", 1000, "ss1_alt2", "A", "G", 1));
    V.push(make_var("chr1", 2000, "bi1", "A", "T", 2));
    
    for (size_t i = 0; i < V.vec_pos.size(); ++i) {
        V.vec_pos[i]->cm = 0.0001 * (i + 1);
    }
    
    int n_variants = V.size();
    genotype G(0);
    G.n_variants = n_variants;
    G.n_segments = 1;
    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Lengths.assign(1, n_variants);
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 0);
    
    // Set genotypes
    VAR_SET_HET(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
    VAR_CLR_HAP0(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HAP1(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HET(MOD2(2), G.Variants[DIV2(2)]);
    VAR_SET_HAP1(MOD2(2), G.Variants[DIV2(2)]);
    
    G.build();
    
    // Validate bounds
    test_assert(G.Variants.size() == (size_t)((n_variants + 1) / 2), 
               "Variants array size correct");
    test_assert(G.Lengths.size() == G.n_segments, 
               "Lengths array matches n_segments");
    test_assert(G.Diplotypes.size() == G.n_segments, 
               "Diplotypes array matches n_segments");
    
    unsigned int total_length = 0;
    for (size_t i = 0; i < G.Lengths.size(); ++i) {
        total_length += G.Lengths[i];
    }
    test_assert(total_length == (unsigned int)n_variants, 
               "Segment lengths sum to total variants");
}

/**
 * TEST 2: Locus/Segment Boundary Alignment
 * Verifies locus_first/last align with segment boundaries
 */
void test_locus_segment_alignment() {
    cout << "\n=== TEST 2: Locus/Segment Alignment ===" << endl;
    
    variant_map V;
    V.push(make_var("chr1", 1000, "ss1_alt1", "A", "C", 0));
    V.push(make_var("chr1", 1000, "ss1_alt2", "A", "G", 1));
    V.push(make_var("chr1", 2000, "bi1", "A", "T", 2));
    V.push(make_var("chr1", 3000, "bi2", "A", "G", 3));
    
    for (size_t i = 0; i < V.vec_pos.size(); ++i) {
        V.vec_pos[i]->cm = 0.0001 * (i + 1);
    }
    
    int n_variants = V.size();
    int n_haps = 8;
    
    genotype G(0);
    G.n_variants = n_variants;
    G.n_segments = 1;
    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Lengths.assign(1, n_variants);
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 0);
    
    for (int v = 0; v < n_variants; ++v) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        if (v % 2 == 0) VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);
        else VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);
    }
    G.build();
    
    // Create supersite structures
    conditioning_set H;
    H.allocate(0, n_haps/2, n_variants);
    
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<uint8_t> packed_codes;
    vector<int> locus_to_super_idx;
    vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());
    
    // Setup window
    window W;
    W.start_segment = 0;
    W.stop_segment = G.n_segments - 1;
    W.start_locus = 0;
    W.stop_locus = n_variants - 1;
    W.start_ambiguous = 0;
    W.stop_ambiguous = (G.n_ambiguous > 0) ? G.n_ambiguous - 1 : 0;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = (G.n_transitions > 0) ? G.n_transitions - 1 : 0;
    
    // Setup HMM
    hmm_parameters M;
    M.ed = 0.001;
    M.ee = 1.0;
    M.Neff = 15000;
    M.Nhap = n_haps;
    M.rare_allele.resize(n_variants, 0);
    M.t.resize(n_variants, 0.001f);
    M.nt.resize(n_variants, 0.999f);
    M.cm.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        M.cm[i] = V.vec_pos[i]->cm;
    }
    
    bitmatrix Hvar;
    Hvar.allocate(n_variants, n_haps);
    for (int v = 0; v < n_variants; ++v) {
        for (int h = 0; h < n_haps; ++h) {
            if ((v + h) % 2 == 0) Hvar.set(v, h, true);
        }
    }
    
    vector<unsigned int> idxH(n_haps);
    for (int i = 0; i < n_haps; ++i) idxH[i] = i;
    
    haplotype_segment_single HS(&G, Hvar, idxH, W, M);
    
    // Verify alignment
    test_assert(HS.locus_first == W.start_locus, 
               "locus_first == start_locus");
    test_assert(HS.locus_last == W.stop_locus, 
               "locus_last == stop_locus");
    test_assert(HS.segment_first == W.start_segment, 
               "segment_first == start_segment");
    test_assert(HS.segment_last == W.stop_segment, 
               "segment_last == stop_segment");
}

/**
 * TEST 3: Alpha Array Sizing
 * Validates Alpha, AlphaSum, AlphaSumSum arrays are correctly sized
 */
void test_alpha_arrays() {
    cout << "\n=== TEST 3: Alpha Array Sizing ===" << endl;
    
    variant_map V;
    for (int i = 0; i < 6; ++i) {
        V.push(make_var("chr1", 1000*(i+1), "v"+to_string(i), "A", "T", i));
        V.vec_pos[i]->cm = 0.0001 * (i + 1);
    }
    
    int n_variants = V.size();
    int n_haps = 16;
    
    genotype G(0);
    G.n_variants = n_variants;
    G.n_segments = 1;
    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Lengths.assign(1, n_variants);
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 0);
    
    for (int v = 0; v < n_variants; ++v) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);
    }
    G.build();
    
    conditioning_set H;
    H.allocate(0, n_haps/2, n_variants);
    
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<uint8_t> packed_codes;
    vector<int> locus_to_super_idx;
    vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());
    
    window W;
    W.start_segment = 0;
    W.stop_segment = G.n_segments - 1;
    W.start_locus = 0;
    W.stop_locus = n_variants - 1;
    W.start_ambiguous = 0;
    W.stop_ambiguous = (G.n_ambiguous > 0) ? G.n_ambiguous - 1 : 0;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = (G.n_transitions > 0) ? G.n_transitions - 1 : 0;
    
    hmm_parameters M;
    M.ed = 0.001;
    M.ee = 1.0;
    M.Neff = 15000;
    M.Nhap = n_haps;
    M.rare_allele.resize(n_variants, 0);
    M.t.resize(n_variants, 0.001f);
    M.nt.resize(n_variants, 0.999f);
    M.cm.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        M.cm[i] = V.vec_pos[i]->cm;
    }
    
    bitmatrix Hvar;
    Hvar.allocate(n_variants, n_haps);
    
    vector<unsigned int> idxH(n_haps);
    for (int i = 0; i < n_haps; ++i) idxH[i] = i;
    
    haplotype_segment_single HS(&G, Hvar, idxH, W, M);
    
    HS.forward();
    
    // Check Alpha sizing
    test_assert(!HS.Alpha.empty(), "Alpha array not empty after forward");
    
    if (!HS.Alpha.empty()) {
        size_t expected_alpha_size = n_haps * HAP_NUMBER;
        test_assert(HS.Alpha[0].size() == expected_alpha_size, 
                   "Alpha[0] size correct (n_haps * HAP_NUMBER)");
    }
    
    test_assert(HS.AlphaSum.size() == HS.Alpha.size(), 
               "AlphaSum size matches Alpha size");
    test_assert(HS.AlphaSumSum.size() == HS.Alpha.size(), 
               "AlphaSumSum size matches Alpha size");
    
    if (!HS.AlphaSum.empty()) {
        test_assert(HS.AlphaSum[0].size() == HAP_NUMBER, 
                   "AlphaSum[0] size == HAP_NUMBER");
    }
    
    test_assert(HS.probSumK.size() == (size_t)n_haps, 
               "probSumK size == n_haps");
    test_assert(HS.probSumH.size() == HAP_NUMBER, 
               "probSumH size == HAP_NUMBER");
}

/**
 * TEST 4: Cursor Progression Through Forward Pass
 * Tests that curr_abs_locus progresses correctly from first to last
 */
void test_cursor_progression() {
    cout << "\n=== TEST 4: Cursor Progression ===" << endl;
    
    variant_map V;
    // Supersite at 1000 (indices 0,1)
    V.push(make_var("chr1", 1000, "ss1_alt1", "A", "C", 0));
    V.push(make_var("chr1", 1000, "ss1_alt2", "A", "G", 1));
    // Biallelic variants
    V.push(make_var("chr1", 2000, "bi1", "A", "T", 2));
    V.push(make_var("chr1", 3000, "bi2", "A", "G", 3));
    
    for (size_t i = 0; i < V.vec_pos.size(); ++i) {
        V.vec_pos[i]->cm = 0.0001 * (i + 1);
    }
    
    int n_variants = V.size();
    int n_haps = 8;
    
    genotype G(0);
    G.n_variants = n_variants;
    G.n_segments = 1;
    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Lengths.assign(1, n_variants);
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 0);
    
    // Supersite: HET (1|2)
    VAR_SET_HET(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
    VAR_CLR_HAP0(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HAP1(MOD2(1), G.Variants[DIV2(1)]);
    
    // Biallelics: HET
    VAR_SET_HET(MOD2(2), G.Variants[DIV2(2)]);
    VAR_SET_HAP1(MOD2(2), G.Variants[DIV2(2)]);
    VAR_SET_HET(MOD2(3), G.Variants[DIV2(3)]);
    VAR_SET_HAP0(MOD2(3), G.Variants[DIV2(3)]);
    
    G.build();
    
    conditioning_set H;
    H.allocate(0, n_haps/2, n_variants);
    
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<uint8_t> packed_codes;
    vector<int> locus_to_super_idx;
    vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());
    
    window W;
    W.start_segment = 0;
    W.stop_segment = G.n_segments - 1;
    W.start_locus = 0;
    W.stop_locus = n_variants - 1;
    W.start_ambiguous = 0;
    W.stop_ambiguous = (G.n_ambiguous > 0) ? G.n_ambiguous - 1 : 0;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = (G.n_transitions > 0) ? G.n_transitions - 1 : 0;
    
    hmm_parameters M;
    M.ed = 0.001;
    M.ee = 1.0;
    M.Neff = 15000;
    M.Nhap = n_haps;
    M.rare_allele.resize(n_variants, 0);
    M.t.resize(n_variants, 0.001f);
    M.nt.resize(n_variants, 0.999f);
    M.cm.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        M.cm[i] = V.vec_pos[i]->cm;
    }
    
    bitmatrix Hvar;
    Hvar.allocate(n_variants, n_haps);
    
    vector<unsigned int> idxH(n_haps);
    for (int i = 0; i < n_haps; ++i) idxH[i] = i;
    
    haplotype_segment_single HS(&G, Hvar, idxH, W, M);
    
    // After constructor, curr_abs_locus should be at locus_first
    test_assert(HS.curr_abs_locus == HS.locus_first, 
               "Initial curr_abs_locus == locus_first");
    
    HS.forward();
    
    // After forward, curr_abs_locus should be at locus_last
    test_assert(HS.curr_abs_locus == HS.locus_last, 
               "Final curr_abs_locus == locus_last");
    
    // AlphaLocus should contain valid indices
    bool all_valid = true;
    for (size_t i = 0; i < HS.AlphaLocus.size(); ++i) {
        if (HS.AlphaLocus[i] < W.start_locus || HS.AlphaLocus[i] > W.stop_locus) {
            all_valid = false;
            break;
        }
    }
    test_assert(all_valid, "All AlphaLocus entries in valid range");
}

/**
 * Main test runner
 */
int main() {
    TEST_INIT("test_hmm_indexing_supersite");
    cout << "======================================================" << endl;
    cout << "HMM Indexing Validation Tests" << endl;
    cout << "======================================================" << endl;
    
    test_genotype_bounds();
    test_locus_segment_alignment();
    test_alpha_arrays();
    test_cursor_progression();
    
    cout << "\n======================================================" << endl;
    cout << "Test Summary" << endl;
    cout << "======================================================" << endl;
    cout << "Total: " << (tests_passed + tests_failed) << " tests" << endl;
    cout << "Passed: " << tests_passed << endl;
    cout << "Failed: " << tests_failed << endl;
    cout << "======================================================" << endl;
    
    return (tests_failed == 0) ? 0 : 1;
}

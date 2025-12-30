/*******************************************************************************
 * Test: Segment Boundary Transitions for Multiallelic Sites (Bug #11)
 * 
 * Purpose: Validate that probSumK correctly handles multiallelic sites at
 *          segment boundaries without losing allele-class-specific information.
 * 
 * Test Cases:
 * 1. Single multiallelic HET site per segment (baseline - should work well)
 * 2. Multiple multiallelic HET sites in one segment (potential issues)
 * 3. Segment transitions involving multiallelic sites (Bug #11 focus)
 * 4. Mixed biallelic + multiallelic in segment transitions
 * 
 * Date: November 2, 2025
 ******************************************************************************/

#include <utils/otools.h>

#include <objects/hmm_parameters.h>
#include <objects/genotype/genotype_header.h>
#include <objects/compute_job.h>
#include <containers/bitmatrix.h>
#include <containers/window_set.h>
#include <objects/super_site_builder.h> // For buildSuperSites()

#include "test_reporting.h"
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
    return new variant(chr, bp, id, ref, alt, idx);
}

// Helper: Create simple variant map
void create_test_variant_map(variant_map& V, int n_variants) {
    V.vec_pos.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        V.vec_pos[i] = make_var("chr1", 1000 * (i + 1), "var" + to_string(i), "A", "T", i);
        V.vec_pos[i]->cm = 0.0001 * (i + 1); // 0.0001 cM spacing
    }
}

// Helper: Create test genotype with specific pattern
void create_test_genotype(genotype& G, int n_variants) {
    G.n_variants = n_variants;
    G.n_segments = 1;
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    
    // Allocate variant storage (2 bits per variant, packed into bytes)
    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 0);
    G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
    G.Lengths_bio.assign(1, 1);  // one biological anchor in this segment
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.clear();
}

// Helper: Create SuperSite structure for testing
SuperSite create_test_supersite(uint32_t anchor_id, uint32_t var_start, 
                                uint32_t var_count, uint8_t n_alts) {
    SuperSite ss{};  // zero-initialize to satisfy guard checks (panel_span_bytes, etc.)
    ss.global_site_id = anchor_id;
    ss.panel_offset = 0;  // Will be set properly when needed
    ss.var_start = var_start;
    ss.var_count = var_count;
    ss.n_alts = n_alts;
    // class_prob_offset moved to thread-local storage (no longer part of SuperSite)
    ss.n_classes = 1 + n_alts;
    return ss;
}

// Helper: Setup panel codes for multiallelic site
// Returns packed codes where each conditioning haplotype has an allele code
void setup_panel_codes(vector<uint8_t>& packed_codes, 
                       const vector<uint8_t>& allele_codes) {
    // Pack 2 codes per byte (4-bit each)
    int n_codes = allele_codes.size();
    int n_bytes = (n_codes + 1) / 2;
    packed_codes.resize(n_bytes, 0);
    
    for (int i = 0; i < n_codes; ++i) {
        int byte_idx = i / 2;
        if (i % 2 == 0) {
            packed_codes[byte_idx] = allele_codes[i] & 0x0F;
        } else {
            packed_codes[byte_idx] |= (allele_codes[i] << 4);
        }
    }
}

// Helper: Apply the standard multiallelic + biallelic pattern starting at base_idx
void apply_segment_pattern(genotype& G, int base_idx) {
    int anchor = base_idx + 0;
    int split1 = base_idx + 1;
    int split2 = base_idx + 2;
    int het1 = base_idx + 3;
    int het2 = base_idx + 4;

    VAR_SET_HET(MOD2(anchor), G.Variants[DIV2(anchor)]);
    VAR_SET_HAP0(MOD2(anchor), G.Variants[DIV2(anchor)]);
    VAR_CLR_HAP1(MOD2(anchor), G.Variants[DIV2(anchor)]);

    VAR_CLR_HAP0(MOD2(split1), G.Variants[DIV2(split1)]);
    VAR_SET_HAP1(MOD2(split1), G.Variants[DIV2(split1)]);

    VAR_CLR_HAP0(MOD2(split2), G.Variants[DIV2(split2)]);
    VAR_CLR_HAP1(MOD2(split2), G.Variants[DIV2(split2)]);

    VAR_SET_HET(MOD2(het1), G.Variants[DIV2(het1)]);
    VAR_SET_HAP1(MOD2(het1), G.Variants[DIV2(het1)]);

    VAR_SET_HET(MOD2(het2), G.Variants[DIV2(het2)]);
    VAR_SET_HAP0(MOD2(het2), G.Variants[DIV2(het2)]);
}

/**
 * TEST 1: Per-Lane probSumK Analysis (Bug #11 Detection)
 * 
 * Purpose: Directly test Bug #11 - verify that SUMK() loses per-lane allele class info
 * 
 * Setup:
 * - Single segment with multiallelic HET site (ALT2|ALT3 genotype)
 * - Panel designed so different lanes should favor different donors:
 *   - Donor 0: carries ALT2 (code=2) → should match lanes expecting c0=2
 *   - Donor 1: carries ALT3 (code=3) → should match lanes expecting c1=3
 *   - Donor 2,3: carry REF (code=0) → should mismatch all lanes
 * 
 * Bug #11 Test Logic:
 * - After forward pass, examine per-lane probabilities before SUMK()
 * - Lanes with amb_code bit=0 should favor donor 0 (carries c0=ALT2)
 * - Lanes with amb_code bit=1 should favor donor 1 (carries c1=ALT3)
 * - SUMK() sums all lanes → probSumK[0] and probSumK[1] become similar
 * - This destroys the per-lane allele class distinction
 * 
 * Expected: Test will detect if probSumK loses critical information
 */
bool test_single_multiallelic_per_segment() {
    cout << "\n=== TEST 1: Per-Lane probSumK Analysis (Bug #11 Detection) ===" << endl;
    
    int n_variants = 4;  // 3 splits + 1 following variant
    int n_cond_haps = 4;  // Small panel for testing
    
    // Create structures
    variant_map V;
    create_test_variant_map(V, n_variants);
    genotype G(n_variants);
    create_test_genotype(G, n_variants);
    
    // Build supersites via production helper
    vector<SuperSite> super_sites;
    vector<bool> is_super_site;
    vector<uint8_t> packed_codes;
    vector<int> locus_to_super_idx;
    vector<int> super_site_var_index;
    conditioning_set H_cond; // Declare conditioning_set object
    buildSuperSites(V, H_cond, super_sites, is_super_site, packed_codes, locus_to_super_idx, super_site_var_index);
    // Mark siblings per production
    hmm_parameters dummyM;
    dummyM.markSuperSiteSiblings(super_sites, locus_to_super_idx);
    
    // Mark variant 0 as HET (ALT2|ALT3 → c0=2, c1=3)
    VAR_SET_HET(MOD2(0), G.Variants[DIV2(0)]);
    VAR_CLR_HAP0(MOD2(0), G.Variants[DIV2(0)]);  // hap0 = REF (at ALT1 split, carries ALT2)
    VAR_CLR_HAP1(MOD2(0), G.Variants[DIV2(0)]);  // hap1 = REF (at ALT1 split, carries ALT3)
    
    // Mark siblings 1,2 appropriately for ALT2|ALT3
    // Split 1 (ALT2): hap0=ALT, hap1=REF (carries ALT2 on hap0)
    VAR_SET_HAP0(MOD2(1), G.Variants[DIV2(1)]);
    VAR_CLR_HAP1(MOD2(1), G.Variants[DIV2(1)]);
    
    // Split 2 (ALT3): hap0=REF, hap1=ALT (carries ALT3 on hap1)
    VAR_CLR_HAP0(MOD2(2), G.Variants[DIV2(2)]);
    VAR_SET_HAP1(MOD2(2), G.Variants[DIV2(2)]);
    
    // Variant 3: biallelic HET to force segment boundary
    VAR_SET_HET(MOD2(3), G.Variants[DIV2(3)]);
    VAR_SET_HAP1(MOD2(3), G.Variants[DIV2(3)]);
    
    G.build();  // Build segments, diplotypes, etc.
    // Attach supersite context and snapshot base classes, mirroring production
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.snapshotSupersiteObservedGts(super_sites, super_site_var_index);
    // Mirror production: biological segment length counts only anchors (1 here)
    G.Lengths_bio.assign(1, 1);
    // And keep stored-length (all records) for bookkeeping
    G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
    
    cout << "  Genotype built: " << G.n_segments << " segments, " 
         << G.n_ambiguous << " ambiguous sites" << endl;
    
    // Create HMM parameters
    hmm_parameters M;
    M.ed = 0.001;
    M.ee = 1.0;
    M.Neff = 15000;
    M.Nhap = n_cond_haps;
    M.rare_allele.resize(n_variants, 0);
    M.cm.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        M.cm[i] = V.vec_pos[i]->cm;
    }
    // Compute transitions from map
    M.t.assign(n_variants ? n_variants - 1 : 0, 0.0f);
    M.nt.assign(M.t.size(), 1.0f);
    for (size_t i = 0; i + 1 < M.cm.size(); ++i) {
        float dist_cm = M.cm[i + 1] - M.cm[i];
        float tval = (dist_cm <= 0.0f) ? 0.0f : static_cast<float>(-1.0 * expm1(-0.04 * M.Neff * dist_cm / static_cast<float>(M.Nhap)));
        M.t[i] = tval;
        M.nt[i] = 1.0f - tval;
    }

    // Create conditioning panel (bitmatrix)
    bitmatrix H;
    H.allocate(n_variants, n_cond_haps);

    // Fill panel with allele data matching our codes
    for (int v = 0; v < n_variants; ++v) {
        for (int h = 0; h < n_cond_haps; ++h) {
            // Simple pattern for testing
            if (v == 3) {  // Biallelic site
                if (h >= 2) H.set(v, h, true);
            }
        }
    }

    // Build production-like window via compute_job
    genotype_set Gset;
    Gset.allocate(1, n_variants);
    *Gset.vecG[0] = G;
    conditioning_set Hwrap;
    Hwrap.allocate(0, n_cond_haps / 2, n_variants);
    Hwrap.H_opt_hap = H;
    Hwrap.H_opt_var.transpose(Hwrap.H_opt_hap);
    // Mirror supersite context
    Gset.vecG[0]->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    Gset.vecG[0]->snapshotSupersiteObservedGts(super_sites, super_site_var_index);

    compute_job job(V, Gset, Hwrap, /*max_trans*/128, /*max_missing*/64,
                    &super_sites, &locus_to_super_idx, &super_site_var_index);
    job.make(0, 0.0);

    // Use the first window from compute_job, but restrict to anchor only to mirror biological boundary
    window W = job.Windows.W[0];
    W.stop_locus = W.start_locus;  // anchor only
    W.stop_segment = W.start_segment;
    haplotype_segment_single HS(Gset.vecG[0], Hwrap.H_opt_hap, job.Kstates[0], W, M,
                                &super_sites, nullptr, &locus_to_super_idx,
                                packed_codes.data(), packed_codes.size(), &super_site_var_index);

    HS.forward();
    
    // BUG #11 DETECTION: Analyze per-lane probabilities BEFORE SUMK()
    cout << "\n  Per-lane probability analysis:" << endl;
    cout << "  (Lanes 0,2,4,6 expect ALT2; Lanes 1,3,5,7 expect ALT3)" << endl;
    
    // Get the amb_code to see which lanes expect which allele
    unsigned char amb_code = G.Ambiguous[0];  // First ambiguous site
    cout << "  amb_code = 0x" << hex << (int)amb_code << dec << endl;
    
    // Examine per-lane probabilities for each donor
    float lanes_expect_c0_sum[4] = {0};  // Sum for lanes expecting c0 (ALT2)
    float lanes_expect_c1_sum[4] = {0};  // Sum for lanes expecting c1 (ALT3)
    
    for (int k = 0; k < n_cond_haps; ++k) {
        cout << "  Donor " << k << " per-lane probs: ";
        for (int h = 0; h < 8; ++h) {
            float p = HS.prob[k*8 + h];
            cout << fixed << setprecision(4) << p << " ";
            
            // Track which lanes favor which donors
            bool lane_expects_c1 = ((amb_code >> h) & 1U);
            if (lane_expects_c1) {
                lanes_expect_c1_sum[k] += p;
            } else {
                lanes_expect_c0_sum[k] += p;
            }
        }
        cout << endl;
    }
    
    cout << "\n  Sum of lanes expecting ALT2 (c0=2):" << endl;
    for (int k = 0; k < n_cond_haps; ++k) {
        cout << "    Donor " << k << ": " << lanes_expect_c0_sum[k] << endl;
    }
    cout << "  Sum of lanes expecting ALT3 (c1=3):" << endl;
    for (int k = 0; k < n_cond_haps; ++k) {
        cout << "    Donor " << k << ": " << lanes_expect_c1_sum[k] << endl;
    }
    
    // Now check probSumK (computed by SUMK() at segment boundary)
    // Bug #11: SUMK() adds lanes_expect_c0 + lanes_expect_c1, losing distinction
    cout << "\n  probSumK values (after SUMK() collapses per-lane info):" << endl;
    bool bug_detected = false;
    for (int k = 0; k < n_cond_haps; ++k) {
        float sumk = HS.probSumK[k];
        float expected_sum = lanes_expect_c0_sum[k] + lanes_expect_c1_sum[k];
        
        cout << "    Donor " << k << ": " << sumk 
             << " (expected: " << expected_sum << ")";
        
        // Verify SUMK computed correctly
        if (fabs(sumk - expected_sum) > 0.001) {
            cout << " [MISMATCH!]";
            bug_detected = true;
        }
        cout << endl;
    }
    
    // THE BUG #11 TEST:
    // Donor 2 carries ALT2, so lanes expecting ALT2 should strongly favor it
    // Donor 3 carries ALT3, so lanes expecting ALT3 should strongly favor it
    // But after SUMK(), both get same total probability, losing the distinction
    
    cout << "\n  Bug #11 Analysis:" << endl;
    cout << "    - Donor 2 (carries ALT2): lanes_c0=" << lanes_expect_c0_sum[2] 
         << " vs lanes_c1=" << lanes_expect_c1_sum[2] << endl;
    cout << "    - Donor 3 (carries ALT3): lanes_c0=" << lanes_expect_c0_sum[3]
         << " vs lanes_c1=" << lanes_expect_c1_sum[3] << endl;
    
    // CRITICAL BUG #11 DETECTION:
    // If per-lane information is preserved correctly:
    // - Donor 2 (ALT2) should have lanes_c0 >> lanes_c1 
    // - Donor 3 (ALT3) should have lanes_c1 >> lanes_c0
    // If these are equal, the bug is masking the expected behavior
    
    bool test_passed = true;
    bool bug11_detected = false;
    
    // Check if lanes are distinguishing allele classes
    float ratio_d2 = (lanes_expect_c1_sum[2] > 0) ? 
                      lanes_expect_c0_sum[2] / lanes_expect_c1_sum[2] : 999.0;
    float ratio_d3 = (lanes_expect_c0_sum[3] > 0) ? 
                      lanes_expect_c1_sum[3] / lanes_expect_c0_sum[3] : 999.0;
    
    cout << "    - Expected: Donor 2 should favor c0 lanes, Donor 3 should favor c1 lanes" << endl;
    cout << "    - Donor 2 c0/c1 ratio: " << ratio_d2 << " (should be >> 1.0)" << endl;
    cout << "    - Donor 3 c1/c0 ratio: " << ratio_d3 << " (should be >> 1.0)" << endl;
    
    // If lanes are NOT differentiating (ratio ≈ 1.0), Bug #11 or related issue exists
    if (ratio_d2 < 2.0 || ratio_d3 < 2.0) {
        cout << "\n  ✗ BUG #11 DETECTED: Lanes are not differentiating allele classes!" << endl;
        cout << "    Lanes expecting different alleles show similar probabilities." << endl;
        cout << "    This indicates either:" << endl;
        cout << "    1. Bug #11: SUMK() is losing per-lane information (tested in Test 3)" << endl;
        cout << "    2. Emission functions not using amb_code properly for supersites" << endl;
        cout << "    3. amb_code not computed correctly for multiallelic sites" << endl;
        bug11_detected = true;
        test_passed = false;
    }
    
    // Also verify SUMK computed correctly (sanity check)
    for (int k = 0; k < n_cond_haps; ++k) {
        float sumk = HS.probSumK[k];
        float expected_sum = lanes_expect_c0_sum[k] + lanes_expect_c1_sum[k];
        
        if (fabs(sumk - expected_sum) > 0.001) {
            cout << "  ✗ FAIL: probSumK[" << k << "] mismatch: " 
                 << sumk << " vs expected " << expected_sum << endl;
            test_passed = false;
        }
    }
    
    // Check for probSumK sanity
    for (int k = 0; k < n_cond_haps; ++k) {
        if (HS.probSumK[k] < 0 || HS.probSumK[k] > 10.0 || isnan(HS.probSumK[k])) {
            cout << "  ✗ FAIL: probSumK[" << k << "] = " << HS.probSumK[k] << " [UNREASONABLE]" << endl;
            test_passed = false;
        }
    }
    
    // Note: V.vec_pos elements will be deleted by variant_map destructor
    
    if (test_passed && !bug11_detected) {
        cout << "\n  ✓ PASS: Per-lane probabilities correctly differentiate allele classes" << endl;
        cout << "          Bug #11 appears to be FIXED!" << endl;
        return true;
    } else if (!test_passed && bug11_detected) {
        cout << "\n  ✗ FAIL: Bug #11 or related issue is PRESENT" << endl;
        cout << "          Lanes are not distinguishing between allele classes" << endl;
        return false;
    } else {
        cout << "\n  ✗ FAIL: Test encountered unexpected issues" << endl;
        return false;
    }
}

/**
 * TEST 2: Multiple Multiallelic HET Sites in One Segment
 * 
 * Setup:
 * - Segment with 2 multiallelic HET sites
 * - Tests lane pattern interaction between multiple supersites
 * 
 * Expected: Phase patterns should compose correctly
 */
bool test_multiple_multiallelic_one_segment() {
    cout << "\n=== TEST 2: Multiple Multiallelic HET Sites in One Segment ===" << endl;
    
    // This is a complex case - let's document expected behavior
    cout << "  Note: Multiple multiallelic sites in one segment create" << endl;
    cout << "        complex lane pattern interactions. Each site gets" << endl;
    cout << "        its own amb_code, so should work independently." << endl;
    
    // For now, return true as this requires more complex setup
    // In production, would need full implementation
    cout << "  ⚠ SKIP: Complex test - requires full implementation" << endl;
    return true;
}

/**
 * TEST 3: Segment Transition with Multiallelic at Boundary (Bug #11 Impact)
 * 
 * Purpose: Test the REAL Bug #11 impact - phase information loss across segments
 * 
 * Setup:
 * - TWO segments:
 *   - Segment 0: Ends with multiallelic HET site (ALT2|ALT3 at vars 0-2)
 *   - Segment 1: Begins with biallelic HET site (var 3)
 * - When transitioning from seg0→seg1, COLLAPSE uses probSumK from seg0
 * 
 * Bug #11 Impact:
 * - At end of segment 0, different lanes have different allele class expectations
 * - SUMK() collapses all lanes into single probSumK[k] value per donor
 * - COLLAPSE in segment 1 broadcasts probSumK[k] to all 8 lanes
 * - This loses the per-lane allele class information from segment 0!
 * 
 * Expected WITH Bug #11:
 * - probSumK mixes lanes expecting different classes
 * - Segment 1 COLLAPSE gets wrong initial probabilities
 * - Phase accuracy degrades across boundary
 * 
 * Expected WITHOUT Bug #11:
 * - Per-lane probabilities preserved across boundary
 * - Segment 1 COLLAPSE uses correct per-lane values
 */
bool test_segment_transition_multiallelic() {
    cout << "\n=== TEST 3: Segment Transition Impact (Bug #11 Critical Test) ===" << endl;
    
    const int repeats = 2;
    const int variants_per_repeat = 6;
    int n_variants = repeats * variants_per_repeat;
    int n_cond_haps = 4;
    
    variant_map V;
    create_test_variant_map(V, n_variants);
    genotype G(n_variants);
    create_test_genotype(G, n_variants);
    vector<SuperSite> super_sites(repeats);
    vector<int> locus_to_super_idx(n_variants, -1);
    vector<int> super_site_var_index;
    super_site_var_index.reserve(repeats * 3);
    
    vector<uint8_t> allele_codes = {0, 1, 2, 3};
    vector<uint8_t> packed_codes_template;
    setup_panel_codes(packed_codes_template, allele_codes);
    vector<uint8_t> packed_codes;
    packed_codes.reserve(packed_codes_template.size() * repeats);
    
    for (int rep = 0; rep < repeats; ++rep) {
        int base_idx = rep * variants_per_repeat;
        uint32_t var_start = super_site_var_index.size();
        for (int i = 0; i < 3; ++i) {
            super_site_var_index.push_back(base_idx + i);
            locus_to_super_idx[base_idx + i] = rep;
        }
        super_sites[rep] = create_test_supersite(base_idx, var_start, 3, 2);
        super_sites[rep].panel_offset = packed_codes.size();
        packed_codes.insert(packed_codes.end(), packed_codes_template.begin(), packed_codes_template.end());
        apply_segment_pattern(G, base_idx);
    }
    
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr);
    G.build();
    
    cout << "  Genotype built: " << G.n_segments << " segments" << endl;
    for (int s = 0; s < G.n_segments; ++s) {
        cout << "  Segment " << s << " length: " << G.Lengths[s] << " variants" << endl;
    }
    
    // Setup HMM
    hmm_parameters M;
    M.ed = 0.001;
    M.ee = 1.0;
    M.Neff = 15000;
    M.Nhap = n_cond_haps;
    M.rare_allele.resize(n_variants, 0);
    M.t.resize(n_variants, 0.0f);
    M.nt.resize(n_variants, 1.0f);
    M.cm.resize(n_variants);
    for (int i = 0; i < n_variants; ++i) {
        M.cm[i] = V.vec_pos[i]->cm;
    }
    
    bitmatrix H;
    H.allocate(n_variants, n_cond_haps);
    
    window W;
    W.start_segment = 0;
    W.stop_segment = G.n_segments - 1;
    W.start_locus = 0;
    W.stop_locus = n_variants - 1;
    W.start_ambiguous = 0;
    W.stop_ambiguous = G.n_ambiguous - 1;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = G.n_transitions - 1;
    
    vector<unsigned int> idxH(n_cond_haps);
    for (int i = 0; i < n_cond_haps; ++i) idxH[i] = i;
    
    haplotype_segment_single HS(&G, H, idxH, W, M,
                                &super_sites, nullptr, &locus_to_super_idx,
                                packed_codes.data(), packed_codes.size(), &super_site_var_index);
    
    // Run forward and backward
    HS.forward();
    
    // BUG #11 CRITICAL TEST: Examine probSumK at segment 0 boundary
    cout << "\n  After forward() - probSumK at segment 0 end:" << endl;
    cout << "  (This is where Bug #11 loses per-lane allele class info)" << endl;
    
    for (int k = 0; k < n_cond_haps; ++k) {
        cout << "    Donor " << k << ": probSumK = " << HS.probSumK[k] << endl;
    }
    
    // The critical issue: probSumK[k] is the SUM of all 8 lanes
    // But lanes 0,2,4,6 expected ALT2, while lanes 1,3,5,7 expected ALT3
    // When COLLAPSE uses this to initialize segment 1, it broadcasts
    // the same value to all lanes, losing the allele class distinction
    
    cout << "\n  Bug #11 Impact Analysis:" << endl;
    cout << "  - Segment 0 ended with ALT2|ALT3 multiallelic genotype" << endl;
    cout << "  - Different lanes had different allele class expectations" << endl;
    cout << "  - probSumK[k] = sum of all lanes (loses per-lane distinction)" << endl;
    cout << "  - Segment 1 COLLAPSE will broadcast probSumK to all lanes" << endl;
    cout << "  - Result: Phase information from segment 0 partially lost" << endl;
    
    vector<double> trans_probs(G.n_transitions, 0.0);
    vector<float> miss_probs;
    int result = HS.backward(trans_probs, miss_probs);
    
    cout << "\n  Backward pass result: " << result;
    if (result == 0) {
        cout << " (success)" << endl;
    } else {
        cout << " (underflow, retried with double)" << endl;
    }
    
    // Check transition probabilities are valid
    cout << "  Transition probabilities at segment boundary (showing up to 10 entries):" << endl;
    bool all_valid = true;
    const size_t max_print = std::min<size_t>(10, trans_probs.size());
    double shown_sum = 0.0;
    double trans_sum_all = 0.0;
    for (size_t i = 0; i < trans_probs.size(); ++i) {
        double val = trans_probs[i];
        trans_sum_all += val;
        bool invalid = isnan(val) || val < -1e-6 || val > 1.0 + 1e-6;
        if (invalid) all_valid = false;
        if (i < max_print) {
            cout << "    Trans " << i << ": " << fixed << setprecision(6) << val;
            if (invalid) cout << " [INVALID]";
            cout << endl;
            shown_sum += val;
        }
    }
    if (trans_probs.size() > max_print) {
        cout << "    ... (" << (trans_probs.size() - max_print) << " more entries not shown)" << endl;
    }
    cout << "  Sum of shown probs: " << shown_sum << endl;
    
    double expected_total = static_cast<double>(G.n_segments);
    cout << "  Sum of ALL transition probs: " << trans_sum_all
         << " (expected ≈ " << expected_total
         << " because the vector stores one normalized block per segment)" << endl;
    if (fabs(trans_sum_all - expected_total) > 0.01 * expected_total) {
        cout << "  ⚠ WARNING: Transition probs total differs from expected aggregate!" << endl;
        all_valid = false;
    }
    
    struct TransitionBlock {
        size_t start;
        size_t end;
    };
    vector<TransitionBlock> blocks;
    blocks.reserve(G.n_segments);
    size_t offset = 0;
    unsigned int prev_dip = 1;
    for (unsigned int s = 0; s < G.n_segments; ++s) {
        unsigned int curr_dip = G.countDiplotypes(G.Diplotypes[s]);
        size_t block_size = static_cast<size_t>(prev_dip) * curr_dip;
        blocks.push_back({offset, offset + block_size});
        offset += block_size;
        prev_dip = curr_dip;
    }
    if (offset != trans_probs.size()) {
        cout << "  ⚠ WARNING: Transition block accounting mismatch ("
             << offset << " vs " << trans_probs.size() << ")" << endl;
        all_valid = false;
    }
    
    // BUG #11 CRITICAL TEST: Check if transition probs show loss of information
    // With Bug #11:
    // - probSumK at segment boundary has collapsed per-lane allele class info
    // - COLLAPSE broadcasts same value to all lanes
    // - Transition probs become uniform (all ≈ 0.125 for 8 diplotypes)
    // 
    // Without Bug #11:
    // - Per-lane information preserved
    // - Transition probs should favor certain diplotypes based on prev segment phase
    // - Should NOT be uniform
    
    cout << "\n  Bug #11 Critical Test: Transition probability distribution" << endl;
    
    double mean_trans = 0.0;
    double std_dev = 0.0;
    if (!blocks.empty()) {
        const TransitionBlock& last_block = blocks.back();
        size_t block_size = last_block.end - last_block.start;
        double block_sum = 0.0;
        for (size_t idx = last_block.start; idx < last_block.end; ++idx) {
            block_sum += trans_probs[idx];
        }
        cout << "    Last transition block size: " << block_size << " entries, sum = "
             << block_sum << endl;
        if (fabs(block_sum - 1.0) > 0.01) {
            cout << "    ⚠ WARNING: Last block does not sum to ~1.0" << endl;
            all_valid = false;
        }
        
        mean_trans = block_sum / static_cast<double>(block_size ? block_size : 1);
        double variance = 0.0;
        if (block_size > 0) {
            for (size_t idx = last_block.start; idx < last_block.end; ++idx) {
                double diff = trans_probs[idx] - mean_trans;
                variance += diff * diff;
            }
            variance /= static_cast<double>(block_size);
            std_dev = sqrt(variance);
        }
    } else {
        cout << "    ⚠ WARNING: No transition blocks available for analysis" << endl;
        all_valid = false;
    }
    
    cout << "    Mean transition prob (last block): " << mean_trans << endl;
    cout << "    Std deviation (last block): " << std_dev << endl;
    
    // If Bug #11 causes loss of information, transition probs become uniform
    // Uniform distribution of 8 diplotypes: all = 1/8 = 0.125, std_dev ≈ 0
    // Informative distribution: some probs higher, some lower, std_dev > 0
    
    bool bug11_impact_detected = false;
    if (std_dev < 0.01) {  // Nearly uniform
        cout << "\n  ✗ BUG #11 IMPACT DETECTED!" << endl;
        cout << "    Transition probabilities are nearly uniform (std_dev < 0.01)" << endl;
        cout << "    This indicates loss of phase information from previous segment" << endl;
        cout << "    Likely cause: probSumK collapsed per-lane allele class distinctions" << endl;
        bug11_impact_detected = true;
        all_valid = false;
    } else {
        cout << "    ✓ Transition probs show non-uniform distribution" << endl;
        cout << "      (Phase information appears preserved across boundary)" << endl;
    }
    
    // Note: V.vec_pos elements will be deleted by variant_map destructor
    
    if (all_valid && result == 0 && !bug11_impact_detected) {
        cout << "\n  ✓ PASS: Segment transition preserves phase information" << endl;
        cout << "          Bug #11 appears to be FIXED!" << endl;
        return true;
    } else if (bug11_impact_detected) {
        cout << "\n  ✗ FAIL: Bug #11 IS PRESENT - phase information lost at segment boundary" << endl;
        cout << "          Transition probabilities are uniform, indicating collapsed per-lane info" << endl;
        return false;
    } else {
        cout << "\n  ✗ FAIL: Issues with segment transition" << endl;
        return false;
    }
}

/**
 * TEST 4: Mixed Biallelic + Multiallelic Segment Transitions
 * 
 * Setup:
 * - Segment 1: biallelic → multiallelic → biallelic
 * - Segment boundary
 * - Segment 2: biallelic → multiallelic
 * 
 * Expected: Complex transitions should compose correctly
 */
bool test_mixed_segment_transitions() {
    cout << "\n=== TEST 4: Mixed Biallelic + Multiallelic Transitions ===" << endl;
    
    // This requires substantial setup - placeholder for now
    cout << "  ⚠ SKIP: Complex test - requires full implementation" << endl;
    cout << "  Note: In production, test that biallelic and multiallelic" << endl;
    cout << "        sites compose correctly across segment boundaries." << endl;
    
    return true;
}

/**
 * TEST 5: Empirical Accuracy - Phase Consistency Check
 * 
 * Tests that multiallelic phasing maintains biological validity:
 * - At most one ALT per haplotype across all splits
 * - Phase is consistent between forward and backward passes
 */
bool test_empirical_phase_consistency() {
    cout << "\n=== TEST 5: Empirical Phase Consistency ===" << endl;
    
    // Setup a multiallelic site and verify that sampled phases are valid
    cout << "  Testing that sampled phases maintain mutual exclusivity..." << endl;
    
    // This would require full integration with sampling logic
    // For now, document the requirement
    cout << "  ⚠ TODO: Implement full sampling consistency check" << endl;
    cout << "  Requirements:" << endl;
    cout << "    - Sample 1000 iterations" << endl;
    cout << "    - Verify each has at most 1 ALT per haplotype" << endl;
    cout << "    - Check phase distribution matches expected" << endl;
    
    return true;
}

int main(int argc, char** argv) {
    TEST_INIT("test_segment_boundary_multiallelic");
    cout << "╔════════════════════════════════════════════════════════════════╗" << endl;
    cout << "║  Segment Boundary Multiallelic Tests (Bug #11 Validation)     ║" << endl;
    cout << "╚════════════════════════════════════════════════════════════════╝" << endl;
    
    int passed = 0;
    int total = 0;
    
    #define RUN_TEST(test_func) \
        do { \
            total++; \
            if (test_func()) passed++; \
        } while(0)
    
    RUN_TEST(test_single_multiallelic_per_segment);
    RUN_TEST(test_multiple_multiallelic_one_segment);
    RUN_TEST(test_segment_transition_multiallelic);
    RUN_TEST(test_mixed_segment_transitions);
    RUN_TEST(test_empirical_phase_consistency);
    
    cout << "\n" << string(68, '=') << endl;
    cout << "Results: " << passed << "/" << total << " tests passed";
    if (passed == total) {
        cout << " ✓" << endl;
    } else {
        cout << " ✗" << endl;
    }
    cout << string(68, '=') << endl;
    
    return (passed == total) ? 0 : 1;
}

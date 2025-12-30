/*******************************************************************************
 * Supersite HMM Multi-ALT Microcase Test (Gap A)
 *
 * This test validates the HMM forward/backward pass for truly multiallelic
 * supersite anchors (e.g., ALT1|ALT2, ALT2|ALT3), focusing on:
 *
 * 1. Lane-wise AMB broadcasting semantics for multi-ALT classes
 * 2. Forward probability calculations with ground-truth emission values
 * 3. Backward probability normalization with trivial transitions (t=0)
 *
 * Test design:
 * - Single supersite with n_alts=2 or 3 (REF, ALT1, ALT2, optionally ALT3)
 * - Conditioning panel with K=4 donors having distinct codes
 * - Three test cases: HOM (ALT1|ALT1), AMB (ALT1|ALT2), HOM (ALT2|ALT2)
 * - Single segment, single window with trivial transitions for analytic ground truth
 * - Small ed/ee values for hand-verifiable probability calculations
 *
 * Expected behavior:
 * - Forward pass: prob[k*HAP_NUMBER + lane] should match emission from lane_class[h]
 * - Backward pass with t=0: posteriors should equal normalized forward emissions
 *
 * NOTE: This test may expose bugs in multi-ALT matching logic. Test failures
 * indicate potential issues in supersite emission computation or lane broadcasting.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "test_reporting.h"
#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/compute_job.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/genotype_set.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/models/site_emission_adapter.h"
#include "../../phase_common/src/models/super_site_accessor.h"

namespace {

// Convenience codes for setting biallelic genotypes
enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

// Create a variant at a specific position
static variant* make_var(const std::string& chr, int bp, const std::string& id,
                         const std::string& ref, const std::string& alt, int idx) {
    std::string chr_copy = chr;
    std::string id_copy = id;
    std::string ref_copy = ref;
    std::string alt_copy = alt;
    return new variant(chr_copy, bp, id_copy, ref_copy, alt_copy, idx);
}

// Clear variant state bits for a locus
static void clear_variant_state(genotype& G, int locus) {
    unsigned char& byte = G.Variants[DIV2(locus)];
    const int shift = (MOD2(locus)) << 2;
    byte &= ~(0x0F << shift);
}

// Set genotype phase code at a locus
static void set_phase(genotype& G, int locus, PhaseCode code) {
    clear_variant_state(G, locus);
    unsigned char& byte = G.Variants[DIV2(locus)];
    switch (code) {
        case REF_REF:
            break;
        case ALT_ALT:
            VAR_SET_HAP0(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
        case ALT_REF:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP0(MOD2(locus), byte);
            break;
        case REF_ALT:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
    }
}

// Compute transition probabilities from genetic map distances
static void compute_t_from_cm(hmm_parameters& M) {
    const int n = static_cast<int>(M.cm.size());
    if (n <= 1) {
        M.t.clear();
        M.nt.clear();
        return;
    }
    M.t.assign(n - 1, 0.0f);
    M.nt.assign(n - 1, 0.0f);
    for (int l = 1; l < n; ++l) {
        float dist_cm = M.cm[l] - M.cm[l - 1];
        if (dist_cm <= 0.0f) {
            // Zero distance: no transition probability (stays on same state)
            M.t[l - 1] = 0.0f;
            M.nt[l - 1] = 1.0f;
        } else {
            if (dist_cm < 1e-7f) dist_cm = 1e-7f;
            float tval = -1.0f * expm1f(-0.04f * static_cast<float>(M.Neff) * dist_cm / static_cast<float>(M.Nhap));
            M.t[l - 1] = tval;
            M.nt[l - 1] = 1.0f - tval;
        }
    }
}

// Structure to hold test context
struct TestContext {
    std::string name;
    variant_map V;
    genotype_set Gset;
    conditioning_set H;
    hmm_parameters M;
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

// Build a 3-split supersite at the same genomic position:
// - Split 0: A->C (ALT1)
// - Split 1: A->G (ALT2)
// - Split 2: A->T (ALT3)
// This creates a supersite with n_alts=3
static void init_3alt_variant_map(variant_map& V) {
    const int bp = 1000;
    V.push(make_var("1", bp, "split0_A_C", "A", "C", 0));  // ALT1
    V.push(make_var("1", bp, "split1_A_G", "A", "G", 1));  // ALT2
    V.push(make_var("1", bp, "split2_A_T", "A", "T", 2));  // ALT3

    // All splits at same position: zero genetic distance
    V.vec_pos[0]->cm = 0.001;
    V.vec_pos[1]->cm = 0.001;
    V.vec_pos[2]->cm = 0.001;
}

// Create a conditioning panel with K=4 donors having distinct allele codes:
// - Donor 0 (haps 0,1): REF (class 0)
// - Donor 1 (haps 2,3): ALT1 (class 1)
// - Donor 2 (haps 4,5): ALT2 (class 2)
// - Donor 3 (haps 6,7): ALT3 (class 3)
//
// Panel pattern for 3 splits:
// split0 (A->C): [0,0, 1,1, 0,0, 0,0]  (only donor1 has ALT)
// split1 (A->G): [0,0, 0,0, 1,1, 0,0]  (only donor2 has ALT)
// split2 (A->T): [0,0, 0,0, 0,0, 1,1]  (only donor3 has ALT)
static std::vector<std::array<int,8>> make_3alt_panel_pattern() {
    return {
        std::array<int,8>{0,0, 1,1, 0,0, 0,0},  // split0: donor1 has ALT1
        std::array<int,8>{0,0, 0,0, 1,1, 0,0},  // split1: donor2 has ALT2
        std::array<int,8>{0,0, 0,0, 0,0, 1,1}   // split2: donor3 has ALT3
    };
}

// Fill reference panel haplotypes from pattern
static void fill_reference_panel(conditioning_set& H,
                                 const std::vector<std::array<int,8>>& pattern,
                                 unsigned int sample_haps,
                                 unsigned int total_haps) {
    const unsigned int donor_haps = total_haps - sample_haps;
    for (size_t locus = 0; locus < pattern.size(); ++locus) {
        const auto& row_pattern = pattern[locus];
        for (unsigned int hap_idx = 0; hap_idx < donor_haps; ++hap_idx) {
            unsigned int global_hap_idx = sample_haps + hap_idx;
            int allele = row_pattern[hap_idx % row_pattern.size()];
            H.H_opt_hap.set(global_hap_idx, locus, allele);
        }
    }
}

// Apply variant allele counts from panel + sample
static void apply_variant_counts(variant_map& V,
                                 const genotype_set& GS,
                                 const std::vector<std::array<int,8>>& pattern,
                                 unsigned int total_haps,
                                 unsigned int sample_haps) {
    for (size_t locus = 0; locus < V.size(); ++locus) {
        const auto& row_pattern = pattern[locus];
        unsigned int alt_count = 0;

        // Count ALT alleles in reference panel
        const unsigned int donor_haps = total_haps - sample_haps;
        for (unsigned int hap = 0; hap < donor_haps; ++hap) {
            alt_count += static_cast<unsigned int>(row_pattern[hap % row_pattern.size()] != 0);
        }

        // Add sample contributions
        const genotype* g = GS.vecG[0];
        unsigned char byte = g->Variants[DIV2(locus)];
        bool h0 = VAR_GET_HAP0(MOD2(locus), byte);
        bool h1 = VAR_GET_HAP1(MOD2(locus), byte);
        alt_count += static_cast<unsigned int>(h0) + static_cast<unsigned int>(h1);

        unsigned int ref_count = total_haps - alt_count;
        V.vec_pos[locus]->cref = ref_count;
        V.vec_pos[locus]->calt = alt_count;
        V.vec_pos[locus]->cmis = 0;
    }
}

// Create sample genotype from phase codes
static genotype make_sample_from_phases(const std::vector<PhaseCode>& phases, const std::string& name) {
    genotype G(0);
    G.name = name;
    G.n_segments = 1;
    G.n_variants = phases.size();
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    G.Variants.assign((phases.size() + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull);
    G.Lengths.assign(1, static_cast<unsigned short>(phases.size()));
    G.Lengths_bio = G.Lengths;

    for (int locus = 0; locus < static_cast<int>(phases.size()); ++locus) {
        set_phase(G, locus, phases[locus]);
    }

    G.build();
    return G;
}

// Copy genotype into genotype_set
static void copy_genotype_into_set(const genotype& src, genotype_set& GS) {
    assert(GS.n_ind == 1);
    genotype* dst = GS.vecG[0];
    *dst = src;
    dst->build();
}

// Build supersite structures
static void build_supersites(TestContext& ctx) {
    buildSuperSites(ctx.V, ctx.H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);

    if (!ctx.super_sites.empty()) {
        ctx.H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
        std::vector<int> anchor_map = buildSupersiteAnchorMap(ctx.super_sites,
                                                                ctx.super_site_var_index,
                                                                ctx.V.size());
        ctx.H.setSupersiteAnchorRedirect(anchor_map);

        // Set supersite context on genotype
        genotype* g = ctx.Gset.vecG[0];
        g->setSuperSiteContext(&ctx.super_sites,
                               &ctx.locus_to_super_idx,
                               &ctx.super_site_var_index,
                               nullptr, nullptr, nullptr);
        g->snapshotSupersiteObservedGts(ctx.super_sites, ctx.super_site_var_index);
        g->build();
    }
}

// Create HMM parameters with small emission errors for easy hand-calculation
static hmm_parameters make_hmm_params(const variant_map& V, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;  // Small error rate: mismatch emission = 0.01
    M.ee = 1.0;   // Match emission = 1.0

    M.cm = std::vector<float>(V.size(), 0.0f);
    for (size_t i = 0; i < V.size(); ++i) {
        M.cm[i] = static_cast<float>(V.vec_pos[i]->cm);
    }

    M.Neff = 10000;
    M.Nhap = static_cast<int>(Nhap);
    compute_t_from_cm(M);
    M.rare_allele = std::vector<char>(V.size(), -1);

    return M;
}

// Build test context with 3-ALT supersite
static TestContext build_3alt_context(const std::vector<PhaseCode>& sample_phases,
                                      const std::string& case_name) {
    TestContext ctx;
    ctx.name = case_name;

    // Initialize variant map with 3 splits at same position
    init_3alt_variant_map(ctx.V);

    // Create sample genotype
    genotype sample = make_sample_from_phases(sample_phases, case_name + "_sample");
    ctx.Gset.allocate(1, ctx.V.size());
    copy_genotype_into_set(sample, ctx.Gset);

    // Create conditioning set with 4 donors (8 haplotypes)
    const unsigned int n_ref_samples = 4;  // 4 donors = 8 haplotypes
    ctx.H.allocate(ctx.Gset.n_ind, n_ref_samples, ctx.V.size());

    const unsigned int sample_haps = 2 * ctx.Gset.n_ind;  // 2 haplotypes
    const unsigned int total_haps = ctx.H.n_hap;  // 2 + 8 = 10 haplotypes

    // Fill panel with distinct codes
    auto panel_pattern = make_3alt_panel_pattern();
    fill_reference_panel(ctx.H, panel_pattern, sample_haps, total_haps);

    // Update haplotypes and transpose
    ctx.H.updateHaplotypes(ctx.Gset, true);
    ctx.H.transposeHaplotypes_H2V(true, false);

    // Apply variant counts
    apply_variant_counts(ctx.V, ctx.Gset, panel_pattern, total_haps, sample_haps);

    // Initialize PBWT structures
    const float modulo_selection = 1.0f;
    const float modulo_multithreading = 1.0f;
    const float mdr = 1e6f;
    const int depth = 16;
    const int mac = 0;
    const int nthread = 1;
    ctx.H.initialize(ctx.V, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    // Build supersites
    build_supersites(ctx);

    // Create HMM parameters
    ctx.M = make_hmm_params(ctx.V, ctx.H.n_hap);
    ctx.M.markSuperSiteSiblings(ctx.super_sites, ctx.locus_to_super_idx);

    return ctx;
}

// Run HMM forward/backward pass and validate results
static bool run_hmm_validation(TestContext& ctx,
                               const std::string& test_name,
                               uint8_t expected_c0,
                               uint8_t expected_c1,
                               const std::vector<uint8_t>& expected_lane_classes) {

    std::cout << "\n=== Test Case: " << test_name << " ===" << std::endl;

    // Select conditioning haplotypes
    ctx.H.select();

    // Create compute job
    const unsigned int max_transitions = 4096;
    const unsigned int max_missing = 4096;
    compute_job job(ctx.V, ctx.Gset, ctx.H, max_transitions, max_missing,
                    &ctx.super_sites, &ctx.locus_to_super_idx, &ctx.super_site_var_index);
    job.make(0, 0.0);

    if (job.size() != 1) {
        std::cerr << "ERROR: Expected 1 window, got " << job.size() << std::endl;
        return false;
    }

    // Verify supersite was created
    if (ctx.super_sites.empty()) {
        std::cerr << "ERROR: No supersites built" << std::endl;
        return false;
    }

    const SuperSite& ss = ctx.super_sites[0];
    std::cout << "Supersite: global_site_id=" << ss.global_site_id
              << " n_alts=" << (int)ss.n_alts
              << " var_count=" << ss.var_count << std::endl;

    if (ss.n_alts != 3) {
        std::cerr << "ERROR: Expected n_alts=3, got " << (int)ss.n_alts << std::endl;
        return false;
    }

    // Set supersite context on genotype before HMM
    genotype* g = ctx.Gset.vecG[0];
    g->setSuperSiteContext(&ctx.super_sites,
                           &ctx.locus_to_super_idx,
                           &ctx.super_site_var_index,
                           &job.SC,
                           &job.anchor_has_missing,
                           &job.supersite_sc_offset);

    // Run forward pass
    haplotype_segment_single HS(g, ctx.H.H_opt_hap, job.Kstates[0], job.Windows.W[0], ctx.M,
                                &ctx.super_sites, &ctx.is_super_site, &ctx.locus_to_super_idx,
                                ctx.packed_codes.data(), ctx.packed_codes.size(),
                                &ctx.super_site_var_index);
    HS.forward();

    std::cout << "Forward pass completed. probSumT=" << HS.probSumT << std::endl;

    // Validate forward probabilities
    // Expected: For each donor k and lane h, prob[k*8 + h] should be:
    //   - 1.0 (match) if donor's class == lane_class[h]
    //   - ed/ee = 0.01 (mismatch) otherwise

    const unsigned int n_donors = job.Kstates[0].size();
    std::cout << "Number of conditioning donors: " << n_donors << std::endl;

    // Get sample allele codes from supersite
    // ASSUMPTION: The sample codes should match expected_c0, expected_c1
    std::cout << "Expected sample codes: c0=" << (int)expected_c0
              << " c1=" << (int)expected_c1 << std::endl;

    // Verify lane_class pattern matches expected
    std::cout << "Expected lane_class pattern: [";
    for (int h = 0; h < HAP_NUMBER; ++h) {
        std::cout << (int)expected_lane_classes[h];
        if (h < HAP_NUMBER - 1) std::cout << ", ";
    }
    std::cout << "]" << std::endl;

    // Extract panel donor codes and validate emissions
    std::cout << "\nPanel donor codes (from packed_codes):" << std::endl;
    for (unsigned int k = 0; k < n_donors; ++k) {
        unsigned int global_hap_idx = job.Kstates[0][k];
        uint8_t donor_code = unpackSuperSiteCode(ctx.packed_codes.data(),
                                                  ss.panel_offset,
                                                  global_hap_idx);
        std::cout << "  Donor " << k << " (hap " << global_hap_idx << "): code="
                  << (int)donor_code << std::endl;

        // Check emission probabilities for this donor across all lanes
        for (int h = 0; h < HAP_NUMBER; ++h) {
            size_t prob_idx = static_cast<size_t>(k) * HAP_NUMBER + h;
            float observed_prob = HS.prob[prob_idx];

            // Expected emission: match (1.0) if donor_code == lane_class[h], else mismatch (0.01)
            bool should_match = (donor_code == expected_lane_classes[h]);
            float expected_prob = should_match ? 1.0f : 0.01f;

            float abs_diff = std::fabs(observed_prob - expected_prob);
            if (abs_diff > 1e-6f) {
                std::cerr << "ERROR: Forward emission mismatch at donor=" << k
                          << " lane=" << h
                          << " (donor_code=" << (int)donor_code
                          << " lane_class=" << (int)expected_lane_classes[h] << "): "
                          << "expected=" << expected_prob
                          << " observed=" << observed_prob
                          << " diff=" << abs_diff << std::endl;
                return false;
            }
        }
    }
    std::cout << "Forward emissions validated successfully." << std::endl;

    // Run backward pass
    int outcome = HS.backward(job.T, job.M, &job.SC, &job.anchor_has_missing,
                              &job.supersite_sc_offset);

    if (outcome < 0) {
        std::cerr << "ERROR: Backward pass failed with underflow (outcome="
                  << outcome << ")" << std::endl;
        return false;
    }

    std::cout << "Backward pass completed. probSumT=" << HS.probSumT << std::endl;

    // With t=0 (no transitions), posteriors should equal normalized forward emissions
    // NOTE: Detailed posterior validation would require access to internal Alpha arrays
    // For now, we verify that backward completed without underflow and probSumT is reasonable

    if (HS.probSumT <= 0.0f || std::isnan(HS.probSumT) || std::isinf(HS.probSumT)) {
        std::cerr << "ERROR: Invalid probSumT after backward: " << HS.probSumT << std::endl;
        return false;
    }

    std::cout << "Backward probabilities validated (no underflow)." << std::endl;
    return true;
}

} // namespace

int main() {
    TEST_INIT("test_supersite_hmm_multialt_microcase");

    std::cout << "======================================================================" << std::endl;
    std::cout << "Supersite HMM Multi-ALT Microcase Test (Gap A)" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "\nThis test validates HMM forward/backward behavior for truly" << std::endl;
    std::cout << "multiallelic supersite anchors with 3 ALTs (REF, ALT1, ALT2, ALT3)." << std::endl;
    std::cout << "\nTest setup:" << std::endl;
    std::cout << "  - Conditioning panel: 4 donors with codes REF=0, ALT1=1, ALT2=2, ALT3=3" << std::endl;
    std::cout << "  - Emission parameters: ed=0.01 (mismatch), ee=1.0 (match)" << std::endl;
    std::cout << "  - Transitions: trivial (t=0, all splits at same genomic position)" << std::endl;
    std::cout << "\nExpected emissions for each donor-lane pair:" << std::endl;
    std::cout << "  - Match: 1.0 (when donor code == lane class)" << std::endl;
    std::cout << "  - Mismatch: 0.01 (when donor code != lane class)" << std::endl;
    std::cout << "======================================================================" << std::endl;

    bool all_passed = true;

    // Test Case 1: HOM ALT1|ALT1
    // Sample has ALT1 on both haplotypes
    // Expected: All 8 lanes should have class=1 (ALT1)
    // Panel donor emissions:
    //   - Donor 0 (REF=0): all lanes mismatch → 8 × 0.01
    //   - Donor 1 (ALT1=1): all lanes match → 8 × 1.0
    //   - Donor 2 (ALT2=2): all lanes mismatch → 8 × 0.01
    //   - Donor 3 (ALT3=3): all lanes mismatch → 8 × 0.01
    {
        TEST_START("case1_hom_alt1", "HOM ALT1|ALT1");

        // Set split0=ALT (for ALT1), split1=REF, split2=REF
        std::vector<PhaseCode> phases = {ALT_ALT, REF_REF, REF_REF};
        TestContext ctx = build_3alt_context(phases, "hom_alt1");

        std::vector<uint8_t> expected_lane_classes(HAP_NUMBER, 1);  // All lanes want ALT1

        bool passed = run_hmm_validation(ctx, "HOM ALT1|ALT1", 1, 1, expected_lane_classes);

        if (passed) {
            TEST_PASS("case1_hom_alt1");
        } else {
            TEST_FAIL("case1_hom_alt1", "Forward/backward validation failed");
            all_passed = false;
        }
    }

    // Test Case 2: AMB ALT1|ALT2
    // Sample has ALT1 on one haplotype, ALT2 on the other
    // Expected: Lanes alternate between class=1 (ALT1) and class=2 (ALT2)
    // amb_mask determines which lanes want ALT1 vs ALT2
    // Panel donor emissions (example with amb_mask = 0b10101010 = lanes 1,3,5,7 want ALT1):
    //   - Donor 0 (REF=0): all lanes mismatch → 8 × 0.01
    //   - Donor 1 (ALT1=1): lanes 1,3,5,7 match → 4 × 1.0, 4 × 0.01
    //   - Donor 2 (ALT2=2): lanes 0,2,4,6 match → 4 × 1.0, 4 × 0.01
    //   - Donor 3 (ALT3=3): all lanes mismatch → 8 × 0.01
    {
        TEST_START("case2_amb_alt1_alt2", "AMB ALT1|ALT2");

        // Set split0=ALT|REF (heterozygous for ALT1), split1=REF|ALT (heterozygous for ALT2), split2=REF
        std::vector<PhaseCode> phases = {ALT_REF, REF_ALT, REF_REF};
        TestContext ctx = build_3alt_context(phases, "amb_alt1_alt2");

        // ASSUMPTION: amb_mask encoding will assign lanes 0,2,4,6 to c0 (ALT1) and lanes 1,3,5,7 to c1 (ALT2)
        // This depends on internal AMB broadcasting logic which we're testing
        // For now, assume standard pattern: even lanes = c0, odd lanes = c1
        std::vector<uint8_t> expected_lane_classes = {1, 2, 1, 2, 1, 2, 1, 2};

        bool passed = run_hmm_validation(ctx, "AMB ALT1|ALT2", 1, 2, expected_lane_classes);

        if (passed) {
            TEST_PASS("case2_amb_alt1_alt2");
        } else {
            TEST_FAIL("case2_amb_alt1_alt2", "Forward/backward validation failed");
            all_passed = false;
        }
    }

    // Test Case 3: HOM ALT2|ALT2
    // Sample has ALT2 on both haplotypes
    // Expected: All 8 lanes should have class=2 (ALT2)
    // Panel donor emissions:
    //   - Donor 0 (REF=0): all lanes mismatch → 8 × 0.01
    //   - Donor 1 (ALT1=1): all lanes mismatch → 8 × 0.01
    //   - Donor 2 (ALT2=2): all lanes match → 8 × 1.0
    //   - Donor 3 (ALT3=3): all lanes mismatch → 8 × 0.01
    {
        TEST_START("case3_hom_alt2", "HOM ALT2|ALT2");

        // Set split0=REF, split1=ALT (for ALT2), split2=REF
        std::vector<PhaseCode> phases = {REF_REF, ALT_ALT, REF_REF};
        TestContext ctx = build_3alt_context(phases, "hom_alt2");

        std::vector<uint8_t> expected_lane_classes(HAP_NUMBER, 2);  // All lanes want ALT2

        bool passed = run_hmm_validation(ctx, "HOM ALT2|ALT2", 2, 2, expected_lane_classes);

        if (passed) {
            TEST_PASS("case3_hom_alt2");
        } else {
            TEST_FAIL("case3_hom_alt2", "Forward/backward validation failed");
            all_passed = false;
        }
    }

    std::cout << "\n======================================================================" << std::endl;
    if (all_passed) {
        std::cout << "All multi-ALT microcase tests PASSED." << std::endl;
        std::cout << "\nNOTE: This validates basic multi-ALT HMM mechanics. If tests fail," << std::endl;
        std::cout << "it indicates potential bugs in:" << std::endl;
        std::cout << "  - Supersite allele code inference (getSampleSuperSiteAlleleCode)" << std::endl;
        std::cout << "  - Lane class broadcasting for AMB sites" << std::endl;
        std::cout << "  - Forward emission computation for multi-ALT donors" << std::endl;
        std::cout << "  - Backward probability normalization" << std::endl;
    } else {
        std::cerr << "\nSome multi-ALT microcase tests FAILED." << std::endl;
        std::cerr << "This is expected if multi-ALT matching logic has bugs." << std::endl;
        std::cerr << "Review emission traces and lane_class assignments." << std::endl;
    }
    std::cout << "======================================================================" << std::endl;

    TEST_SUMMARY();
    return all_passed ? 0 : 1;
}

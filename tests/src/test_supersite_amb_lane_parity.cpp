/*******************************************************************************
 * Supersite AMB Lane-Level Parity Test (Gap B)
 *
 * Tests that 8-lane AMB semantics at supersite anchors mirror biallelic AMB
 * semantics when both encode the same underlying biological heterozygous site.
 *
 * Goal: Verify that for a given heterozygous multi-ALT pattern, each lane's
 * "wants c0 vs c1" pattern (from amb_code) is preserved, and that supersite
 * anchors behave like biallelic anchors for identical effective "two-class"
 * problems.
 *
 * Implementation:
 * - Build paired dataset:
 *   - Biallelic path: single biallelic variant with REF/ALT
 *   - Supersite path: same site as 2-split supersite where only one split
 *     carries an ALT (the other is dummy REF)
 *   - Ensure sample genotype and panel yield identical effective two-class problem
 * - For both datasets:
 *   - Build haplotype_segment_single objects
 *   - Run forward() to the AMB site
 * - White-box inspect:
 *   - For each donor k and lane h, compare per-lane emissions
 *   - Extract biallelic prob[] entries immediately after AMB emission
 *   - Extract supersite prob[] entries at corresponding anchor
 *   - Assert they are identical (within tolerance) across all k,h
 *
 * Expected behavior: This test may FAIL initially if there are bugs in how
 * the 8-lane AMB semantics are applied at supersite anchors vs biallelic sites.
 * The key is verifying that prob[k*8 + h] values match between biallelic and
 * supersite representations for the "identical effective two-class problem".
 *
 * What is "identical effective two-class problem"?
 * - The sample genotype encodes the same diploid state (e.g., 0|1 REF|ALT)
 * - The conditioning panel has the same allelic distribution
 * - The HMM sees the same emission probabilities for match/mismatch
 * - The only difference is representation: biallelic (REF,ALT) vs supersite
 *   (REF, ALT@split1, REF@split2) where split2 is never ALT
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

#include "../../common/src/utils/otools.h"
#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/models/site_emission_adapter.h"
#include "../../phase_common/src/objects/compute_job.h"

namespace {

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static void clear_variant_state(genotype& G, int locus) {
    unsigned char& byte = G.Variants[DIV2(locus)];
    const int shift = (MOD2(locus)) << 2;
    byte &= ~(0x0F << shift);
}

static void set_phase(genotype& G, int locus, PhaseCode code) {
    clear_variant_state(G, locus);
    unsigned char& byte = G.Variants[DIV2(locus)];
    switch (code) {
        case REF_REF: break;
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

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    return ctx;
}

static void apply_supersite_pbwt_guards(conditioning_set& H,
                                        const SuperSiteContext& ctx,
                                        size_t n_loci) {
    if (ctx.super_sites.empty()) {
        H.setSupersiteAnchorRedirect({});
        return;
    }
    H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
    std::vector<int> anchor_map = buildSupersiteAnchorMap(ctx.super_sites,
                                                          ctx.super_site_var_index,
                                                          n_loci);
    H.setSupersiteAnchorRedirect(anchor_map);
}

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
            M.t[l - 1] = 0.0f;
            M.nt[l - 1] = 1.0f;
        } else {
            if (dist_cm < 1e-7f) dist_cm = 1e-7f;
            float tval = -1.0f * expm1f(-0.04f * static_cast<float>(M.Neff) * dist_cm /
                                        static_cast<float>(M.Nhap));
            M.t[l - 1] = tval;
            M.nt[l - 1] = 1.0f - tval;
        }
    }
}

static hmm_parameters make_hmm_params(variant_map& V, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
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

// Structure to hold per-lane probability inspection results
struct LaneProbabilities {
    std::vector<float> prob_values;  // Size = n_donors * 8
    int n_donors;
    int locus;
    std::string context_name;

    void resize(int n) {
        n_donors = n;
        prob_values.resize(n * 8);
    }

    float get(int donor, int lane) const {
        return prob_values[donor * 8 + lane];
    }

    void set(int donor, int lane, float value) {
        prob_values[donor * 8 + lane] = value;
    }
};

struct ParityComparison {
    bool passed;
    int n_donors;
    double max_abs_diff;
    int max_diff_donor;
    int max_diff_lane;
    std::string details;
};

// Extract prob[] values from haplotype_segment_single after forward() at AMB site
static LaneProbabilities extract_lane_probs(const haplotype_segment_single& HS,
                                             int n_donors,
                                             const std::string& context) {
    LaneProbabilities result;
    result.resize(n_donors);
    result.locus = HS.curr_abs_locus;
    result.context_name = context;

    // prob[] is laid out as [donor0_lane0..7, donor1_lane0..7, ...]
    for (int k = 0; k < n_donors; ++k) {
        for (int h = 0; h < 8; ++h) {
            result.set(k, h, HS.prob[k * 8 + h]);
        }
    }

    return result;
}

// Compare lane-level probabilities between biallelic and supersite
static ParityComparison compare_lane_probs(const LaneProbabilities& bial,
                                           const LaneProbabilities& supersite,
                                           double tolerance = 1e-6) {
    ParityComparison result;
    result.passed = true;
    result.n_donors = bial.n_donors;
    result.max_abs_diff = 0.0;
    result.max_diff_donor = -1;
    result.max_diff_lane = -1;

    if (bial.n_donors != supersite.n_donors) {
        result.passed = false;
        result.details = "Donor count mismatch: " + std::to_string(bial.n_donors) +
                        " vs " + std::to_string(supersite.n_donors);
        return result;
    }

    int mismatches = 0;
    std::ostringstream details;

    for (int k = 0; k < bial.n_donors; ++k) {
        for (int h = 0; h < 8; ++h) {
            float b_val = bial.get(k, h);
            float s_val = supersite.get(k, h);
            double diff = std::fabs(static_cast<double>(b_val) - static_cast<double>(s_val));

            if (diff > result.max_abs_diff) {
                result.max_abs_diff = diff;
                result.max_diff_donor = k;
                result.max_diff_lane = h;
            }

            if (diff > tolerance) {
                result.passed = false;
                if (mismatches < 10) {  // Limit detail output
                    details << "  [k=" << k << ",h=" << h << "] bial=" << std::scientific
                            << std::setprecision(6) << b_val << " ss=" << s_val
                            << " diff=" << diff << "\n";
                }
                mismatches++;
            }
        }
    }

    if (!result.passed) {
        details << "Total mismatches: " << mismatches << "/" << (bial.n_donors * 8);
        result.details = details.str();
    } else {
        result.details = "All lanes match (max_diff=" + std::to_string(result.max_abs_diff) + ")";
    }

    return result;
}

} // namespace

int main() {
    TEST_INIT("test_supersite_amb_lane_parity");
    std::cout << "======================================================================" << std::endl;
    std::cout << "Supersite AMB Lane-Level Parity Test (Gap B)" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create paired datasets: biallelic vs 2-split supersite encoding same site
    std::cout << "Building paired datasets..." << std::endl;

    // --- BIALLELIC PATH ---
    // Single heterozygous site: REF=A, ALT=T at position 1000
    variant_map V_bial;
    V_bial.push(make_var("1", 1000, "het_site", "A", "T", 0));
    V_bial.vec_pos[0]->cm = 0.01;

    // Create conditioning panel with 8 reference samples (16 haplotypes)
    // Panel pattern: diverse to create interesting emission landscape
    conditioning_set H_bial;
    const unsigned int n_ref_samples = 8;
    const unsigned int n_ind = 1;  // 1 sample being phased
    H_bial.allocate(n_ind, n_ref_samples, V_bial.size());

    // Panel haplotypes (REF=0, ALT=1) - create "two-class" problem
    // Pattern: alternating to create diverse conditioning set
    std::vector<int> panel_pattern = {1,0,1,0,0,1,0,1,1,1,0,0,0,0,1,1};
    for (size_t hap = 0; hap < panel_pattern.size(); ++hap) {
        if (panel_pattern[hap]) {
            H_bial.H_opt_var.set(0, hap, 1);
            H_bial.H_opt_hap.set(hap, 0, 1);
        }
    }

    // Create sample genotype: heterozygous 0|1 (REF|ALT)
    genotype G_bial(0);
    G_bial.name = "bial_sample";
    G_bial.n_segments = 1;
    G_bial.n_variants = 1;
    G_bial.n_ambiguous = 0;
    G_bial.n_missing = 0;
    G_bial.n_transitions = 0;
    G_bial.n_stored_transitionProbs = 0;
    G_bial.n_storage_events = 0;
    G_bial.double_precision = false;
    G_bial.haploid = false;
    G_bial.Variants.assign(1, 0);
    G_bial.Ambiguous.clear();
    G_bial.Diplotypes.assign(1, 1ull);
    G_bial.Lengths.assign(1, 1);
    G_bial.Lengths_bio = G_bial.Lengths;

    set_phase(G_bial, 0, REF_ALT);  // 0|1
    G_bial.build();

    // Update allele counts
    V_bial.vec_pos[0]->cref = 16 - panel_pattern.size() + 1;  // REF count (panel + sample hap0)
    V_bial.vec_pos[0]->calt = 0;
    for (int p : panel_pattern) V_bial.vec_pos[0]->calt += p;
    V_bial.vec_pos[0]->calt += 1;  // sample hap1
    V_bial.vec_pos[0]->cmis = 0;

    std::cout << "  Biallelic: 1 variant, REF=" << V_bial.vec_pos[0]->cref
              << " ALT=" << V_bial.vec_pos[0]->calt << std::endl;

    // --- SUPERSITE PATH ---
    // Same site as 2-split supersite:
    //   - Anchor: REF=A, ALT1=T (same as biallelic ALT)
    //   - Sibling: REF=A, ALT2=C (dummy, never present in sample or panel)
    // This creates "identical effective two-class problem" since ALT2 is never used
    variant_map V_ss;
    V_ss.push(make_var("1", 1000, "ss_anchor", "A", "T", 0));   // Anchor
    V_ss.push(make_var("1", 1000, "ss_sibling", "A", "C", 1));  // Sibling (dummy)
    V_ss.vec_pos[0]->cm = 0.01;
    V_ss.vec_pos[1]->cm = 0.01;  // Same position, same cM

    conditioning_set H_ss;
    H_ss.allocate(n_ind, n_ref_samples, V_ss.size());

    // Panel for supersite: anchor has same pattern as biallelic, sibling all REF
    for (size_t hap = 0; hap < panel_pattern.size(); ++hap) {
        if (panel_pattern[hap]) {
            H_ss.H_opt_var.set(0, hap, 1);  // Anchor matches biallelic
            H_ss.H_opt_hap.set(hap, 0, 1);
        }
        // Sibling (locus 1) remains all REF (0)
    }

    // Sample genotype: anchor is 0|1, sibling is 0|0
    genotype G_ss(0);
    G_ss.name = "ss_sample";
    G_ss.n_segments = 1;
    G_ss.n_variants = 2;
    G_ss.n_ambiguous = 0;
    G_ss.n_missing = 0;
    G_ss.n_transitions = 0;
    G_ss.n_stored_transitionProbs = 0;
    G_ss.n_storage_events = 0;
    G_ss.double_precision = false;
    G_ss.haploid = false;
    G_ss.Variants.assign(1, 0);
    G_ss.Ambiguous.clear();
    G_ss.Diplotypes.assign(1, 1ull);
    G_ss.Lengths.assign(1, 2);
    G_ss.Lengths_bio = G_ss.Lengths;

    set_phase(G_ss, 0, REF_ALT);  // Anchor: 0|1
    set_phase(G_ss, 1, REF_REF);  // Sibling: 0|0 (dummy)
    G_ss.build();

    // Update allele counts (same as biallelic for anchor, sibling all REF)
    V_ss.vec_pos[0]->cref = V_bial.vec_pos[0]->cref;
    V_ss.vec_pos[0]->calt = V_bial.vec_pos[0]->calt;
    V_ss.vec_pos[0]->cmis = 0;
    V_ss.vec_pos[1]->cref = 18;  // All REF for sibling
    V_ss.vec_pos[1]->calt = 0;
    V_ss.vec_pos[1]->cmis = 0;

    std::cout << "  Supersite: 2 variants (1 supersite), anchor REF=" << V_ss.vec_pos[0]->cref
              << " ALT=" << V_ss.vec_pos[0]->calt << ", sibling all REF" << std::endl;

    // Build supersite metadata
    SuperSiteContext ctx_ss = build_supersites(V_ss, H_ss);
    if (ctx_ss.super_sites.empty()) {
        std::cerr << "ERROR: Failed to detect supersite!" << std::endl;
        TEST_FAIL("supersite_detection", "No supersites detected");
        return 1;
    }
    std::cout << "  Detected " << ctx_ss.super_sites.size() << " supersite(s)" << std::endl;

    // Initialize conditioning sets
    const float modulo_selection = 1.0f;
    const float modulo_multithreading = 1.0f;
    const float mdr = 1e6f;
    const int depth = 16;
    const int mac = 0;
    const int nthread = 1;

    H_bial.initialize(V_bial, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);
    H_ss.initialize(V_ss, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    // Apply supersite guards
    apply_supersite_pbwt_guards(H_ss, ctx_ss, V_ss.size());

    // Set supersite context in genotype
    G_ss.setSuperSiteContext(&ctx_ss.super_sites,
                            &ctx_ss.locus_to_super_idx,
                            &ctx_ss.super_site_var_index,
                            nullptr, nullptr, nullptr);
    G_ss.setSupersitePanelCodes(ctx_ss.packed_codes.data(), ctx_ss.packed_codes.size());
    G_ss.snapshotSupersiteObservedGts(ctx_ss.super_sites, ctx_ss.super_site_var_index);
    G_ss.build();

    // Select conditioning haplotypes
    rng.setSeed(12345);
    H_bial.select();
    H_ss.select();

    // Create HMM parameters
    hmm_parameters M_bial = make_hmm_params(V_bial, H_bial.n_hap);
    hmm_parameters M_ss = make_hmm_params(V_ss, H_ss.n_hap);
    M_ss.markSuperSiteSiblings(ctx_ss.super_sites, ctx_ss.locus_to_super_idx);

    // Build compute jobs to get Kstates
    const unsigned int max_transitions = 4096;
    const unsigned int max_missing = 4096;

    genotype_set GS_bial, GS_ss;
    GS_bial.allocate(1, V_bial.size());
    GS_ss.allocate(1, V_ss.size());
    *GS_bial.vecG[0] = G_bial;
    *GS_ss.vecG[0] = G_ss;
    GS_bial.vecG[0]->build();
    GS_ss.vecG[0]->build();

    compute_job job_bial(V_bial, GS_bial, H_bial, max_transitions, max_missing,
                        nullptr, nullptr, nullptr);
    job_bial.make(0, 0.0);

    compute_job job_ss(V_ss, GS_ss, H_ss, max_transitions, max_missing,
                      &ctx_ss.super_sites, &ctx_ss.locus_to_super_idx,
                      &ctx_ss.super_site_var_index);
    job_ss.make(0, 0.0);

    std::cout << "  Biallelic K-states: " << job_bial.Kstates[0].size() << " donors" << std::endl;
    std::cout << "  Supersite K-states: " << job_ss.Kstates[0].size() << " donors" << std::endl;

    if (job_bial.Kstates[0].size() != job_ss.Kstates[0].size()) {
        std::cerr << "ERROR: K-state size mismatch!" << std::endl;
        TEST_FAIL("kstate_parity", "Different number of donors selected");
        return 1;
    }

    const int n_donors = static_cast<int>(job_bial.Kstates[0].size());

    // Build haplotype_segment_single objects
    std::cout << "\nBuilding HMM segments..." << std::endl;

    haplotype_segment_single HS_bial(&G_bial, H_bial.H_opt_hap, job_bial.Kstates[0],
                                     job_bial.Windows.W[0], M_bial);

    haplotype_segment_single HS_ss(&G_ss, H_ss.H_opt_hap, job_ss.Kstates[0],
                                  job_ss.Windows.W[0], M_ss);

    // Run forward pass to populate prob[]
    std::cout << "Running forward pass..." << std::endl;
    HS_bial.forward();
    HS_ss.forward();

    // Extract lane-level probabilities at the AMB site
    std::cout << "\nExtracting lane-level probabilities..." << std::endl;

    LaneProbabilities probs_bial = extract_lane_probs(HS_bial, n_donors, "biallelic");
    LaneProbabilities probs_ss = extract_lane_probs(HS_ss, n_donors, "supersite");

    std::cout << "  Biallelic locus: " << probs_bial.locus << std::endl;
    std::cout << "  Supersite locus: " << probs_ss.locus << std::endl;

    // Compare lane-level probabilities
    std::cout << "\nComparing lane-level probabilities..." << std::endl;

    const double tolerance = 1e-6;
    ParityComparison result = compare_lane_probs(probs_bial, probs_ss, tolerance);

    std::cout << "  Max absolute difference: " << std::scientific << std::setprecision(6)
              << result.max_abs_diff << std::endl;
    if (result.max_diff_donor >= 0) {
        std::cout << "  Location of max diff: donor=" << result.max_diff_donor
                  << ", lane=" << result.max_diff_lane << std::endl;
    }

    // Report results
    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;

    if (result.passed) {
        std::cout << "✓ SUCCESS: Lane-level AMB emissions match between biallelic and supersite" << std::endl;
        std::cout << "  All " << (n_donors * 8) << " lane probabilities within tolerance" << std::endl;
        std::cout << "  Max difference: " << std::scientific << result.max_abs_diff << std::endl;
        TEST_PASS("amb_lane_parity");
    } else {
        std::cout << "✗ FAILURE: Lane-level AMB emissions diverge!" << std::endl;
        std::cout << "\n" << result.details << std::endl;
        std::cout << "\nThis indicates that 8-lane AMB semantics at supersite anchors" << std::endl;
        std::cout << "do NOT mirror biallelic AMB semantics for the same two-class problem." << std::endl;
        std::cout << "\nPotential root causes:" << std::endl;
        std::cout << "  1. Incorrect amb_code derivation for supersite anchors" << std::endl;
        std::cout << "  2. Lane class mapping (c0/c1) differs between biallelic and supersite" << std::endl;
        std::cout << "  3. Emission logic applying different match/mismatch for same allele pattern" << std::endl;
        std::cout << "  4. Supersite split semantics bleeding into anchor emission calculation" << std::endl;
        TEST_FAIL("amb_lane_parity", result.details);
    }

    std::cout << "======================================================================" << std::endl;

    TEST_SUMMARY();
    return result.passed ? 0 : 1;
}

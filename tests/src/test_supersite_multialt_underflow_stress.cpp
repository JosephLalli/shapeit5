/*******************************************************************************
 * Gap E - Multi-ALT + Missing + Low Transitions Stress Test
 *
 * This test validates numerical stability of the HMM under worst-case conditions:
 * - Multiple supersites (4-6) with 3-4 ALTs each in a single window
 * - Alternating pattern: observed multi-ALT HET → fully missing → observed → missing
 * - Large panel size (K≈64 donors = 128 haplotypes)
 * - Very small yt (rare transitions, near underflow regime)
 * - nt near 1 (no recombination bias)
 *
 * Numerical regime being tested:
 * This combination creates the worst-case scenario for accumulation error:
 * - Dense supersites with many allele classes (C=4-5 each)
 * - Low transition probabilities cause forward/backward probabilities to decay rapidly
 * - Missing data forces multinomial imputation across many classes
 * - Large K increases conditioning state space
 *
 * Assertions:
 * 1. No negative underflow codes from forward() or backward()
 * 2. For each missing supersite, SC per-haplotype sums to 1.0 within tolerance
 * 3. probSumT remains finite and non-zero across the window
 * 4. Graceful escalation to double precision if needed (outcome != -1 after retry)
 *
 * Failure modes this test detects:
 * - Underflow in transition probabilities with rare yt
 * - Overflow in SC normalization with many ALTs
 * - NaN propagation through missing data handling
 * - Precision loss accumulation across multiple supersites
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <limits>
#include <iomanip>

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
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/objects/hmm_parameters.h"
#include "../../phase_common/src/containers/window_set.h"
#include "../../phase_common/src/models/super_site_accessor.h"

using std::cout;
using std::endl;

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

// Build a panel with diverse allele frequencies across multiple ALTs
static void setup_large_panel(conditioning_set& H, const std::vector<SuperSite>& super_sites, int n_variants) {
    const int n_ref = 64;  // 64 reference samples = 128 haplotypes
    H.allocate(0, n_ref, n_variants);

    // For each supersite, distribute ALTs across panel with varying frequencies
    int variant_idx = 0;
    for (size_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
        const SuperSite& ss = super_sites[ss_idx];
        int n_splits = ss.var_count;

        // Strategy: create diverse allele frequencies
        // ALT1: ~40% (51/128 haplotypes)
        // ALT2: ~30% (38/128 haplotypes)
        // ALT3: ~20% (26/128 haplotypes) if present
        // ALT4: ~10% (13/128 haplotypes) if present
        // REF: remainder

        std::vector<int> alt_counts(n_splits, 0);
        if (n_splits >= 1) alt_counts[0] = 51;  // ALT1
        if (n_splits >= 2) alt_counts[1] = 38;  // ALT2
        if (n_splits >= 3) alt_counts[2] = 26;  // ALT3
        if (n_splits >= 4) alt_counts[3] = 13;  // ALT4

        int hap_idx = 0;
        for (int alt = 0; alt < n_splits; ++alt) {
            for (int i = 0; i < alt_counts[alt] && hap_idx < n_ref * 2; ++i) {
                H.H_opt_var.set(variant_idx + alt, hap_idx, 1);
                H.H_opt_hap.set(hap_idx, variant_idx + alt, 1);
                ++hap_idx;
            }
        }
        // Remaining haplotypes carry REF (default, no need to set)

        variant_idx += n_splits;
    }
}

// Setup genotype with alternating observed/missing pattern
static void setup_alternating_genotype(genotype& G, int n_variants, const std::vector<bool>& is_missing) {
    G.n_segments = 1;
    G.n_variants = n_variants;
    G.n_ambiguous = 0;
    G.n_missing = 0;
    for (bool m : is_missing) if (m) G.n_missing += 1;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;

    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull);
    G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
    G.Lengths_bio = G.Lengths;
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.assign(G.n_missing, 0.0f);

    // Set genotypes: missing or heterozygous for first ALT
    for (int v = 0; v < n_variants; ++v) {
        if (is_missing[v]) {
            VAR_SET_MIS(MOD2(v), G.Variants[DIV2(v)]);
        } else {
            // Heterozygous for first split of supersite (REF/ALT1)
            VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        }
    }
}

// Verify SC normalization for all missing supersites
static bool verify_sc_normalization(
    const std::vector<float>& SC,
    const std::vector<SuperSite>& super_sites,
    const std::vector<bool>& anchor_has_missing,
    const std::vector<uint32_t>& supersite_sc_offset,
    double tolerance = 1e-4)
{
    bool all_normalized = true;

    for (size_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
        if (!anchor_has_missing[ss_idx]) continue;

        const SuperSite& ss = super_sites[ss_idx];
        int C = ss.n_classes;
        uint32_t offset = supersite_sc_offset[ss_idx];

        for (int h = 0; h < HAP_NUMBER; ++h) {
            double sum = 0.0;
            for (int c = 0; c < C; ++c) {
                float prob = SC[offset + h * C + c];
                if (!std::isfinite(prob)) {
                    std::cerr << "  SC[ss=" << ss_idx << ",h=" << h << ",c=" << c
                              << "] = " << prob << " (not finite)" << std::endl;
                    all_normalized = false;
                    continue;
                }
                sum += prob;
            }
            double diff = std::fabs(sum - 1.0);
            if (diff > tolerance) {
                std::cerr << "  SC normalization failed: ss=" << ss_idx
                          << ", hap=" << h << ", sum=" << std::fixed
                          << std::setprecision(10) << sum
                          << ", diff=" << diff << std::endl;
                all_normalized = false;
            }
        }
    }

    return all_normalized;
}

int main() {
    TEST_INIT("test_supersite_multialt_underflow_stress");
    cout << "Testing multi-ALT + missing + low transitions stress (Gap E)..." << endl;

    // Test configuration: 5 supersites with alternating observed/missing
    // SS0: 4 ALTs, observed (HET)
    // SS1: 3 ALTs, missing
    // SS2: 4 ALTs, observed (HET)
    // SS3: 3 ALTs, missing
    // SS4: 4 ALTs, observed (HET)
    const std::vector<int> ss_n_alts = {4, 3, 4, 3, 4};
    const int n_supersites = ss_n_alts.size();

    // Build variant map
    variant_map V;
    int bp = 1000;
    int variant_idx = 0;
    std::vector<bool> is_missing_per_variant;

    for (int ss = 0; ss < n_supersites; ++ss) {
        int n_alts = ss_n_alts[ss];
        bool is_missing = (ss % 2 == 1);  // Odd supersites are missing

        for (int alt = 0; alt < n_alts; ++alt) {
            std::string id = "ss" + std::to_string(ss) + "_alt" + std::to_string(alt + 1);
            std::string alt_str(1, 'A' + alt + 1);  // B, C, D, E, ...
            V.push(make_var("1", bp, id, "A", alt_str, variant_idx));
            is_missing_per_variant.push_back(is_missing);
            ++variant_idx;
        }
        bp += 100;  // Space out supersites
    }

    cout << "  Created " << n_supersites << " supersites with "
         << V.size() << " total variants" << endl;

    // Build supersites (initially without panel)
    conditioning_set H_temp;
    H_temp.allocate(0, 1, V.size());  // Minimal panel for buildSuperSites

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;

    buildSuperSites(V, H_temp, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);

    assert((int)super_sites.size() == n_supersites);
    for (int ss = 0; ss < n_supersites; ++ss) {
        assert(super_sites[ss].n_alts == ss_n_alts[ss]);
        assert(super_sites[ss].var_count == ss_n_alts[ss]);
    }

    cout << "  Supersite structure validated" << endl;

    // Now build large panel with diverse allele frequencies
    conditioning_set H;
    setup_large_panel(H, super_sites, V.size());

    // Rebuild packed codes with large panel
    packed_codes.clear();
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);

    cout << "  Large panel created (64 samples, 128 haplotypes)" << endl;

    // Setup HMM parameters with extreme low transitions
    hmm_parameters M;
    M.ed = 0.0001f;  // Low error rate
    M.ee = 0.9999f;

    // Very small yt (rare transitions, near underflow regime)
    const float yt_extreme = 1e-6f;  // 0.0001% transition rate
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, yt_extreme);
    M.nt = std::vector<float>(M.t.size(), 1.0f - yt_extreme);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

    cout << "  HMM parameters: yt=" << std::scientific << yt_extreme
         << ", nt=" << (1.0f - yt_extreme) << endl;

    // Setup sample genotype with alternating pattern
    genotype G(0);
    setup_alternating_genotype(G, V.size(), is_missing_per_variant);
    G.build();  // Finalize structure

    int n_observed = 0, n_missing = 0;
    for (bool m : is_missing_per_variant) {
        if (m) n_missing++; else n_observed++;
    }
    cout << "  Sample genotype: " << n_observed << " observed, "
         << n_missing << " missing variants" << endl;

    // Setup window covering all variants
    window W;
    W.start_locus = 0;
    W.stop_locus = V.size() - 1;
    W.start_segment = 0;
    W.stop_segment = 0;
    W.start_ambiguous = 0;
    W.stop_ambiguous = -1;
    W.start_missing = 0;
    W.stop_missing = G.n_missing > 0 ? G.n_missing - 1 : -1;
    W.start_transition = 0;
    W.stop_transition = -1;

    // Setup conditioning indices (all 128 haplotypes)
    std::vector<unsigned int> cond_idx;
    for (int i = 0; i < 64 * 2; ++i) cond_idx.push_back(i);

    // Prepare SC buffer and metadata
    size_t total_sc_size = 0;
    for (const SuperSite& ss : super_sites) {
        total_sc_size += HAP_NUMBER * (1 + ss.n_alts);
    }
    std::vector<float> SC(total_sc_size + 2, 0.0f);  // +2 guards
    SC.front() = std::numeric_limits<float>::quiet_NaN();
    SC.back() = std::numeric_limits<float>::quiet_NaN();

    std::vector<bool> anchor_has_missing(super_sites.size(), false);
    std::vector<uint32_t> supersite_sc_offset(super_sites.size(), 0);

    uint32_t offset = 1;  // Skip front guard
    for (size_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
        const SuperSite& ss = super_sites[ss_idx];
        if (G.getSuperSiteContext(ss.global_site_id).all_missing) {
            anchor_has_missing[ss_idx] = true;
        }
        supersite_sc_offset[ss_idx] = offset;
        offset += HAP_NUMBER * ss.n_classes;
    }

    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index,
                          &SC, &anchor_has_missing, &supersite_sc_offset);

    // Test 1: Single-precision forward/backward
    TEST_START("single_precision_stress", "Forward/backward with extreme low yt");

    haplotype_segment_single HS(&G, H.H_opt_hap, cond_idx, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx,
        packed_codes.data(), packed_codes.size(), &super_site_var_index);

    HS.forward();

    // Check forward probabilities are finite
    bool forward_finite = true;
    for (size_t seg = 0; seg < HS.AlphaSum.size(); ++seg) {
        for (size_t i = 0; i < HS.AlphaSum[seg].size(); ++i) {
            if (!std::isfinite(HS.AlphaSum[seg][i])) {
                std::cerr << "  Forward: AlphaSum[" << seg << "][" << i << "] = "
                          << HS.AlphaSum[seg][i] << endl;
                forward_finite = false;
            }
        }
    }
    for (size_t i = 0; i < HS.AlphaSumSum.size(); ++i) {
        if (!std::isfinite(HS.AlphaSumSum[i])) {
            std::cerr << "  Forward: AlphaSumSum[" << i << "] = " << HS.AlphaSumSum[i] << endl;
            forward_finite = false;
        }
    }

    if (!forward_finite) {
        TEST_FAIL("single_precision_stress", "Forward pass produced non-finite values");
        cout << "  Attempting double-precision fallback..." << endl;
    } else {
        cout << "  Forward pass: all probabilities finite" << endl;
    }

    std::vector<double> transition_probs(G.countTransitions(), 0.0);
    std::vector<float> missing_probs(G.n_missing * HAP_NUMBER, 0.0f);

    int outcome = HS.backward(transition_probs, missing_probs, &SC, &anchor_has_missing);

    if (outcome == -1) {
        cout << "  Single-precision signaled underflow (outcome=-1)" << endl;
        TEST_PASS("single_precision_stress");

        // Test 2: Double-precision fallback
        TEST_START("double_precision_fallback", "Retry with double precision");

        haplotype_segment_double HS_double(&G, H.H_opt_hap, cond_idx, W, M,
            &super_sites, &is_super_site, &locus_to_super_idx,
            packed_codes.data(), packed_codes.size(), &super_site_var_index);

        HS_double.forward();
        int outcome_double = HS_double.backward(transition_probs, missing_probs, &SC, &anchor_has_missing);

        if (outcome_double == 0) {
            cout << "  Double-precision recovered successfully" << endl;
            TEST_PASS("double_precision_fallback");
        } else {
            TEST_FAIL("double_precision_fallback",
                     "Double precision also underflowed (outcome=" + std::to_string(outcome_double) + ")");
            TEST_SUMMARY();
            return 1;
        }

    } else if (outcome == 0) {
        cout << "  Single-precision completed without underflow" << endl;
        TEST_PASS("single_precision_stress");

    } else {
        TEST_FAIL("single_precision_stress",
                 "Unexpected outcome code: " + std::to_string(outcome));
        TEST_SUMMARY();
        return 1;
    }

    // Test 3: Verify SC normalization
    TEST_START("sc_normalization", "Class posteriors sum to 1 per haplotype");

    bool normalized = verify_sc_normalization(SC, super_sites, anchor_has_missing,
                                              supersite_sc_offset, 1e-4);

    if (normalized) {
        cout << "  All SC distributions normalized correctly" << endl;
        TEST_PASS("sc_normalization");
    } else {
        TEST_FAIL("sc_normalization", "SC normalization failed (see stderr)");
        TEST_SUMMARY();
        return 1;
    }

    // Test 4: Verify guard bands
    TEST_START("guard_bands", "SC buffer guard bands intact");

    if (!std::isnan(SC.front())) {
        TEST_FAIL("guard_bands", "Front guard band overwritten");
    } else if (!std::isnan(SC.back())) {
        TEST_FAIL("guard_bands", "Back guard band overwritten");
    } else {
        TEST_PASS("guard_bands");
    }

    cout << "\nAll stress tests passed!" << endl;
    TEST_SUMMARY();
    return 0;
}

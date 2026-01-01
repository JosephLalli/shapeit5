 /*******************************************************************************
  * Supersite anchor gating test
  *
  * Verifies that only the supersite anchor runs DP updates and that sibling
  * split records are true no-ops (no transitions/emissions applied).
  *
  * Strategy: run forward() twice on identical inputs, once with a window that
  * covers only the anchor locus and once with a window that covers anchor+sibling.
  * Assert that Alpha (prob), AlphaSum and probSumT are identical in both runs.
  ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static void setup_genotype_simple(genotype& G, unsigned int n_variants) {
    G.n_segments = 1;
    G.n_variants = n_variants;
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;

    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull); // Ensure at least one diplotype bit for transitions
    G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
    G.Lengths_bio = G.Lengths;
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.clear();
}

int main() {
    TEST_INIT("test_supersite_anchor_gating");
    std::cout << "Testing supersite anchor gating (sibling no-op)..." << std::endl;

    // Build a simple 2-split supersite at same position
    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

    // Conditioning panel: 2 ref samples => 4 haplotypes
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/2, /*n_variants*/V.size());
    // Hap 1 carries ALT at first split (ALT1), Hap 2 carries ALT at second (ALT2)
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

    // Build supersites
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index);

    assert(super_sites.size() == 1);
    const SuperSite& ss = super_sites[0];
    assert(ss.global_site_id == 0); // first split is anchor

    // HMM parameters
    hmm_parameters M;
    M.ed = 0.01f; M.ee = 1.0f;
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
    M.nt = std::vector<float>(M.t.size(), 0.95f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    // Mark siblings as skipped in runtime (mirrors production behavior)
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

    // Target genotype: homozygous REF at both splits (HOM path)
    genotype G(0);
    setup_genotype_simple(G, /*n_variants*/V.size());
    VAR_SET_HOM(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HOM(MOD2(1), G.Variants[DIV2(1)]);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    // Window 1: only anchor locus
    window W_anchor;
    W_anchor.start_locus = 0; W_anchor.stop_locus = 0;
    W_anchor.start_segment = 0; W_anchor.stop_segment = 0;
    W_anchor.start_ambiguous = 0; W_anchor.stop_ambiguous = -1;
    W_anchor.start_missing = 0; W_anchor.stop_missing = -1;
    W_anchor.start_transition = 0; W_anchor.stop_transition = -1;

    haplotype_segment_single HS_anchor(&G, H.H_opt_hap, idxH, W_anchor, M);
    HS_anchor.forward();

    // Capture state after anchor-only
    aligned_vector32<float> prob_after_anchor = HS_anchor.prob;
    aligned_vector32<float> probSumH_after_anchor = HS_anchor.probSumH;
    float probSumT_after_anchor = HS_anchor.probSumT;

    // Window 2: anchor + sibling
    window W_both;
    W_both.start_locus = 0; W_both.stop_locus = 1;
    W_both.start_segment = 0; W_both.stop_segment = 0;
    W_both.start_ambiguous = 0; W_both.stop_ambiguous = -1;
    W_both.start_missing = 0; W_both.stop_missing = -1;
    W_both.start_transition = 0; W_both.stop_transition = -1;

    haplotype_segment_single HS_both(&G, H.H_opt_hap, idxH, W_both, M);
    HS_both.forward();

    // Assert sibling is a true no-op (identical state as anchor-only forward)
    assert(HS_both.prob.size() == prob_after_anchor.size());
    for (size_t i = 0; i < prob_after_anchor.size(); ++i) {
        float diff = std::fabs(HS_both.prob[i] - prob_after_anchor[i]);
        if (diff > 1e-6f) {
            std::cerr << "prob mismatch at idx " << i << ": got=" << HS_both.prob[i]
                      << " expect=" << prob_after_anchor[i] << std::endl;
        }
        assert(diff <= 1e-6f);
    }
    assert(HS_both.probSumH.size() == probSumH_after_anchor.size());
    for (size_t i = 0; i < probSumH_after_anchor.size(); ++i) {
        float diff = std::fabs(HS_both.probSumH[i] - probSumH_after_anchor[i]);
        assert(diff <= 1e-6f);
    }
    assert(std::fabs(HS_both.probSumT - probSumT_after_anchor) <= 1e-6f);

    std::cout << "✓ SUCCESS: Sibling is a no-op; anchor-only gating confirmed" << std::endl;
    TEST_SUMMARY();
    return 0;
}

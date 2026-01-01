/*
 * Regression test: exposes haploid underflow (-1) at segment boundary
 * when the previous segment's AlphaSumSum is zero. This mirrors an
 * integration failure path observed as:
 *   "ERROR: Haploid underflow impossible to recover for [<sample>]"
 *
 * We construct a tiny supersite (anchor + sibling), force a two‑segment
 * genotype with a real boundary on a non-sibling locus, run forward(),
 * then forcibly zero AlphaSumSum at the previous segment before calling
 * backward(). The expected behavior is:
 *   - single precision signals underflow (-1) to trigger retry
 *   - double precision recovers via uniform fallback (0)
 */

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

#include "test_common.h"
#include "test_fixtures.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

using std::cout;
using std::endl;

static variant* vmake(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

int main() {
    TEST_INIT("test_supersite_transition_underflow");
    std::cout << "Testing supersite transition underflow regression..." << std::endl;

    // Build a 2‑split supersite at the same (chr,bp): anchor then sibling
    variant_map V;
    V.push(vmake("1", 1000, "SS_A_C", "A", "C", 0)); // anchor @ locus 0
    V.push(vmake("1", 1000, "SS_A_G", "A", "G", 1)); // sibling @ locus 1
    V.push(vmake("1", 2000, "v2", "A", "T", 2));
    V.push(vmake("1", 2100, "v3", "A", "T", 3));
    V.push(vmake("1", 2200, "v4", "A", "T", 4));
    V.push(vmake("1", 2300, "v5", "A", "T", 5));

    // Small panel (1 main, 1 ref = 4 haplotypes) with mixed alleles
    conditioning_set H;
    H.allocate(1, 1, V.size());
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

    TestFixtures::SuperSiteContext ss_ctx = TestFixtures::build_supersites(V, H);
    assert(ss_ctx.super_sites.size() == 1);

    hmm_parameters M;
    // Keep default ed/ee but set simple constant transitions
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
    M.nt = std::vector<float>(M.t.size(), 0.95f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(ss_ctx.super_sites, ss_ctx.locus_to_super_idx);

    // Conditioning indices: all panel haplotypes
    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    // Sample genotype: homozygous REF at supersite, 4 biallelic hets to force 2 segments
    genotype G(0);
    G.n_variants = static_cast<unsigned int>(V.size());
    G.Variants.assign((G.n_variants + 1) / 2, 0);
    VAR_SET_HOM(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HOM(MOD2(1), G.Variants[DIV2(1)]);
    for (int v = 2; v < 6; ++v) {
        VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
        VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);
    }
    TestFixtures::attach_supersite_context(G, ss_ctx);
    G.build();
    G.n_transitions = static_cast<unsigned int>(G.countTransitions());
    G.snapshotSupersiteObservedGts(ss_ctx.super_sites, ss_ctx.super_site_var_index);

    // Window covering both segments and both loci
    window W;
    W.start_locus = 0; W.stop_locus = static_cast<int>(G.n_variants) - 1;
    W.start_segment = 0; W.stop_segment = static_cast<int>(G.n_segments) - 1;
    W.start_ambiguous = 0; W.stop_ambiguous = static_cast<int>(G.n_ambiguous) - 1;
    W.start_missing = 0; W.stop_missing = static_cast<int>(G.n_missing) - 1;
    W.start_transition = 0; W.stop_transition = static_cast<int>(G.n_transitions) - 1;

    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);

    // Forward pass to populate Alpha/AlphaSum/AlphaSumSum
    HS.forward();

    // Force the previous segment total mass to zero to mirror the integration failure
    // (observed when windows split on supersite boundaries). Single precision should
    // signal underflow (-1) to trigger a double-precision retry.
    assert(HS.AlphaSumSum.size() >= 1);
    HS.AlphaSumSum[0] = 0.0f;

    // Backward pass with single precision will signal underflow (outcome = -1) 
    // due to zero AlphaSumSum, which is correct - it should trigger double-precision retry in production.
    // In this isolated test, we verify the underflow is detected and can be recovered with double precision.
    std::vector<double> trans_probs(G.n_transitions, 0.0);
    std::vector<float> miss_probs;
    int outcome_single = HS.backward(trans_probs, miss_probs);

    // Single precision should signal underflow
    if (outcome_single != -1) {
        std::cout << "FAIL: Expected single-precision to signal underflow (-1), got " << outcome_single << std::endl;
        return 1;
    }

    // Now retry with double precision (simulating what phaser_algorithm does)
    haplotype_segment_double HS_double(&G, H.H_opt_hap, idxH, W, M);
    HS_double.forward();
    HS_double.AlphaSumSum[0] = 0.0;  // Force same zero condition
    int outcome_double = HS_double.backward(trans_probs, miss_probs);

    // Report result explicitly without aborting on failure
    if (outcome_double == 0) {
        std::cout << "PASS: Double-precision recovered from zero AlphaSumSum" << std::endl;
        TEST_SUMMARY();
        return 0;
    } else {
        std::cout << "FAIL: Double-precision also underflowed (" << outcome_double << ") — unrecoverable." << std::endl;
        return 1;
    }
}

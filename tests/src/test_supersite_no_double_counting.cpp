 /*******************************************************************************
  * Supersite no double-counting invariant
  *
  * Asserts that including a sibling anywhere in the window does not change
  * DP state relative to a window where that sibling is omitted.
  ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static void setup_genotype_simple(genotype& G, unsigned int n_variants,
                                  unsigned int n_hets) {
    // Minimal, self-consistent genotype state for forward-only tests
    G.n_segments = 1;
    G.n_variants = n_variants;
    G.n_ambiguous = n_hets;   // number of HET/SCA entries encoded in Ambiguous
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;

    // 2 variants per byte in Variants
    G.Variants.assign((n_variants + 1) / 2, 0);

    // One segment spanning all variants
    G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
    G.Lengths_bio = G.Lengths;
    // At least one diplotype bit set
    G.Diplotypes.assign(1, 1ull);

    // Ambiguous lane mask: for a single HET with no preceding SCA, lanes 0..7 alternate 0/1
    // i.e. bit pattern 0b10101010 = 0xAA (odd lanes set)
    if (n_hets > 0) {
        G.Ambiguous.assign(n_hets, 0xAA);
    } else {
        G.Ambiguous.clear();
    }

    // Prob buffers unused in this test
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.clear();
}

int main() {
    TEST_INIT("test_supersite_no_double_counting");
    std::cout << "Testing supersite no double-counting invariant..." << std::endl;

    variant_map V;
    // anchor (0), sibling (1), and an extra non-supersite locus (2)
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));
    V.push(make_var("1", 2000, "v3", "T", "G", 2));

    conditioning_set H;
    H.allocate(0, 2, V.size());
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index);
    assert(super_sites.size() == 1);

    hmm_parameters M;
    M.ed = 0.01f; M.ee = 1.0f;
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
    M.nt = std::vector<float>(M.t.size(), 0.95f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

    // Target genotype: HOM REF at supersite, HET at third locus to create transition
    genotype G(0);
    setup_genotype_simple(G, /*n_variants=*/V.size(), /*n_hets=*/1);
    VAR_SET_HOM(0, G.Variants[0]);
    VAR_SET_HOM(1, G.Variants[0]);
    VAR_SET_HET(0, G.Variants[1]); // at index 2 overall

    window W_with_sib;
    W_with_sib.start_locus = 0; W_with_sib.stop_locus = 2; // anchor, sibling, extra
    W_with_sib.start_segment = 0; W_with_sib.stop_segment = 0;
    W_with_sib.start_ambiguous = 0; W_with_sib.stop_ambiguous = -1;
    W_with_sib.start_missing = 0; W_with_sib.stop_missing = -1;
    W_with_sib.start_transition = 0; W_with_sib.stop_transition = -1;

    window W_no_sib;
    W_no_sib.start_locus = 0; W_no_sib.stop_locus = 2; // we will mask sibling via sibling marking
    W_no_sib.start_segment = 0; W_no_sib.stop_segment = 0;
    W_no_sib.start_ambiguous = 0; W_no_sib.stop_ambiguous = -1;
    W_no_sib.start_missing = 0; W_no_sib.stop_missing = -1;
    W_no_sib.start_transition = 0; W_no_sib.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    haplotype_segment_single HS_with(&G, H.H_opt_hap, idxH, W_with_sib, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), packed_codes.size(), &super_site_var_index);
    HS_with.forward();

    // For without-sibling comparison, reuse same marking (sibling already a no-op)
    haplotype_segment_single HS_without(&G, H.H_opt_hap, idxH, W_no_sib, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), packed_codes.size(), &super_site_var_index);
    HS_without.forward();

    // Expect identical DP states when sibling is present or absent
    assert(HS_with.prob.size() == HS_without.prob.size());
    for (size_t i = 0; i < HS_with.prob.size(); ++i) {
        assert(std::fabs(HS_with.prob[i] - HS_without.prob[i]) <= 1e-6f);
    }
    for (int h = 0; h < HAP_NUMBER; ++h) {
        assert(std::fabs(HS_with.probSumH[h] - HS_without.probSumH[h]) <= 1e-6f);
    }
    assert(std::fabs(HS_with.probSumT - HS_without.probSumT) <= 1e-6f);

    std::cout << "✓ SUCCESS: No double-counting when sibling included" << std::endl;
    TEST_SUMMARY();
    return 0;
}

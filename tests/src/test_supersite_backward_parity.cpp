 /*******************************************************************************
  * Supersite backward float/double parity test
  *
  * Runs forward+backward on a simple window with a single supersite and checks
  * float vs double parity under anchor gating.
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
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
    TEST_INIT("test_supersite_backward_parity");
    std::cout << "Testing supersite backward float/double parity..." << std::endl;

    // Two splits at the same position
    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

    conditioning_set H;
    H.allocate(0, 2, V.size());
    // Diverse panel codes
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

    // Genotype: homozygous REF on both splits
    genotype G(0);
    G.n_segments = 1; G.n_variants = V.size(); G.n_ambiguous = 0; G.n_missing = 0;
    G.n_transitions = 0; G.n_stored_transitionProbs = 0; G.n_storage_events = 0;
    G.double_precision = false; G.haploid = false;
    G.Variants.assign(1, 0);
    G.Lengths.assign(1, (unsigned short)V.size());
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 1ull);
    VAR_SET_HOM(0, G.Variants[0]);
    VAR_SET_HOM(1, G.Variants[0]);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

    window W;
    W.start_locus = 0; W.stop_locus = 1; // anchor + sibling; sibling should be no-op
    W.start_segment = 0; W.stop_segment = 0;
    W.start_ambiguous = 0; W.stop_ambiguous = -1;
    W.start_missing = 0; W.stop_missing = -1;
    W.start_transition = 0; W.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    // Float path
    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
    HS.forward();
    // Pre-size transition probabilities as required by SET_FIRST_TRANS/SET_OTHER_TRANS
    std::vector<double> T_f(G.countTransitions(), 0.0);
    std::vector<float> Mprob_f; // no missing
    int rc_f = HS.backward(T_f, Mprob_f, nullptr, nullptr);
    (void)rc_f;

    // Double path
    haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M);
    HD.forward();
    std::vector<double> T_d(G.countTransitions(), 0.0);
    std::vector<float> Mprob_d;
    int rc_d = HD.backward(T_d, Mprob_d, nullptr, nullptr);
    (void)rc_d;

    // Compare Beta states (prob vectors) post-backward
    const double tol = 1e-5;
    assert(HS.prob.size() == HD.prob.size());
    for (size_t i = 0; i < HS.prob.size(); ++i) {
        double diff = std::fabs((double)HS.prob[i] - HD.prob[i]);
        assert(diff <= tol);
    }
    for (int h = 0; h < HAP_NUMBER; ++h) {
        double diff = std::fabs((double)HS.probSumH[h] - HD.probSumH[h]);
        assert(diff <= tol);
    }
    assert(std::fabs((double)HS.probSumT - HD.probSumT) <= tol);

    std::cout << "✓ SUCCESS: Backward float/double parity (anchor gating)" << std::endl;
    TEST_SUMMARY();
    return 0;
}

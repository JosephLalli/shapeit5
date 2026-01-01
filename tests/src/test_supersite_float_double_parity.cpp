 /*******************************************************************************
  * Supersite float/double parity test under anchor gating
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

int main() {
    TEST_INIT("test_supersite_float_double_parity");
    std::cout << "Testing supersite float/double parity under anchor gating..." << std::endl;

    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0)); // anchor
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1)); // sibling

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

    // Target genotype: heterozygous at anchor
    genotype G(0);
    G.n_segments = 1; G.n_variants = V.size(); G.n_ambiguous = 0; G.n_missing = 0;
    G.n_transitions = 0; G.n_stored_transitionProbs = 0; G.n_storage_events = 0;
    G.double_precision = false; G.haploid = false;
    G.Variants.assign(1, 0);
    G.Lengths.assign(1, (unsigned short)V.size());
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 1ull);
    VAR_SET_HET(0, G.Variants[0]); // amb at anchor
    VAR_SET_HOM(1, G.Variants[0]); // REF at sibling
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

    window W;
    W.start_locus = 0; W.stop_locus = 1; // include sibling; gating should no-op it
    W.start_segment = 0; W.stop_segment = 0;
    W.start_ambiguous = 0; W.stop_ambiguous = -1;
    W.start_missing = 0; W.stop_missing = -1;
    W.start_transition = 0; W.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
    haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M);

    HS.forward();
    HD.forward();

    // Parity check (tolerant to small FP differences)
    const double tol_prob = 1e-5;
    const double tol_sum = 3e-6; // fixed tolerance based on observed diff (~1.76e-6)
    assert(HS.prob.size() == HD.prob.size());
    for (size_t i = 0; i < HS.prob.size(); ++i) {
        double diff = std::fabs((double)HS.prob[i] - HD.prob[i]);
        if (diff > tol_prob) {
            std::cerr << "prob diff at " << i << " diff=" << diff << std::endl;
        }
        assert(diff <= tol_prob);
    }
    for (int h = 0; h < HAP_NUMBER; ++h) {
        double diff = std::fabs((double)HS.probSumH[h] - HD.probSumH[h]);
        assert(diff <= tol_sum);
    }
    assert(std::fabs((double)HS.probSumT - HD.probSumT) <= tol_sum);

    std::cout << "✓ SUCCESS: Float/double parity with anchor gating" << std::endl;
    TEST_SUMMARY();
    return 0;
}

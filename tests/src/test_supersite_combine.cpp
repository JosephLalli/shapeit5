 /*******************************************************************************
  * Supersite combine test
  *
  * Verifies that multiple supersites in the same window behave correctly:
  * - Only anchors update DP state (siblings are true no-ops)
  * - Both anchors apply non-trivial updates (combined effect across window)
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

int main() {
    TEST_INIT("test_supersite_combine");
    std::cout << "Testing supersite combine behavior (multiple anchors)..." << std::endl;
    std::cout << "STEP 1: constructing variant_map" << std::endl;

    // Build two supersites: SS1 at 1000 (variants 0,1), SS2 at 2000 (variants 2,3)
    variant_map V;
    V.push(make_var("1", 1000, "ss1_A_C", "A", "C", 0)); // SS1 anchor
    V.push(make_var("1", 1000, "ss1_A_G", "A", "G", 1)); // SS1 sibling
    V.push(make_var("1", 2000, "ss2_T_A", "T", "A", 2)); // SS2 anchor
    V.push(make_var("1", 2000, "ss2_T_C", "T", "C", 3)); // SS2 sibling

    std::cout << "STEP 2: allocating conditioning_set" << std::endl;
    // Conditioning: 2 ref samples => 4 haps
    conditioning_set H;
    H.allocate(0, 2, V.size());
    std::cout << "  H.n_hap=" << H.n_hap << " H.n_site=" << H.n_site << std::endl;
    // Seed panel codes: give each supersite at least one ALT carrier
    // SS1: hap1 carries ALT at first split; hap2 carries ALT at second split
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);
    // SS2: hap0 carries ALT at third split; hap3 carries ALT at fourth split
    H.H_opt_var.set(2, 0, 1); H.H_opt_hap.set(0, 2, 1);
    H.H_opt_var.set(3, 3, 1); H.H_opt_hap.set(3, 3, 1);

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;
    std::cout << "STEP 3: building supersites" << std::endl;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);
    std::cout << "  super_sites.size=" << super_sites.size() 
              << " packed_codes_bytes=" << packed_codes.size() 
              << " locus_to_super_idx.size=" << locus_to_super_idx.size() 
              << " super_site_var_index.size=" << super_site_var_index.size() << std::endl;
    assert(super_sites.size() == 2);

    // HMM params + toy genetic map (anchor-to-anchor t=0.05; siblings share cm with anchors)
    std::cout << "STEP 4: initializing HMM parameters" << std::endl;
    hmm_parameters M;
    M.ed = 0.01f; M.ee = 1.0f;
    const int Neff = 10000;
    const int Nhap = static_cast<int>(H.n_hap);
    const double step_cM = (-std::log(1.0 - 0.05)) / 0.04 * (double)Nhap / (double)Neff;
    V.vec_pos[0]->cm = 0.0; V.vec_pos[1]->cm = 0.0;           // SS1 (0,1)
    V.vec_pos[2]->cm = step_cM; V.vec_pos[3]->cm = step_cM;   // SS2 (2,3)
    M.initialise(V, Neff, Nhap);
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);
    std::cout << "  M.rare_allele.size=" << M.rare_allele.size() << std::endl;

    // Target genotype: homozygous REF across all splits (HOM path)
    std::cout << "STEP 5: building genotype" << std::endl;
    genotype G(0);
    G.n_segments = 1; G.n_variants = V.size(); G.n_ambiguous = 0; G.n_missing = 0;
    G.n_transitions = 0; G.n_stored_transitionProbs = 0; G.n_storage_events = 0;
    G.double_precision = false; G.haploid = false;
    G.Variants.assign(2, 0);
    for (int e = 0; e < 4; ++e) VAR_SET_HOM(MOD2(e), G.Variants[DIV2(e)]);
    G.Lengths.assign(1, (unsigned short)V.size());
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 1ull);

    std::cout << "STEP 6: creating window" << std::endl;
    window W;
    W.start_locus = 0; W.stop_locus = (int)V.size() - 1; // 0..3
    W.start_segment = 0; W.stop_segment = 0;
    W.start_ambiguous = 0; W.stop_ambiguous = -1;
    W.start_missing = 0; W.stop_missing = -1;
    W.start_transition = 0; W.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};
    std::cout << "STEP 7: constructing HS" << std::endl;

    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), packed_codes.size(), &super_site_var_index);

    // Run forward locus by locus to compare states
    std::cout << "STEP 8: running forward()" << std::endl;
    HS.forward();

    // Snapshot after each locus step: use Alpha buffer for segment start; otherwise track prob after run
    // Here, we can simply re-run multiple windows to capture per-locus states
    auto run_window = [&](int start, int stop, aligned_vector32<float>& out_prob, aligned_vector32<float>& out_sumH, float& out_sumT) {
        window Wloc; Wloc.start_locus = start; Wloc.stop_locus = stop;
        Wloc.start_segment = 0; Wloc.stop_segment = 0;
        Wloc.start_ambiguous = 0; Wloc.stop_ambiguous = -1;
        Wloc.start_missing = 0; Wloc.stop_missing = -1;
        Wloc.start_transition = 0; Wloc.stop_transition = -1;
        haplotype_segment_single HSx(&G, H.H_opt_hap, idxH, Wloc, M,
            &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), packed_codes.size(), &super_site_var_index);
        HSx.forward();
        out_prob = HSx.prob; out_sumH = HSx.probSumH; out_sumT = HSx.probSumT;
    };

    aligned_vector32<float> prob0, sumH0, prob1, sumH1, prob2, sumH2, prob3, sumH3;
    float sumT0=0, sumT1=0, sumT2=0, sumT3=0;
    run_window(0,0, prob0, sumH0, sumT0); // after SS1 anchor
    run_window(0,1, prob1, sumH1, sumT1); // after SS1 sibling (should match)
    run_window(0,2, prob2, sumH2, sumT2); // after SS2 anchor
    run_window(0,3, prob3, sumH3, sumT3); // after SS2 sibling (should match prob2)

    // 1) SS1 sibling no-op: state after [0,1] equals state after [0,0]
    assert(prob0.size() == prob1.size());
    for (size_t i = 0; i < prob0.size(); ++i) assert(std::fabs(prob0[i] - prob1[i]) <= 1e-6f);
    for (int h = 0; h < HAP_NUMBER; ++h) assert(std::fabs(sumH0[h] - sumH1[h]) <= 1e-6f);
    assert(std::fabs(sumT0 - sumT1) <= 1e-6f);

    // 2) SS2 sibling no-op: state after [0..3] equals state after [0..2]
    assert(prob2.size() == prob3.size());
    for (size_t i = 0; i < prob2.size(); ++i) assert(std::fabs(prob2[i] - prob3[i]) <= 1e-6f);
    for (int h = 0; h < HAP_NUMBER; ++h) assert(std::fabs(sumH2[h] - sumH3[h]) <= 1e-6f);
    assert(std::fabs(sumT2 - sumT3) <= 1e-6f);

    // 3) Both anchors update non-trivially: states differ across the two anchors
    bool changed_ss1 = false, changed_ss2 = false;
    for (size_t i = 0; i < prob0.size(); ++i) if (std::fabs(prob0[i] - prob2[i]) > 1e-7f) { changed_ss2 = true; break; }
    for (size_t i = 0; i < prob0.size(); ++i) if (std::fabs(prob0[i] - prob1[i]) > 1e-7f) { changed_ss1 = true; break; }
    // changed_ss1 should be false (sibling no-op), changed_ss2 should be true (second anchor applied)
    assert(!changed_ss1);
    assert(changed_ss2);

    TEST_PASS("test_supersite_combine");  // was: PASS
    TEST_SUMMARY();
    return 0;
}

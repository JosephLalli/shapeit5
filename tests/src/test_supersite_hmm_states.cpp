#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <limits>

#include "../../common/src/utils/otools.h"

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

static void init_genotype(genotype& G,
                          unsigned int n_variants,
                          const std::vector<unsigned short>& segment_lengths,
                          unsigned int n_ambiguous,
                          unsigned int n_missing) {
    G.n_segments = static_cast<unsigned int>(segment_lengths.size());
    G.n_variants = n_variants;
    G.n_ambiguous = n_ambiguous;
    G.n_missing = n_missing;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;

    G.Variants.assign((n_variants + 1) / 2, 0);
    G.Ambiguous.assign(n_ambiguous, 0);
    G.Diplotypes.assign(G.n_segments, 0);
    G.Lengths.assign(G.n_segments, 0);
    for (size_t i = 0; i < segment_lengths.size(); ++i) {
        G.Lengths[i] = segment_lengths[i];
    }
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.assign(n_missing, 0.0f);
}

static void set_cond_alt(conditioning_set& H, int variant_idx, int hap_idx) {
    H.H_opt_var.set(static_cast<unsigned int>(variant_idx), static_cast<unsigned int>(hap_idx), 1);
    H.H_opt_hap.set(static_cast<unsigned int>(hap_idx), static_cast<unsigned int>(variant_idx), 1);
}

int main() {
    std::cout << "Testing supersite HMM state transitions (float vs double)..." << std::endl;

    // Build synthetic variant map: three supersites, each split into two biallelic records.
    variant_map V;
    V.push(make_var("1", 1000, "ssA_alt1", "A", "C", 0)); // index 0
    V.push(make_var("1", 1000, "ssA_alt2", "A", "G", 1)); // index 1
    V.push(make_var("1", 2000, "ssB_alt1", "A", "T", 2)); // index 2
    V.push(make_var("1", 2000, "ssB_alt2", "A", "G", 3)); // index 3
    V.push(make_var("1", 3000, "ssC_alt1", "A", "T", 4)); // index 4
    V.push(make_var("1", 3000, "ssC_alt2", "A", "G", 5)); // index 5

    const unsigned int n_variants = static_cast<unsigned int>(V.size());

    // Conditioning panel: one reference sample => two haplotypes.
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V.size());

    // Encode panel haplotypes so each supersite has distinct conditioning codes.
    set_cond_alt(H, 0, 0); // Supersite A, hap0 carries ALT1
    set_cond_alt(H, 2, 1); // Supersite B, hap1 carries ALT1
    set_cond_alt(H, 3, 0); // Supersite B, hap0 carries ALT2
    // Supersite C remains REF for both haps (all zeros).

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;

    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, locus_to_super_idx, super_site_var_index, unused_sample_codes);
    assert(super_sites.size() == 3);

    // Genotype: two segments (lengths 3 and 3) touching INIT/RUN/COLLAPSE with hom/amb/mis cases.
    genotype G(0);
    init_genotype(G, n_variants, {3, 3}, /*n_ambiguous*/2, /*n_missing*/2);

    // Supersite A (indices 0,1): homozygous ALT
    VAR_SET_HOM(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP1(MOD2(0), G.Variants[DIV2(0)]);

    VAR_SET_HOM(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HAP0(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HAP1(MOD2(1), G.Variants[DIV2(1)]);

    // Supersite B (indices 2,3): ambiguous - hap0 carries ALT1, hap1 carries ALT2
    VAR_SET_HET(MOD2(2), G.Variants[DIV2(2)]);
    VAR_SET_HAP0(MOD2(2), G.Variants[DIV2(2)]);
    VAR_CLR_HAP1(MOD2(2), G.Variants[DIV2(2)]);

    VAR_SET_HET(MOD2(3), G.Variants[DIV2(3)]);
    VAR_CLR_HAP0(MOD2(3), G.Variants[DIV2(3)]);
    VAR_SET_HAP1(MOD2(3), G.Variants[DIV2(3)]);

    // Ambiguous encodings (toy patterns sufficient for exercising AMB code paths)
    G.Ambiguous[0] = 0xAA; // alternate hap on odd lanes
    G.Ambiguous[1] = 0x55; // alternate hap on even lanes

    // Supersite C (indices 4,5): missing at both components
    VAR_SET_MIS(MOD2(4), G.Variants[DIV2(4)]);
    VAR_SET_MIS(MOD2(5), G.Variants[DIV2(5)]);

    // HMM parameters
    hmm_parameters M;
    M.ed = 0.1;
    M.ee = 1.0;
    M.t = std::vector<float>(n_variants ? n_variants - 1 : 0, 0.2f);
    M.nt = std::vector<float>(M.t.size(), 0.8f);
    M.rare_allele = std::vector<char>(n_variants, -1);

    window W;
    W.start_locus = 0;
    W.stop_locus = static_cast<int>(n_variants) - 1;
    W.start_segment = 0;
    W.stop_segment = 1;
    W.start_ambiguous = 0;
    W.stop_ambiguous = static_cast<int>(G.n_ambiguous) - 1;
    W.start_missing = 0;
    W.stop_missing = static_cast<int>(G.n_missing) - 1;
    W.start_transition = 0;
    W.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u};

    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), &super_site_var_index);
    haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), &super_site_var_index);

    HS.forward();
    HD.forward();

    auto assert_close = [](double lhs, double rhs, double tol, const char* label) {
        double diff = std::fabs(lhs - rhs);
        if (diff > tol) {
            std::cerr << label << " diff=" << diff << " lhs=" << lhs << " rhs=" << rhs << " tol=" << tol << std::endl;
            assert(false);
        }
    };

    const double tol = 1e-6;

    const std::array<double,16> expected_prob = {
        0.022937, 0.022937, 0.022937, 0.022937,
        0.151169, 0.151169, 0.151169, 0.151169,
        0.0392858, 0.0392858, 0.0392858, 0.0392858,
        0.0366085, 0.0366085, 0.0366085, 0.0366085
    };
    assert(HD.prob.size() == expected_prob.size());
    for (size_t i = 0; i < expected_prob.size(); ++i) {
        assert_close(HD.prob[i], expected_prob[i], tol, "HD.prob");
    }

    const std::array<double,8> expected_probSumH = {
        0.0622228, 0.0622228, 0.0622228, 0.0622228,
        0.187777, 0.187777, 0.187777, 0.187777
    };
    assert(HD.probSumH.size() == expected_probSumH.size());
    for (size_t i = 0; i < expected_probSumH.size(); ++i) {
        assert_close(HD.probSumH[i], expected_probSumH[i], tol, "HD.probSumH");
    }
    assert_close(HD.probSumT, 1.0, tol, "HD.probSumT");

    assert(HD.Alpha.size() == 2);
    const std::array<double,16> expected_alpha0 = {
        0.0110455, 0.0110455, 0.0110455, 0.0110455,
        0.110455, 0.110455, 0.110455, 0.110455,
        0.0145452, 0.0145452, 0.0145452, 0.0145452,
        0.00145452, 0.00145452, 0.00145452, 0.00145452
    };
    const std::array<double,16> expected_alpha1 = expected_prob;
    for (size_t i = 0; i < expected_alpha0.size(); ++i) {
        assert_close(HD.Alpha[0][i], expected_alpha0[i], tol, "HD.Alpha[0]");
        assert_close(HD.Alpha[1][i], expected_alpha1[i], tol, "HD.Alpha[1]");
    }

    const std::array<double,8> expected_alpha_sum0 = {
        0.0255907, 0.0255907, 0.0255907, 0.0255907,
        0.111909, 0.111909, 0.111909, 0.111909
    };
    const std::array<double,8> expected_alpha_sum1 = expected_probSumH;
    for (size_t i = 0; i < expected_alpha_sum0.size(); ++i) {
        assert_close(HD.AlphaSum[0][i], expected_alpha_sum0[i], tol, "HD.AlphaSum[0]");
        assert_close(HD.AlphaSum[1][i], expected_alpha_sum1[i], tol, "HD.AlphaSum[1]");
    }

    const std::array<double,2> expected_alpha_sum_sum = {0.55, 1.0};
    assert(HD.AlphaSumSum.size() == expected_alpha_sum_sum.size());
    for (size_t i = 0; i < expected_alpha_sum_sum.size(); ++i) {
        assert_close(HD.AlphaSumSum[i], expected_alpha_sum_sum[i], tol, "HD.AlphaSumSum");
    }

    assert(HD.AlphaMissing.size() == 2);
    const std::array<double,16> expected_alpha_missing0 = {
        0.0208934, 0.0208934, 0.0208934, 0.0208934,
        0.165489, 0.165489, 0.165489, 0.165489,
        0.0413294, 0.0413294, 0.0413294, 0.0413294,
        0.0222884, 0.0222884, 0.0222884, 0.0222884
    };
    const std::array<double,16> expected_alpha_missing1 = expected_prob;
    for (size_t i = 0; i < expected_alpha_missing0.size(); ++i) {
        assert_close(HD.AlphaMissing[0][i], expected_alpha_missing0[i], tol, "HD.AlphaMissing[0]");
        assert_close(HD.AlphaMissing[1][i], expected_alpha_missing1[i], tol, "HD.AlphaMissing[1]");
    }

    assert(HD.AlphaSumMissing.size() == 2);
    const std::array<double,8> expected_alpha_sum_missing0 = expected_probSumH;
    const std::array<double,8> expected_alpha_sum_missing1 = expected_probSumH;
    for (size_t i = 0; i < expected_alpha_sum_missing0.size(); ++i) {
        assert_close(HD.AlphaSumMissing[0][i], expected_alpha_sum_missing0[i], tol, "HD.AlphaSumMissing[0]");
        assert_close(HD.AlphaSumMissing[1][i], expected_alpha_sum_missing1[i], tol, "HD.AlphaSumMissing[1]");
    }

    for (float v : HS.prob) assert(std::isfinite(v));
    for (float v : HS.probSumH) assert(std::isfinite(v));
    assert(std::isfinite(HS.probSumT));

    std::vector<double> trans_single;
    std::vector<float> missing_single(G.n_missing, 0.0f);
    int recover_single = HS.backward(trans_single, missing_single);

    std::vector<double> trans_double;
    std::vector<float> missing_double(G.n_missing, 0.0f);
    int recover_double = HD.backward(trans_double, missing_double);

    assert(recover_single == recover_double);
    assert(trans_single.size() == trans_double.size());
    for (size_t i = 0; i < trans_single.size(); ++i) {
        assert_close(trans_single[i], trans_double[i], tol, "transition_prob");
    }
    assert(missing_single.size() == missing_double.size());
    for (size_t i = 0; i < missing_single.size(); ++i) {
        assert_close(static_cast<double>(missing_single[i]), static_cast<double>(missing_double[i]), tol, "missing_prob");
    }

    std::cout << "  OK" << std::endl;
    return 0;
}

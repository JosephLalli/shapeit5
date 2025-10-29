#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

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

static void setup_genotype(genotype& G, unsigned int n_variants, unsigned int segment_length) {
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
    G.Diplotypes.assign(1, 0);
    G.Lengths.assign(1, static_cast<unsigned short>(segment_length));
    G.ProbMask.clear();
    G.ProbStored.clear();
    G.ProbMissing.clear();
}

int main() {
    std::cout << "Testing supersite HMM harness..." << std::endl;

    // Synthetic variant map: two split records at the same position (supersite)
    variant_map V;
    V.push(make_var("1", 1000, "ss_1_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_1_A_G", "A", "G", 1));

    // Conditioning panel: 1 reference sample => 2 haplotypes
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V.size());

    // Hap0 carries ALT for the first split record, Hap1 stays REF
    H.H_opt_var.set(/*variant*/0, /*hap*/0, /*bit*/1);
    H.H_opt_hap.set(/*hap*/0, /*variant*/0, /*bit*/1);

    // Build supersites
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;

    std::cout << "  Building supersites" << std::endl;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, locus_to_super_idx, super_site_var_index, unused_sample_codes);

    assert(super_sites.size() == 1);
    assert(!packed_codes.empty());
    const SuperSite& ss = super_sites[0];
    assert(ss.var_count == 2);
    assert(locus_to_super_idx.size() == V.size());
    assert(locus_to_super_idx[0] == 0 && locus_to_super_idx[1] == 0);

    // Sample genotype: homozygous ALT at the first member variant, REF at the second
    genotype G(0);
    setup_genotype(G, /*n_variants*/V.size(), /*segment_length*/V.size());
    VAR_SET_HOM(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_HAP1(MOD2(0), G.Variants[DIV2(0)]);

    // HMM parameters with deterministic transition rate
    hmm_parameters M;
    M.ed = 0.1;
    M.ee = 1.0;
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.2f);
    M.nt = std::vector<float>(M.t.size(), 0.8f);
    M.rare_allele = std::vector<char>(V.size(), -1);

    // Single window covering both loci
    window W;
    W.start_locus = 0;
    W.stop_locus = static_cast<int>(V.size()) - 1;
    W.start_segment = 0;
    W.stop_segment = 0;
    W.start_ambiguous = 0;
    W.stop_ambiguous = -1;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = -1;

    std::vector<unsigned int> idxH = {0u, 1u};

    // Expected emissions after INIT (locus 0) and RUN (locus 1)
    const double mismatch = M.ed / M.ee; // 0.1
    const double sum_emit = 1.0 + mismatch;
    const double probSumT_init = HAP_NUMBER * sum_emit;
    const double yt = 0.2;
    const double nt = 1.0 - yt;
    const double tFreq = sum_emit * (yt / (static_cast<double>(idxH.size()) * probSumT_init));
    const double nt_norm = nt / probSumT_init;
    const double expect_cond0 = (1.0 * nt_norm + tFreq) * 1.0;
    const double expect_cond1 = (mismatch * nt_norm + tFreq) * mismatch;
    const double expect_sum = expect_cond0 + expect_cond1;
    const double expect_total = expect_sum * HAP_NUMBER;

    // Single-precision path
    std::cout << "  Running single-precision forward" << std::endl;
    haplotype_segment_single HS(
        &G,
        H.H_opt_hap,
        idxH,
        W,
        M,
        &super_sites,
        &is_super_site,
        &locus_to_super_idx,
        packed_codes.data(),
        &super_site_var_index
    );

    HS.forward();
    std::cout << "    single forward complete" << std::endl;

    assert(std::fabs(HS.prob[0] - static_cast<float>(expect_cond0)) < 1e-5f);
    assert(std::fabs(HS.prob[HAP_NUMBER] - static_cast<float>(expect_cond1)) < 1e-6f);
    assert(std::fabs(HS.probSumH[0] - static_cast<float>(expect_sum)) < 1e-6f);
    assert(std::fabs(HS.probSumT - static_cast<float>(expect_total)) < 1e-5f);

    // Double-precision path
    std::cout << "  Running double-precision forward" << std::endl;
    haplotype_segment_double HD(
        &G,
        H.H_opt_hap,
        idxH,
        W,
        M,
        &super_sites,
        &is_super_site,
        &locus_to_super_idx,
        packed_codes.data(),
        &super_site_var_index
    );

    HD.forward();
    std::cout << "    double forward complete" << std::endl;

    assert(std::fabs(HD.prob[0] - expect_cond0) < 1e-8);
    assert(std::fabs(HD.prob[HAP_NUMBER] - expect_cond1) < 1e-8);
    assert(std::fabs(HD.probSumH[0] - expect_sum) < 1e-8);
    assert(std::fabs(HD.probSumT - expect_total) < 1e-8);

    std::cout << "  OK" << std::endl;
    return 0;
}

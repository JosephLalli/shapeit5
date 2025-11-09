#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/containers/variant_map.h"

using std::cout;
using std::endl;

static variant *vmake(const std::string &chr, int bp, const std::string &id,
                      const std::string &ref, const std::string &alt, int idx) {
    std::string chr_copy(chr);
    std::string id_copy(id);
    std::string ref_copy(ref);
    std::string alt_copy(alt);
    return new variant(chr_copy, bp, id_copy, ref_copy, alt_copy, idx);
}

// Helper to seed a minimal panel with four haplotypes and two variants (anchor + sibling)
static void build_minimal_panel(conditioning_set &H, variant_map &V) {
    H.allocate(0, 2, V.size());
    // Simple pattern: toggle alleles across donors so row sums are non-uniform
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);
}

int main() {
    cout << "Testing supersite outer-product seeding..." << endl;

    variant_map V;
    V.push(vmake("1", 1000, "SS_A_C", "A", "C", 0)); // anchor (segment 0)
    V.push(vmake("1", 1000, "SS_A_G", "A", "G", 1)); // sibling (segment 1)

    conditioning_set H;
    build_minimal_panel(H, V);

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;

    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);
    assert(super_sites.size() == 1);

    hmm_parameters M;
    M.initialise(V, /*Neff=*/10000, /*Nhap=*/4);

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    genotype G(0);
    G.n_segments = 2;
    G.n_variants = V.size();
    G.Lengths.assign(2, 1);
    G.Variants.assign((G.n_variants + 1) / 2, 0u);
    G.Diplotypes.assign(2, 1ull);

    window W;
    W.start_locus = 0; W.stop_locus = 1;
    W.start_segment = 0; W.stop_segment = 1;
    W.start_ambiguous = 0; W.stop_ambiguous = -1;
    W.start_missing = 0; W.stop_missing = -1;
    W.start_transition = 0; W.stop_transition = 0;

    haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx,
        packed_codes.data(), packed_codes.size(), &super_site_var_index);

    haplotype_segment_double HD(&G, H.H_opt_hap, idxH, W, M,
        &super_sites, &is_super_site, &locus_to_super_idx,
        packed_codes.data(), packed_codes.size(), &super_site_var_index);

    const int prev_seg = 0;
    const float prev_total_f = 2.0f;
    const double prev_total_d = 2.0;
    const float nt_f = 0.8f;
    const float yt_f = 0.2f;
    const double nt_d = 0.8;
    const double yt_d = 0.2;

    // Seed lane marginals (column sums) and totals for previous segment
    std::fill(HS.AlphaLaneSum[prev_seg].begin(), HS.AlphaLaneSum[prev_seg].end(), 0.0f);
    std::fill(HS.AlphaSum[prev_seg].begin(), HS.AlphaSum[prev_seg].end(), 0.0f);
    HS.AlphaLaneSum[prev_seg][0] = 1.5f;
    HS.AlphaLaneSum[prev_seg][1] = 0.5f;
    HS.AlphaSum[prev_seg][0] = 1.5f;
    HS.AlphaSum[prev_seg][1] = 0.5f;
    HS.AlphaSumSum[prev_seg] = prev_total_f;

    std::fill(HD.AlphaLaneSum[prev_seg].begin(), HD.AlphaLaneSum[prev_seg].end(), 0.0);
    std::fill(HD.AlphaSum[prev_seg].begin(), HD.AlphaSum[prev_seg].end(), 0.0);
    HD.AlphaLaneSum[prev_seg][0] = 1.5;
    HD.AlphaLaneSum[prev_seg][1] = 0.5;
    HD.AlphaSum[prev_seg][0] = 1.5;
    HD.AlphaSum[prev_seg][1] = 0.5;
    HD.AlphaSumSum[prev_seg] = prev_total_d;

    // Seed donor marginals (row sums)
    const float row_sums_f[4] = {1.2f, 0.4f, 0.2f, 0.2f};
    const double row_sums_d[4] = {1.2, 0.4, 0.2, 0.2};
    for (int k = 0; k < 4; ++k) {
        HS.probSumK[k] = row_sums_f[k];
        HD.probSumK[k] = row_sums_d[k];
    }

    HS.nt = nt_f; HS.yt = yt_f; HS.probSumT = 1.0f;
    HD.nt = nt_d; HD.yt = yt_d; HD.probSumT = 1.0;
    HS.curr_segment_index = HS.segment_first + 1;
    HD.curr_segment_index = HD.segment_first + 1;

    // --- Check column mixing ---
    __m256 col_mix_f;
    float row_stay_f = 0.0f, row_switch_f = 0.0f;
    bool use_outer_f = HS.prepare_outer_product_mix(prev_seg, col_mix_f, row_stay_f, row_switch_f);
    assert(use_outer_f);

    alignas(32) float col_vals[HAP_NUMBER];
    _mm256_store_ps(col_vals, col_mix_f);

    const float expected_row_stay = nt_f / prev_total_f;
    const float expected_row_switch = yt_f / static_cast<float>(HS.n_cond_haps);
    assert(std::fabs(row_stay_f - expected_row_stay) < 1e-6f);
    assert(std::fabs(row_switch_f - expected_row_switch) < 1e-6f);

    const float inv_lanes = yt_f / static_cast<float>(HAP_NUMBER);
    const float expected_col0 = expected_row_stay * 1.5f + inv_lanes;
    const float expected_col1 = expected_row_stay * 0.5f + inv_lanes;
    const float expected_col_other = inv_lanes;
    assert(std::fabs(col_vals[0] - expected_col0) < 1e-6f);
    assert(std::fabs(col_vals[1] - expected_col1) < 1e-6f);
    for (int h = 2; h < HAP_NUMBER; ++h) {
        assert(std::fabs(col_vals[h] - expected_col_other) < 1e-6f);
    }

    // --- Run COLLAPSE_MIS in single precision and validate full outer product ---
    HS.COLLAPSE_MIS();
    const float row_mix_expected_f[4] = {
        expected_row_stay * row_sums_f[0] + expected_row_switch,
        expected_row_stay * row_sums_f[1] + expected_row_switch,
        expected_row_stay * row_sums_f[2] + expected_row_switch,
        expected_row_stay * row_sums_f[3] + expected_row_switch
    };
    for (int k = 0; k < 4; ++k) {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const float expected = row_mix_expected_f[k] * col_vals[h];
            const float observed = HS.prob[k * HAP_NUMBER + h];
            assert(std::fabs(observed - expected) < 1e-6f);
        }
    }

    // --- Double precision path should mirror single precision ---
    HD.COLLAPSE_MIS();
    const double expected_row_stay_d = nt_d / prev_total_d;
    const double expected_row_switch_d = yt_d / static_cast<double>(HD.n_cond_haps);
    const double col_expected_d[HAP_NUMBER] = {
        expected_row_stay_d * 1.5 + yt_d / static_cast<double>(HAP_NUMBER),
        expected_row_stay_d * 0.5 + yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER),
        yt_d / static_cast<double>(HAP_NUMBER)
    };
    for (int k = 0; k < 4; ++k) {
        const double row_mix = expected_row_stay_d * row_sums_d[k] + expected_row_switch_d;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const double expected = row_mix * col_expected_d[h];
            const double observed = HD.prob[k * HAP_NUMBER + h];
            assert(std::fabs(observed - expected) < 1e-9);
        }
    }

    cout << "PASS: outer-product seeding preserved in single and double precision" << endl;
    return 0;
}

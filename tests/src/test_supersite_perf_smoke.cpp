 /*******************************************************************************
  * Supersite performance smoke test (non-failing)
  *
  * Runs a number of forward() passes with and without supersites and prints
  * elapsed times. No assertions, diagnostic only.
  ******************************************************************************/

#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

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
    TEST_INIT("test_supersite_perf_smoke");
    using clock = std::chrono::high_resolution_clock;
    std::cout << "Supersite performance smoke test..." << std::endl;

    // Build a window with 100 loci; every 10th locus is a 2-split supersite
    variant_map V;
    int idx = 0;
    for (int pos = 1000; pos < 2000; ++pos) {
        if ((pos % 10) == 0) {
            V.push(make_var("1", pos, "ss_A_C", "A", "C", idx++));
            V.push(make_var("1", pos, "ss_A_G", "A", "G", idx++));
        } else {
            V.push(make_var("1", pos, "b", "A", "T", idx++));
        }
    }

    conditioning_set H;
    H.allocate(0, 2, V.size());
    // Seed some panel alts
    for (size_t v = 0; v < V.size(); ++v) if ((v % 7) == 0) { H.H_opt_var.set(v, 1, 1); H.H_opt_hap.set(1, v, 1); }

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index);

    hmm_parameters M;
    M.ed = 0.01f; M.ee = 1.0f;
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
    M.nt = std::vector<float>(M.t.size(), 0.95f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

    genotype G(0);
    G.n_variants = static_cast<unsigned int>(V.size());
    G.Variants.assign((V.size() + 1) / 2, 0);
    for (size_t v = 0; v < V.size(); ++v) {
        if ((v % 5) == 0) {
            VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
            if (((v / 5) % 2) == 0) {
                VAR_CLR_HAP0(MOD2(v), G.Variants[DIV2(v)]);
                VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);
            } else {
                VAR_SET_HAP0(MOD2(v), G.Variants[DIV2(v)]);
                VAR_CLR_HAP1(MOD2(v), G.Variants[DIV2(v)]);
            }
        } else {
            VAR_SET_HOM(MOD2(v), G.Variants[DIV2(v)]);
        }
    }
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());
    G.build();

    window W;
    W.start_locus = 0; W.stop_locus = (int)V.size() - 1;
    W.start_segment = 0; W.stop_segment = static_cast<int>(G.n_segments) - 1;
    W.start_ambiguous = 0; W.stop_ambiguous = G.n_ambiguous ? static_cast<int>(G.n_ambiguous) - 1 : -1;
    W.start_missing = 0; W.stop_missing = G.n_missing ? static_cast<int>(G.n_missing) - 1 : -1;
    W.start_transition = 0; W.stop_transition = G.n_transitions ? static_cast<int>(G.n_transitions) - 1 : -1;
    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    // Measure supersite-enabled
    auto t1 = clock::now();
    for (int r = 0; r < 10; ++r) {
        haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
        HS.forward();
    }
    auto t2 = clock::now();
    auto ms_ss = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

    // Measure biallelic-only (no supersite pointers)
    G.setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(nullptr, 0);
    t1 = clock::now();
    for (int r = 0; r < 10; ++r) {
        haplotype_segment_single HS(&G, H.H_opt_hap, idxH, W, M);
        HS.forward();
    }
    t2 = clock::now();
    auto ms_bi = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

    std::cout << "Elapsed (supersites on):  " << ms_ss << " ms\n";
    std::cout << "Elapsed (biallelic only): " << ms_bi << " ms\n";
    std::cout << "(Perf smoke; no assertions)\n";
    TEST_SUMMARY();
    return 0;
}

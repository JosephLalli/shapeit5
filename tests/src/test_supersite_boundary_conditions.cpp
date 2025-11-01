 /*******************************************************************************
  * Supersite boundary condition tests
  *
  * - INIT boundary: window starts on a sibling → neutral (uniform) initialization
  * - COLLAPSE boundary: window ends on a sibling → sibling is no-op; same state
  ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"

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

static void setup_genotype(genotype& G, unsigned int n_variants) {
    G.n_segments = 1;
    G.n_variants = n_variants;
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    G.Variants.assign((n_variants + 1)/2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull);
    G.Lengths.assign(1, (unsigned short)n_variants);
    G.ProbMask.clear();
    G.ProbStored.clear();
}

int main() {
    std::cout << "Testing supersite boundary conditions (INIT/COLLAPSE)..." << std::endl;

    // Build 2-split supersite
    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0)); // anchor
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1)); // sibling

    conditioning_set H;
    H.allocate(0, 2, V.size()); // 4 haplotypes
    // Diverse panel to avoid degeneracy
    H.H_opt_var.set(0, 1, 1); H.H_opt_hap.set(1, 0, 1);
    H.H_opt_var.set(1, 2, 1); H.H_opt_hap.set(2, 1, 1);

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> unused_sample_codes;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, unused_sample_codes);
    assert(super_sites.size() == 1);

    // HMM params and sibling marking
    hmm_parameters M;
    M.ed = 0.01f; M.ee = 1.0f;
    M.t = std::vector<float>(V.size() ? V.size() - 1 : 0, 0.05f);
    M.nt = std::vector<float>(M.t.size(), 0.95f);
    M.rare_allele = std::vector<char>(V.size(), -1);
    M.markSuperSiteSiblings(super_sites, locus_to_super_idx);

    std::vector<unsigned int> idxH = {0u, 1u, 2u, 3u};

    // Target genotype: homozygous REF
    genotype G(0);
    setup_genotype(G, V.size());
    VAR_SET_HOM(0, G.Variants[0]);
    VAR_SET_HOM(1, G.Variants[0]);

    // 1) INIT boundary: window starts on sibling only
    window W_init;
    W_init.start_locus = 1; W_init.stop_locus = 1; // sibling only
    W_init.start_segment = 0; W_init.stop_segment = 0;
    W_init.start_ambiguous = 0; W_init.stop_ambiguous = -1;
    W_init.start_missing = 0; W_init.stop_missing = -1;
    W_init.start_transition = 0; W_init.stop_transition = -1;

    haplotype_segment_single HS_init(&G, H.H_opt_hap, idxH, W_init, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), &super_site_var_index);
    HS_init.forward();
    // Expect neutral (uniform) init to prevent underflow
    assert(std::fabs(HS_init.probSumT - 1.0f) < 1e-6f);

    // 2) COLLAPSE boundary: compare anchor-only vs anchor+sibling windows
    window W_anchor;
    W_anchor.start_locus = 0; W_anchor.stop_locus = 0; // anchor only
    W_anchor.start_segment = 0; W_anchor.stop_segment = 0;
    W_anchor.start_ambiguous = 0; W_anchor.stop_ambiguous = -1;
    W_anchor.start_missing = 0; W_anchor.stop_missing = -1;
    W_anchor.start_transition = 0; W_anchor.stop_transition = -1;

    haplotype_segment_single HS_anchor(&G, H.H_opt_hap, idxH, W_anchor, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), &super_site_var_index);
    HS_anchor.forward();

    window W_both;
    W_both.start_locus = 0; W_both.stop_locus = 1; // anchor + sibling
    W_both.start_segment = 0; W_both.stop_segment = 0;
    W_both.start_ambiguous = 0; W_both.stop_ambiguous = -1;
    W_both.start_missing = 0; W_both.stop_missing = -1;
    W_both.start_transition = 0; W_both.stop_transition = -1;

    haplotype_segment_single HS_both(&G, H.H_opt_hap, idxH, W_both, M,
        &super_sites, &is_super_site, &locus_to_super_idx, packed_codes.data(), &super_site_var_index);
    HS_both.forward();

    // Sibling at end is a no-op → identical state
    assert(HS_both.prob.size() == HS_anchor.prob.size());
    for (size_t i = 0; i < HS_anchor.prob.size(); ++i) {
        assert(std::fabs(HS_both.prob[i] - HS_anchor.prob[i]) <= 1e-6f);
    }
    for (int h = 0; h < HAP_NUMBER; ++h) {
        assert(std::fabs(HS_both.probSumH[h] - HS_anchor.probSumH[h]) <= 1e-6f);
    }
    assert(std::fabs(HS_both.probSumT - HS_anchor.probSumT) <= 1e-6f);

    std::cout << "✓ SUCCESS: Boundary conditions validated (INIT neutral, COLLAPSE no-op)" << std::endl;
    return 0;
}


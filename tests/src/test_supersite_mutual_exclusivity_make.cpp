 /*******************************************************************************
  * Supersite mutual exclusivity projection test
  *
  * Verifies that genotype::make() projects multivariant classes to split records
  * such that exactly one ALT is set per haplotype across all member variants.
  ******************************************************************************/

#include <cassert>
#include <cmath>
#include <limits>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
    std::cout << "Testing supersite mutual exclusivity in genotype::make()..." << std::endl;

    // Create a 2-split supersite
    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V.size());

    // Build supersite metadata
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> sample_codes_unused;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, sample_codes_unused);
    assert(super_sites.size() == 1);
    SuperSite& ss = super_sites[0];

    // Genotype with both splits missing
    genotype G(0);
    G.n_segments = 1;
    G.n_variants = V.size();
    G.n_ambiguous = 0;
    G.n_missing = V.size();
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    G.Variants.assign((V.size() + 1)/2, 0);
    G.Lengths.assign(1, (unsigned short)V.size());
    G.Diplotypes.assign(1, 0);
    VAR_SET_MIS(0, G.Variants[0]);
    VAR_SET_MIS(1, G.Variants[0]);

    // Prepare SC: C = 1 + n_alts = 3; offset=0; prefer ALT1 for hap0, ALT2 for hap1
    ss.n_classes = 1 + ss.n_alts; // 3
    uint32_t offset = 1;  // Local offset for test (+1 guard)
    std::vector<float> SC(HAP_NUMBER * ss.n_classes + 2, 0.0f);
    SC.front() = std::numeric_limits<float>::quiet_NaN();
    SC.back() = std::numeric_limits<float>::quiet_NaN();
    std::vector<uint32_t> supersite_sc_offset(1, offset);  // Thread-local offset vector for test
    auto set_row = [&](int hap, float p_ref, float p_alt1, float p_alt2) {
        SC[offset + hap*ss.n_classes + 0] = p_ref;
        SC[offset + hap*ss.n_classes + 1] = p_alt1;
        SC[offset + hap*ss.n_classes + 2] = p_alt2;
    };
    // hap0 → ALT1, hap1 → ALT2, others uniform
    set_row(0, 0.0f, 1.0f, 0.0f);
    set_row(1, 0.0f, 0.0f, 1.0f);
    for (int h = 2; h < HAP_NUMBER; ++h) set_row(h, 1.0f, 0.0f, 0.0f);

    // Mark supersite as missing at anchor
    std::vector<bool> anchor_has_missing(1, true);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, &SC, &anchor_has_missing, &supersite_sc_offset);

    // Choose hap indices for this segment: hap0=0, hap1=1 → diplotype index d = 0*8 + 1 = 1
    std::vector<unsigned char> DipSampled(1, 1);
    std::vector<float> Mprob; // unused for supersite path

    // Project multivariant sample to splits
    G.make(DipSampled, Mprob);

    // Validate mutual exclusivity: exactly one ALT per hap across the two splits
    unsigned char vbyte = G.Variants[0];
    bool h0_v0 = VAR_GET_HAP0(0, vbyte);
    bool h0_v1 = VAR_GET_HAP0(1, vbyte);
    bool h1_v0 = VAR_GET_HAP1(0, vbyte);
    bool h1_v1 = VAR_GET_HAP1(1, vbyte);

    // hap0 should carry ALT at split0 only (ALT1)
    assert(h0_v0 == true);
    assert(h0_v1 == false);
    // hap1 should carry ALT at split1 only (ALT2)
    assert(h1_v0 == false);
    assert(h1_v1 == true);

    std::cout << "✓ SUCCESS: Mutual exclusivity enforced across split records" << std::endl;
    return 0;
}


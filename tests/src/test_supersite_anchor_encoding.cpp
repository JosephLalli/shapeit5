/*******************************************************************************
 * Test supersite anchor encoding update
 * Verifies that updateSuperSiteAnchorEncoding() marks the supersite anchor as
 * HET when the multiallelic site is heterozygous across splits (e.g., ALT1|ALT2),
 * while leaving HOM (or non-heterozygous) anchors unchanged.
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/containers/genotype_set.h"



#include "test_common.h"
using namespace std;

// Minimal helper to build a SuperSite spanning two split variants at positions 0 and 1
static SuperSite make_test_supersite(uint32_t anchor_id, uint32_t var_start, uint16_t var_count, uint8_t n_alts) {
    SuperSite ss{};
    ss.global_site_id = anchor_id;
    ss.chr = 0;
    ss.bp = 0;
    ss.n_alts = n_alts;
    ss.panel_offset = 0; // not used by updateSuperSiteAnchorEncoding
    ss.var_start = var_start;
    ss.var_count = var_count;
    ss.n_classes = static_cast<uint8_t>(1 + n_alts);
    return ss;
}

int main() {
    TEST_INIT("test_supersite_anchor_encoding");
    // One sample, two split variants (indices 0 and 1) forming one supersite
    genotype_set G;
    G.allocate(/*n_main_samples*/1, /*n_variants*/2);
    genotype* g = G.vecG[0];

    // Sanity: anchor (variant 0) starts as HOM 0/0
    unsigned char v0 = g->Variants[DIV2(0)];
    assert(VAR_GET_HOM(MOD2(0), v0));
    assert(!VAR_GET_HET(MOD2(0), v0));

    // Build supersite metadata: one supersite covering variants {0,1}, anchor at 0
    vector<SuperSite> super_sites(1);
    super_sites[0] = make_test_supersite(/*anchor_id*/0, /*var_start*/0, /*var_count*/2, /*n_alts*/2);
    vector<int> super_site_var_index = {0, 1};
    vector<int> locus_to_super_idx = {0, 0}; // both 0 and 1 map to supersite 0

    // Set supersite context on genotype (required by updateSuperSiteAnchorEncoding assertions)
    g->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index,
                           /*SC*/nullptr, /*anchor_has_missing*/nullptr, /*offsets*/nullptr);

    // Case 1: Heterozygous across splits (ALT1 on hap0 at split 0, ALT2 on hap1 at split 1)
    // Encode: variant 0 → hap0=1, variant 1 → hap1=1; anchor remains 0/0 pre-update
    VAR_SET_HAP0(MOD2(0), g->Variants[DIV2(0)]); // split0: hap0 carries ALT
    VAR_SET_HAP1(MOD2(1), g->Variants[DIV2(1)]); // split1: hap1 carries ALT

    // Verify classification helper detects heterozygous supersite
    assert(isSuperSiteHeterozygous(g, super_sites[0], super_site_var_index));

    // Run update – should mark anchor (variant 0) as HET
    updateSuperSiteAnchorEncoding(G, super_sites, super_site_var_index);
    v0 = g->Variants[DIV2(0)];
    assert(VAR_GET_HET(MOD2(0), v0));

    // Case 2: Non-heterozygous (HOM ALT1|ALT1) – anchor should not be forced to HET
    // Reset variants
    G.allocate(/*n_main_samples*/1, /*n_variants*/2);
    g = G.vecG[0];
    g->setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index,
                           /*SC*/nullptr, /*anchor_has_missing*/nullptr, /*offsets*/nullptr);
    // Both haplotypes carry ALT at split 0; no ALT at split 1
    VAR_SET_HAP0(MOD2(0), g->Variants[DIV2(0)]);
    VAR_SET_HAP1(MOD2(0), g->Variants[DIV2(0)]);
    // Pre-check: supersite is not heterozygous (c0==c1)
    assert(!isSuperSiteHeterozygous(g, super_sites[0], super_site_var_index));
    // Update should not set HET at anchor
    updateSuperSiteAnchorEncoding(G, super_sites, super_site_var_index);
    v0 = g->Variants[DIV2(0)];
    assert(!VAR_GET_HET(MOD2(0), v0));

    cout << "test_supersite_anchor_encoding: PASS" << endl;
    TEST_SUMMARY();
    return 0;
}

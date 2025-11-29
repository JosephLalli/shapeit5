 /*******************************************************************************
  * Supersite class microcases (HOM / AMB / MIS) classification tests
  ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"


#include "test_reporting.h"
static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
    TEST_INIT("test_supersite_class_microcases");
    std::cout << "Testing supersite class microcases (HOM/AMB/MIS)..." << std::endl;

    variant_map V;
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

    conditioning_set H;
    H.allocate(0, 1, V.size());

    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> sample_codes_unused;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
                    locus_to_super_idx, super_site_var_index, sample_codes_unused);
    assert(super_sites.size() == 1);
    const SuperSite& ss = super_sites[0];

    // Genotype setup
    genotype G(0);
    G.n_segments = 1; G.n_variants = V.size(); G.Variants.assign(1, 0);

    // Case 1: HOM (both haps same class) → ALT1 on both haps at first split
    G.Variants[0] = 0;
    VAR_SET_HAP0(0, G.Variants[0]); VAR_SET_HAP1(0, G.Variants[0]);
    VAR_CLR_HAP0(1, G.Variants[0]); VAR_CLR_HAP1(1, G.Variants[0]);
    {
        uint8_t c0, c1; SSClass cls = classify_supersite(&G, ss, super_site_var_index, c0, c1);
        assert(cls == SSClass::HOM);
    }

    // Case 2: AMB (c0 != c1): hap0 ALT1, hap1 ALT2
    G.Variants[0] = 0;
    VAR_SET_HAP0(0, G.Variants[0]); VAR_CLR_HAP1(0, G.Variants[0]);
    VAR_CLR_HAP0(1, G.Variants[0]); VAR_SET_HAP1(1, G.Variants[0]);
    {
        uint8_t c0, c1; SSClass cls = classify_supersite(&G, ss, super_site_var_index, c0, c1);
        assert(cls == SSClass::AMB);
    }

    // Case 3: MIS (both missing across splits)
    G.Variants[0] = 0;
    VAR_SET_MIS(0, G.Variants[0]); VAR_SET_MIS(1, G.Variants[0]);
    {
        uint8_t c0, c1; SSClass cls = classify_supersite(&G, ss, super_site_var_index, c0, c1);
        assert(cls == SSClass::MIS);
    }

    std::cout << "✓ SUCCESS: Supersite classification microcases OK" << std::endl;
    TEST_SUMMARY();
    return 0;
}


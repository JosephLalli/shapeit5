#include <cassert>
#include <iostream>
#include <vector>

#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"

int main() {
    std::cout << "Testing getSampleSuperSiteAlleleCode (per-hap)...\n" << std::flush;

    // Supersite of 3 split biallelic records at the same position
    SuperSite ss;
    ss.global_site_id = 0;
    ss.chr = 1;
    ss.bp = 1000;
    ss.n_alts = 3;
    ss.panel_offset = 0;
    ss.var_start = 0;
    ss.var_count = 3;

    std::vector<int> super_site_var_index = {0,1,2};

    // Hardcode a genotype with 3 variants
    genotype G(0);
    G.n_variants = 3;
    G.Variants.assign((G.n_variants + 1)/2, 0);

    // Hap0 carries ALT at variant 1, hap1 is REF everywhere
    VAR_SET_HET(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_HAP0(MOD2(1), G.Variants[DIV2(1)]);
    VAR_CLR_HAP1(MOD2(1), G.Variants[DIV2(1)]);

    uint8_t c0 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 0);
    uint8_t c1 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 1);
    // ALT index 1 -> class code 2
    assert(c0 == 2);
    assert(c1 == 0);

    // All missing across all three -> missing code
    G.Variants.assign((G.n_variants + 1)/2, 0);
    VAR_SET_MIS(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_MIS(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_MIS(MOD2(2), G.Variants[DIV2(2)]);
    c0 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 0);
    c1 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 1);
    assert(c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING);

    std::cout << "  OK\n";
    return 0;
}

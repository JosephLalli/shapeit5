#ifndef SHAPEIT5_TEST_FIXTURES_H
#define SHAPEIT5_TEST_FIXTURES_H

#include <vector>

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"

namespace TestFixtures {

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

inline SuperSiteContext build_supersites(variant_map& V, conditioning_set& H, bool clear_rare_mask = false) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    if (clear_rare_mask) {
        for (auto& ss : ctx.super_sites) {
            ss.rare_code_mask = {0u, 0u, 0u, 0u};
        }
    }
    return ctx;
}

inline void attach_supersite_context(genotype& G, const SuperSiteContext& ctx) {
    G.setSuperSiteContext(&ctx.super_sites, &ctx.locus_to_super_idx, &ctx.super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(ctx.packed_codes.data(), ctx.packed_codes.size());
}

inline void apply_supersite_pbwt_context(conditioning_set& H,
                                        const SuperSiteContext& ctx,
                                        size_t n_loci) {
    if (ctx.super_sites.empty()) {
        H.setSupersiteAnchorRedirect({});
        H.clearSupersitePBWTContext();
        return;
    }
    H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
    std::vector<int> anchor_map = buildSupersiteAnchorMap(ctx.super_sites, ctx.super_site_var_index, n_loci);
    H.setSupersiteAnchorRedirect(anchor_map);
    H.setSupersitePBWTContext(&ctx.super_sites, &ctx.locus_to_super_idx, &ctx.packed_codes);
}

} // namespace TestFixtures

#endif // SHAPEIT5_TEST_FIXTURES_H

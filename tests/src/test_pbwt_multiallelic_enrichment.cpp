/*******************************************************************************
 * PBWT multiallelic enrichment test
 *
 * Builds a small synthetic panel with one supersite (3 ALT classes) and tests
 * whether PBWT neighbors for a target hap are enriched for the same allele
 * class. Donor classes are decoded from packed_allele_codes (exercise packing).
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>
#include <set>

#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"

#include "test_reporting.h"

namespace {

enum AlleleClass : uint8_t { REF = 0, ALT1 = 1, ALT2 = 2, ALT3 = 3 };

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

struct PanelLayout {
    std::vector<int> alt1_markers;
    std::vector<int> alt2_markers;
    int anchor_alt1 = -1;
    int anchor_alt2 = -1;
    int anchor_alt3 = -1;
};

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx, double cm) {
    variant* v = new variant(chr, bp, id, ref, alt, idx);
    v->cm = cm;
    return v;
}

// Build a longer map with two class-specific marker blocks followed by a supersite anchor
// (ALT1/ALT2/ALT3 splits). cM increases linearly to keep PBWT grouping sensible.
static PanelLayout build_variant_map(variant_map& V) {
    PanelLayout layout;
    V.vec_pos.clear();

    const int n_alt1_markers = 6;
    const int n_alt2_markers = 20;
    const double cm_step = 0.001; // ~1 cM/Mb if bp step is 1kb
    int idx = 0;
    int bp = 1000000;

    auto add_block = [&](int count, const std::string& prefix, std::vector<int>& out_idx) {
        for (int i = 0; i < count; ++i) {
            double cm = (idx + i) * cm_step;
            V.push(make_var("1", bp + (idx + i) * 1000, prefix + std::to_string(i), "A", "C", idx + i, cm));
            out_idx.push_back(idx + i);
        }
        idx += count;
    };

    add_block(n_alt1_markers, "alt1_block_", layout.alt1_markers);
    add_block(n_alt2_markers, "alt2_block_", layout.alt2_markers);

    const int anchor_bp = bp + (idx * 1000);
    const double anchor_cm = idx * cm_step;
    layout.anchor_alt1 = idx++;
    V.push(make_var("1", anchor_bp, "ss_alt1", "A", "G", layout.anchor_alt1, anchor_cm));
    layout.anchor_alt2 = idx++;
    V.push(make_var("1", anchor_bp, "ss_alt2", "A", "T", layout.anchor_alt2, anchor_cm));
    layout.anchor_alt3 = idx++;
    V.push(make_var("1", anchor_bp, "ss_alt3", "A", "C", layout.anchor_alt3, anchor_cm));

    return layout;
}

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    return ctx;
}

// Decode donor class from packed codes at this supersite
static AlleleClass donor_class_from_packed(const std::vector<uint8_t>& packed_codes, const SuperSite& ss, unsigned hap) {
    uint32_t byte_idx = ss.panel_offset + hap / 2;
    uint8_t byte = packed_codes[byte_idx];
    uint8_t shift = (hap % 2) ? 4 : 0;
    return static_cast<AlleleClass>((byte >> shift) & 0x0F);
}

// Extract neighbors for a specific locus/target hap using the PBWT layout.
static std::vector<unsigned int> extract_neighbors(const conditioning_set& H, int locus, int target_hap) {
    std::vector<unsigned int> donors;
    if (locus >= (int)H.sites_pbwt_selection.size() || !H.sites_pbwt_selection[locus]) return donors;
    unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
    unsigned long base_idx = static_cast<unsigned long>(target_hap) * H.sites_pbwt_ngroups + H.sites_pbwt_grouping[locus];
    for (int d = 0; d < H.depth && base_idx < H.indexes_pbwt_neighbour.size(); ++d) {
        unsigned long neighbor_idx = static_cast<unsigned long>(d) * addr_offset + base_idx;
        if (neighbor_idx < H.indexes_pbwt_neighbour.size()) {
            int donor_hap = H.indexes_pbwt_neighbour[neighbor_idx];
            if (donor_hap >= 0) donors.push_back((unsigned int)donor_hap);
        }
    }
    return donors;
}

struct TargetResult {
    AlleleClass target;
    double alt1_frac;
    double alt2_frac;
    double alt3_frac;
    double ref_frac;
};

static TargetResult run_case(AlleleClass target_class) {
    rng.setSeed(12345); // deterministic selection when randomness is used
    const bool trace = std::getenv("SHAPEIT5_TEST_TRACE");

    variant_map V;
    PanelLayout layout = build_variant_map(V);

    const int n_individuals = 10; // 1 target + 9 donors
    conditioning_set H;
    H.allocate(/*n_main_samples=*/n_individuals, /*n_ref_samples=*/0, /*n_variants=*/V.size());

    // Assign classes per hap: indices 0/1 belong to target; then donors in blocks.
    std::vector<AlleleClass> hap_class(H.n_hap, REF);
    hap_class[0] = target_class; // target hap0
    hap_class[1] = REF;          // target hap1
    int h = 2;
    for (int i = 0; i < 3; ++i, h += 2) hap_class[h] = hap_class[h+1] = ALT1; // 3 ALT1 donors
    for (int i = 0; i < 3; ++i, h += 2) hap_class[h] = hap_class[h+1] = ALT2; // 3 ALT2 donors
    for (int i = 0; i < 3; ++i, h += 2) hap_class[h] = hap_class[h+1] = REF;  // 3 REF donors
    assert(h == (int)H.n_hap);

    // Populate hap-major matrix with class-specific marker blocks and anchor splits
    for (int hap = 0; hap < (int)H.n_hap; ++hap) {
        AlleleClass cls = hap_class[hap];
        bool is_alt1 = (cls == ALT1);
        bool is_alt2 = (cls == ALT2);
        bool is_alt3 = (cls == ALT3);
        for (int idx : layout.alt1_markers) H.H_opt_hap.set(hap, idx, is_alt1);
        for (int idx : layout.alt2_markers) H.H_opt_hap.set(hap, idx, is_alt2);
        H.H_opt_hap.set(hap, layout.anchor_alt1, is_alt1);
        H.H_opt_hap.set(hap, layout.anchor_alt2, is_alt2);
        H.H_opt_hap.set(hap, layout.anchor_alt3, is_alt3);
    }

    // Variant-major view for supersite packing/PBWT
    H.transposeHaplotypes_H2V(true, false);

    SuperSiteContext ctx = build_supersites(V, H);
    assert(ctx.super_sites.size() == 1);
    const SuperSite& ss = ctx.super_sites[0];
    const int anchor = ss.global_site_id;
    assert(ss.var_count == 3);

    // Initialize PBWT (all sites selected) and mask supersite siblings
    H.initialize_for_test(V, /*depth=*/8, /*nthread=*/1);
    H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
    std::vector<int> anchor_map = buildSupersiteAnchorMap(ctx.super_sites, ctx.super_site_var_index, V.size());
    H.setSupersiteAnchorRedirect(anchor_map);

    H.select();

    std::vector<unsigned int> neigh = extract_neighbors(H, anchor, /*target_hap=*/0);
    if (trace) {
        std::cout << "Anchor " << anchor << " selected=" << H.sites_pbwt_selection[anchor]
                  << " neighbors:";
        for (auto hidx : neigh) {
            std::cout << " " << hidx;
        }
        std::cout << " | raw:";
        unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
        unsigned long base_idx = static_cast<unsigned long>(0) * H.sites_pbwt_ngroups + H.sites_pbwt_grouping[anchor];
        for (int d = 0; d < H.depth; ++d) {
            unsigned long neighbor_idx = static_cast<unsigned long>(d) * addr_offset + base_idx;
            if (neighbor_idx < H.indexes_pbwt_neighbour.size()) {
                std::cout << " " << H.indexes_pbwt_neighbour[neighbor_idx];
            }
        }
        std::cout << std::endl;
    }
    int alt1_count=0, alt2_count=0, alt3_count=0, ref_count=0;
    for (unsigned int hidx : neigh) {
        if (hidx < 2) continue; // skip target's own haps
        AlleleClass cls = donor_class_from_packed(ctx.packed_codes, ss, hidx);
        if (trace) {
            std::cout << "  neigh hap=" << hidx << " class=" << (int)cls << std::endl;
        }
        if (cls == ALT1) alt1_count++;
        else if (cls == ALT2) alt2_count++;
        else if (cls == ALT3) alt3_count++;
        else ref_count++;
    }
    const int counted = alt1_count + alt2_count + alt3_count + ref_count;
    const double total = counted > 0 ? (double)counted : 1.0;
    return TargetResult{
        target_class,
        alt1_count / total,
        alt2_count / total,
        alt3_count / total,
        ref_count  / total
    };
}

} // namespace

int main() {
    TEST_INIT("test_pbwt_multiallelic_enrichment");
    std::cout << "Testing PBWT enrichment for supersite allele classes..." << std::endl;

    const double strict = 0.50; // require majority of neighbors to match target class
    const double loose = 0.30;  // REF case has lower signal because donors split across classes

    TargetResult t1 = run_case(ALT1);
    std::cout << "ALT1 target neighbor fractions: ALT1=" << t1.alt1_frac
              << " ALT2=" << t1.alt2_frac << " ALT3=" << t1.alt3_frac
              << " REF=" << t1.ref_frac << std::endl;
    assert(t1.alt1_frac >= strict && t1.alt1_frac > t1.alt2_frac && t1.alt1_frac > t1.ref_frac);

    TargetResult t2 = run_case(ALT2);
    std::cout << "ALT2 target neighbor fractions: ALT1=" << t2.alt1_frac
              << " ALT2=" << t2.alt2_frac << " ALT3=" << t2.alt3_frac
              << " REF=" << t2.ref_frac << std::endl;
    assert(t2.alt2_frac >= strict && t2.alt2_frac > t2.alt1_frac && t2.alt2_frac > t2.ref_frac);

    TargetResult t3 = run_case(REF);
    std::cout << "REF target neighbor fractions: ALT1=" << t3.alt1_frac
              << " ALT2=" << t3.alt2_frac << " ALT3=" << t3.alt3_frac
              << " REF=" << t3.ref_frac << std::endl;
    assert(t3.ref_frac >= loose && t3.ref_frac >= t3.alt1_frac && t3.ref_frac >= t3.alt2_frac);

    std::cout << "✓ SUCCESS: PBWT enriches for matching supersite allele class" << std::endl;
    TEST_SUMMARY();
    return 0;
}

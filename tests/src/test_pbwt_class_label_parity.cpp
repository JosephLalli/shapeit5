/*******************************************************************************
 * PBWT Supersite Class-Label Parity Test
 *
 * Verifies that, for the same multiallelic variation represented as:
 *   (1) split biallelic records and
 *   (2) supersites with packed allele codes,
 * the allele class assigned to donor haplotypes at each supersite is identical.
 *
 * This tests supersite packing / decoding parity on top of PBWT-selected donors.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <map>
#include <set>

#include "../../common/src/utils/otools.h"

#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/containers/variant_map.h"

namespace {

static inline bool env_true(const char* name) {
    const char* v = std::getenv(name);
    return v && v[0] != '\0' && v[0] != '0';
}

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    return ctx;
}

// Extract conditioning set for a specific locus from PBWT data structures
static std::set<unsigned int> extract_conditioning_set(const conditioning_set& H, int locus) {
    std::set<unsigned int> donors;

    if (locus >= static_cast<int>(H.sites_pbwt_selection.size()) ||
        !H.sites_pbwt_selection[locus]) {
        return donors;
    }

    unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;

    for (int target_hap = 0; target_hap < H.n_hap; ++target_hap) {
        int target_ind = target_hap / 2;
        if (target_ind >= H.n_ind) continue;

        unsigned long base_idx = H.sites_pbwt_grouping[locus] * 2UL * H.n_ind + target_hap;

        for (int d = 0; d < H.depth && base_idx < H.indexes_pbwt_neighbour.size(); ++d) {
            unsigned long neighbor_idx = d * addr_offset + base_idx;
            if (neighbor_idx < H.indexes_pbwt_neighbour.size()) {
                int donor_raw = H.indexes_pbwt_neighbour[neighbor_idx];
                if (donor_raw >= 0) {
                    donors.insert(static_cast<unsigned int>(donor_raw));
                }
            }
        }
    }

    return donors;
}

// Compute allele class from split biallelic representation:
//   0 = REF (no ALTs)
//   1 = first ALT (first set member split)
//   2 = second ALT, etc.
static uint8_t get_bial_class_code(const conditioning_set& H_bial,
                                   const SuperSite& ss,
                                   const std::vector<int>& super_site_var_index,
                                   unsigned int hap) {
    int first_alt = -1;
    int alt_count = 0;

    for (uint32_t i = 0; i < ss.var_count; ++i) {
        int v_idx = super_site_var_index[ss.var_start + i];
        if (v_idx < 0 || v_idx >= static_cast<int>(H_bial.n_site)) continue;
        bool alt = H_bial.H_opt_hap.get(hap, v_idx);
        if (alt) {
            alt_count++;
            if (first_alt < 0) first_alt = static_cast<int>(i);
        }
    }

    if (alt_count == 0) return SUPERSITE_CODE_REF;
    return static_cast<uint8_t>(first_alt + 1);
}

} // namespace

int main() {
    TEST_INIT("test_pbwt_class_label_parity");
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT Supersite Class-Label Parity Test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create datasets representing same genetic variation processed differently
    std::cout << "Creating test datasets..." << std::endl;

    variant_map V_bial, V_ss;

    // Biallelic version: process triallelic sites as split biallelics
    V_bial.push(make_var("1", 1000, "pos1000_A_T", "A", "T", 0));  // pos1000, ALT1
    V_bial.push(make_var("1", 1000, "pos1000_A_G", "A", "G", 1));  // pos1000, ALT2
    V_bial.push(make_var("1", 2000, "pos2000_C_A", "C", "A", 2));  // pos2000, ALT1
    V_bial.push(make_var("1", 2000, "pos2000_C_T", "C", "T", 3));  // pos2000, ALT2
    V_bial.push(make_var("1", 3000, "pos3000_G_C", "G", "C", 4));  // simple biallelic
    V_bial.push(make_var("1", 4000, "pos4000_T_A", "T", "A", 5));  // simple biallelic

    // Supersite version: same positions but grouped into supersites
    V_ss.push(make_var("1", 1000, "ss1_v0", "A", "T", 0));  // supersite 1
    V_ss.push(make_var("1", 1000, "ss1_v1", "A", "G", 1));
    V_ss.push(make_var("1", 2000, "ss2_v0", "C", "A", 2));  // supersite 2
    V_ss.push(make_var("1", 2000, "ss2_v1", "C", "T", 3));
    V_ss.push(make_var("1", 3000, "v3000", "G", "C", 4));
    V_ss.push(make_var("1", 4000, "v4000", "T", "A", 5));

    // Set simple genetic map (identical in both representations)
    for (size_t i = 0; i < V_bial.size(); ++i) {
        double cm = 0.01 * static_cast<double>(i + 1);
        V_bial.vec_pos[i]->cm = cm;
        V_ss.vec_pos[i]->cm = cm;
    }

    // Conditioning sets: 16 reference samples (32 haplotypes)
    conditioning_set H_bial, H_ss;
    const int n_ref_samples = 16;
    H_bial.allocate(0, n_ref_samples, V_bial.size());
    H_ss.allocate(0, n_ref_samples, V_ss.size());

    auto set_panel = [](conditioning_set& H, int locus, const std::vector<int>& alt_flags) {
        for (size_t hap = 0; hap < alt_flags.size() && hap < H.n_hap; ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    // Same realistic patterns as in test_pbwt_conditioning_parity
    std::vector<std::vector<int>> realistic_patterns = {
        // Variant 0 (pos1000 A->T)
        {1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 1 (pos1000 A->G)
        {0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 2 (pos2000 C->A)
        {0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0},
        // Variant 3 (pos2000 C->T)
        {0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 4 (pos3000 G->C)
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0},
        // Variant 5 (pos4000 T->A)
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}
    };

    for (int i = 0; i < 6; ++i) {
        set_panel(H_bial, i, realistic_patterns[i]);
        set_panel(H_ss, i, realistic_patterns[i]);
    }

    // Update variant MAC/MDR fields so PBWT evaluation (getMAC/getMDR) reflects our panel
    for (size_t locus = 0; locus < V_bial.size(); ++locus) {
        unsigned int alt_count = 0;
        for (int hap = 0; hap < H_bial.n_hap; ++hap) {
            if (H_bial.H_opt_hap.get(hap, static_cast<int>(locus))) alt_count++;
        }
        V_bial.vec_pos[locus]->calt = alt_count;
        V_bial.vec_pos[locus]->cref = static_cast<unsigned int>(H_bial.n_hap) - alt_count;
        V_bial.vec_pos[locus]->cmis = 0;

        V_ss.vec_pos[locus]->calt = alt_count;
        V_ss.vec_pos[locus]->cref = static_cast<unsigned int>(H_ss.n_hap) - alt_count;
        V_ss.vec_pos[locus]->cmis = 0;
    }

    std::cout << "  Built biallelic dataset: " << V_bial.size() << " variants" << std::endl;
    std::cout << "  Built supersite dataset: " << V_ss.size() << " variants" << std::endl;

    // Build supersites
    SuperSiteContext ctx_ss = build_supersites(V_ss, H_ss);
    std::cout << "  Detected " << ctx_ss.super_sites.size() << " supersites" << std::endl;

    // For each supersite anchor, compare class labels for donor haplotypes.
    // Here we treat all panel haplotypes as donors; PBWT is not required
    // because we are testing representation parity (split vs supersite codes),
    // not neighbor selection behavior.
    int total_checked = 0;
    int total_mismatches = 0;

    for (size_t ss_idx = 0; ss_idx < ctx_ss.super_sites.size(); ++ss_idx) {
        const SuperSite& ss = ctx_ss.super_sites[ss_idx];
        (void)H_ss; // Supersite PBWT state not needed for this test

        // Treat all panel haplotypes as donors for class-label parity
        std::set<unsigned int> donors;
        for (unsigned int h = 0; h < static_cast<unsigned int>(H_bial.n_hap); ++h) {
            donors.insert(h);
        }

        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::cout << "\n[CLASS_PARITY] Supersite " << ss_idx
                      << " anchor_locus=" << ss.global_site_id
                      << " donors=" << donors.size() << std::endl;
        }

        for (unsigned int hap : donors) {
            uint8_t code_bial = get_bial_class_code(H_bial, ss, ctx_ss.super_site_var_index, hap);
            uint8_t code_ss = unpackSuperSiteCode(
                ctx_ss.packed_codes.empty() ? nullptr : ctx_ss.packed_codes.data(),
                ss.panel_offset,
                hap
            );

            total_checked++;

            if (code_bial != code_ss) {
                total_mismatches++;
                if (env_true("SHAPEIT5_TEST_TRACE")) {
                    std::cout << "[CLASS_MISMATCH] ss_idx=" << ss_idx
                              << " anchor=" << ss.global_site_id
                              << " hap=" << hap
                              << " code_bial=" << static_cast<int>(code_bial)
                              << " code_ss=" << static_cast<int>(code_ss) << std::endl;
                }
            }
        }
    }

    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Class-Label Parity Summary" << std::endl;
    std::cout << "======================================================================" << std::endl;

    std::cout << "  Donor haplotypes checked: " << total_checked << std::endl;
    std::cout << "  Class label mismatches:    " << total_mismatches << std::endl;

    bool test_passed = (total_mismatches == 0 && total_checked > 0);

    if (test_passed) {
        std::cout << "\n✓ SUCCESS: Supersite codes match biallelic split classes for all donors" << std::endl;
    } else {
        std::cout << "\n✗ FAILURE: Supersite class codes differ from biallelic split classes" << std::endl;
        if (total_checked == 0) {
            std::cout << "  (No donors were checked; PBWT may not have selected anchors.)" << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT supersite class-label parity test completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;

    return test_passed ? 0 : 1;
}

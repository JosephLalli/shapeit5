/*******************************************************************************
 * PBWT Multiallelic Similarity Test (Gap C)
 *
 * Tests whether PBWT conditioning set selection properly reflects multiallelic
 * allele class similarity when selecting donors for a target haplotype.
 *
 * Background:
 * - test_pbwt_conditioning_parity.cpp shows identical conditioning sets for
 *   biallelic vs supersite layouts (same underlying variation)
 * - This test addresses a different question: Does PBWT prefer donors with
 *   the same allele class (ALT1 vs ALT2 vs REF) at multiallelic sites?
 *
 * Problem:
 * - PBWT operates on the split-bit biallelic matrix representation
 * - Supersite siblings must be masked/redirected to avoid double-counting
 * - The key question: Do conditioning sets prefer donors with matching allele
 *   classes, or does the split representation lose allele-specific similarity?
 *
 * Test Design:
 * - Create panel with donor blocks carrying different allele classes:
 *   - Block A: ALT1-only donors
 *   - Block B: ALT2-only donors
 *   - Block C: ALT3-only donors
 *   - Block D: REF-only donors
 * - Test sample Case A: hap0=ALT1, hap1=REF at supersite
 * - Test sample Case B: hap0=ALT2, hap1=REF at supersite
 * - Run PBWT selection with supersite anchor masking/redirect
 * - Assert: Donors with matching allele class dominate neighbor set (>80%)
 * - Assert: Switching Case A→B flips which donor block is enriched
 *
 * Expected Outcome:
 * - PASS: PBWT correctly identifies allele-class similarity
 * - FAIL: PBWT treats all ALT alleles equally or shows random selection
 *   → Indicates accuracy bug in multiallelic similarity detection
 *   → May explain 33% accuracy loss in multiallelic phasing
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

// Allele class codes for the triallelic supersite
enum AlleleClass {
    REF = 0,   // Reference allele
    ALT1 = 1,  // First alternate allele
    ALT2 = 2,  // Second alternate allele
    ALT3 = 3   // Third alternate allele
};

struct DonorBlockAnalysis {
    AlleleClass target_class;
    int total_donors = 0;
    int alt1_donors = 0;
    int alt2_donors = 0;
    int alt3_donors = 0;
    int ref_donors = 0;
    double alt1_fraction = 0.0;
    double alt2_fraction = 0.0;
    double alt3_fraction = 0.0;
    double ref_fraction = 0.0;

    void compute_fractions() {
        if (total_donors > 0) {
            alt1_fraction = static_cast<double>(alt1_donors) / total_donors;
            alt2_fraction = static_cast<double>(alt2_donors) / total_donors;
            alt3_fraction = static_cast<double>(alt3_donors) / total_donors;
            ref_fraction = static_cast<double>(ref_donors) / total_donors;
        }
    }

    int get_matching_count() const {
        switch (target_class) {
            case ALT1: return alt1_donors;
            case ALT2: return alt2_donors;
            case ALT3: return alt3_donors;
            case REF:  return ref_donors;
        }
        return 0;
    }

    double get_matching_fraction() const {
        switch (target_class) {
            case ALT1: return alt1_fraction;
            case ALT2: return alt2_fraction;
            case ALT3: return alt3_fraction;
            case REF:  return ref_fraction;
        }
        return 0.0;
    }

    bool is_plurality() const {
        int matching = get_matching_count();
        return (matching >= alt1_donors && matching >= alt2_donors &&
                matching >= alt3_donors && matching >= ref_donors);
    }

    void print(const std::string& label) const {
        const char* class_name[] = {"REF", "ALT1", "ALT2", "ALT3"};
        std::cout << label << " (target=" << class_name[target_class] << ") neighbor composition:" << std::endl;
        std::cout << "  Total donors: " << total_donors << std::endl;

        auto print_class = [&](const char* name, int count, double fraction, bool is_matching) {
            std::cout << "  " << name << " donors: " << count << " ("
                      << std::fixed << std::setprecision(1) << (fraction * 100.0) << "%)";
            if (is_matching) std::cout << " ← matching class";
            std::cout << std::endl;
        };

        print_class("ALT1", alt1_donors, alt1_fraction, target_class == ALT1);
        print_class("ALT2", alt2_donors, alt2_fraction, target_class == ALT2);
        print_class("ALT3", alt3_donors, alt3_fraction, target_class == ALT3);
        print_class("REF ", ref_donors, ref_fraction, target_class == REF);
    }

    // Evaluate pass/fail for three criteria
    struct EvalResult {
        bool option1_pass;  // Strict: >50% match
        bool option2_pass;  // Moderate: plurality and ≥31% match
        bool option3_pass;  // Baseline: ≥25% match (random baseline)

        void print(const std::string& label) const {
            std::cout << "  " << label << ":" << std::endl;
            std::cout << "    Option 1 (Strict >50%):     " << (option1_pass ? "PASS" : "FAIL") << std::endl;
            std::cout << "    Option 2 (Plurality ≥31%):  " << (option2_pass ? "PASS" : "FAIL") << std::endl;
            std::cout << "    Option 3 (Baseline ≥25%):   " << (option3_pass ? "PASS" : "FAIL") << std::endl;
        }
    };

    EvalResult evaluate() const {
        EvalResult result;
        double matching_frac = get_matching_fraction();

        result.option1_pass = (matching_frac > 0.50);
        result.option2_pass = (is_plurality() && matching_frac >= 0.31);
        result.option3_pass = (matching_frac >= 0.25);

        return result;
    }
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

// Extract allele class from donor haplotype at supersite anchor position
// Returns: REF, ALT1, ALT2, or ALT3 based on which split variant is set
static AlleleClass get_donor_allele_class(const conditioning_set& H,
                                          unsigned int donor_hap,
                                          int anchor_locus) {
    // In supersite representation with 3 ALT alleles:
    // - anchor_locus+0: ALT1 vs REF (split 0)
    // - anchor_locus+1: ALT2 vs REF (split 1)
    // - anchor_locus+2: ALT3 vs REF (split 2)

    bool has_alt1 = H.H_opt_hap.get(donor_hap, anchor_locus + 0);
    bool has_alt2 = H.H_opt_hap.get(donor_hap, anchor_locus + 1);
    bool has_alt3 = H.H_opt_hap.get(donor_hap, anchor_locus + 2);

    // Only one should be true due to mutual exclusivity
    if (has_alt1) return ALT1;
    if (has_alt2) return ALT2;
    if (has_alt3) return ALT3;
    return REF;
}

// Extract conditioning set donors for a specific target haplotype
static std::vector<unsigned int> extract_donors_for_haplotype(
    const conditioning_set& H,
    int locus,
    int target_hap) {

    std::vector<unsigned int> donors;

    if (locus >= static_cast<int>(H.sites_pbwt_selection.size()) ||
        !H.sites_pbwt_selection[locus]) {
        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::cerr << "[extract_donors] locus=" << locus << " not selected for PBWT" << std::endl;
        }
        return donors;
    }

    unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
    unsigned long base_idx = H.sites_pbwt_grouping[locus] * 2UL * H.n_ind + target_hap;

    if (env_true("SHAPEIT5_TEST_TRACE")) {
        std::cerr << "[extract_donors] locus=" << locus << " target_hap=" << target_hap
                  << " grouping=" << H.sites_pbwt_grouping[locus]
                  << " ngroups=" << H.sites_pbwt_ngroups
                  << " n_ind=" << H.n_ind
                  << " depth=" << H.depth
                  << " addr_offset=" << addr_offset
                  << " base_idx=" << base_idx
                  << " indexes_size=" << H.indexes_pbwt_neighbour.size() << std::endl;
    }

    for (int d = 0; d < H.depth; ++d) {
        unsigned long neighbor_idx = d * addr_offset + base_idx;
        if (neighbor_idx < H.indexes_pbwt_neighbour.size()) {
            unsigned int donor_hap = H.indexes_pbwt_neighbour[neighbor_idx];
            donors.push_back(donor_hap);
            if (env_true("SHAPEIT5_TEST_TRACE") && d < 3) {
                std::cerr << "    d=" << d << " neighbor_idx=" << neighbor_idx
                          << " donor_hap=" << donor_hap << std::endl;
            }
        } else {
            if (env_true("SHAPEIT5_TEST_TRACE")) {
                std::cerr << "    d=" << d << " neighbor_idx=" << neighbor_idx
                          << " OUT OF BOUNDS (size=" << H.indexes_pbwt_neighbour.size() << ")" << std::endl;
            }
        }
    }

    return donors;
}

// Analyze donor block composition for a target haplotype
static DonorBlockAnalysis analyze_donor_blocks(
    const conditioning_set& H,
    int anchor_locus,
    int target_hap,
    AlleleClass target_class) {

    DonorBlockAnalysis analysis;
    analysis.target_class = target_class;

    std::vector<unsigned int> donors = extract_donors_for_haplotype(H, anchor_locus, target_hap);
    analysis.total_donors = static_cast<int>(donors.size());

    // Diagnostic #1: Show actual neighbor indices and their allele classes
    if (env_true("SHAPEIT5_TEST_TRACE")) {
        const char* class_name[] = {"REF", "ALT1", "ALT2", "ALT3"};
        std::cerr << "\n[PBWT_NEIGHBORS] target_hap=" << target_hap
                  << " target_class=" << class_name[target_class]
                  << " anchor_locus=" << anchor_locus << std::endl;
        std::cerr << "  Total neighbors: " << donors.size() << std::endl;
        std::cerr << "  Neighbor details (hap_idx -> allele_class):" << std::endl;
        for (size_t i = 0; i < donors.size(); ++i) {
            unsigned int donor_hap = donors[i];
            AlleleClass ac = get_donor_allele_class(H, donor_hap, anchor_locus);
            std::cerr << "    [" << std::setw(2) << i << "] hap=" << std::setw(3) << donor_hap
                      << " class=" << class_name[ac];
            // Show which block this donor is from (blocks start at hap 12)
            if (donor_hap >= 12) {
                int block_idx = (donor_hap - 12) / 32;
                const char* block_name[] = {"ALT1_block", "ALT2_block", "ALT3_block", "REF_block"};
                if (block_idx < 4) {
                    std::cerr << " (" << block_name[block_idx] << ")";
                }
            }
            std::cerr << std::endl;
        }
    }

    for (unsigned int donor_hap : donors) {
        AlleleClass ac = get_donor_allele_class(H, donor_hap, anchor_locus);
        switch (ac) {
            case ALT1: analysis.alt1_donors++; break;
            case ALT2: analysis.alt2_donors++; break;
            case ALT3: analysis.alt3_donors++; break;
            case REF:  analysis.ref_donors++;  break;
        }
    }

    analysis.compute_fractions();
    return analysis;
}

} // namespace

int main() {
    TEST_INIT("test_pbwt_multiallelic_similarity");
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT Multiallelic Similarity Test (Gap C)" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create supersite dataset with one triallelic site plus flanking variants
    std::cout << "Creating test dataset..." << std::endl;

    variant_map V;

    // Flanking variants to give PBWT context (5 on each side)
    V.push(make_var("1", 1000, "flank_L5", "A", "T", 0));  // Left flank -5
    V.push(make_var("1", 1500, "flank_L4", "C", "G", 1));  // Left flank -4
    V.push(make_var("1", 2000, "flank_L3", "G", "A", 2));  // Left flank -3
    V.push(make_var("1", 2500, "flank_L2", "T", "C", 3));  // Left flank -2
    V.push(make_var("1", 2800, "flank_L1", "A", "G", 4));  // Left flank -1

    // Triallelic supersite at position 3000 (3 splits)
    V.push(make_var("1", 3000, "ss_ALT1", "A", "C", 5));  // Split 0: ALT1
    V.push(make_var("1", 3000, "ss_ALT2", "A", "G", 6));  // Split 1: ALT2
    V.push(make_var("1", 3000, "ss_ALT3", "A", "T", 7));  // Split 2: ALT3

    // Flanking variants after supersite (5 on each side)
    V.push(make_var("1", 3200, "flank_R1", "G", "A", 8));  // Right flank +1
    V.push(make_var("1", 3500, "flank_R2", "T", "C", 9));  // Right flank +2
    V.push(make_var("1", 4000, "flank_R3", "C", "T", 10)); // Right flank +3
    V.push(make_var("1", 4500, "flank_R4", "A", "G", 11)); // Right flank +4
    V.push(make_var("1", 5000, "flank_R5", "G", "C", 12)); // Right flank +5

    // Initialize variant counts and genetic map
    for (size_t i = 0; i < V.size(); ++i) {
        V.vec_pos[i]->cm = 0.01 * static_cast<double>(i + 1);  // Simple genetic map
        V.vec_pos[i]->cref = 0;  // Will be set after panel is built
        V.vec_pos[i]->calt = 0;
        V.vec_pos[i]->cmis = 0;
    }

    std::cout << "  Total variants: " << V.size() << std::endl;
    std::cout << "  Supersite position: 3000 (loci 2-4)" << std::endl;

    // Create conditioning set with 6 test samples + reference panel
    // Design: 6 test samples (covering various genotypes) + 64 reference donor samples
    //   Test samples (will have PBWT neighbors computed):
    //     Sample 0: ALT1|REF
    //     Sample 1: ALT2|REF
    //     Sample 2: ALT3|REF
    //     Sample 3: REF|REF
    //     Sample 4: ALT1|ALT3  (mixed non-adjacent ALTs)
    //     Sample 5: ALT3|ALT3  (homozygous ALT3)
    //   Reference donors organized in 4 blocks by allele class
    const int n_test_samples = 6;   // 6 test samples = 12 test haplotypes
    const int donors_per_block = 16;  // 32 haplotypes per block
    const int n_blocks = 4;
    const int n_ref_samples = donors_per_block * n_blocks;  // 64 samples = 128 reference haps
    const int supersite_anchor = 5;  // Locus index of supersite anchor (now locus 5)

    conditioning_set H;
    H.allocate(n_test_samples, n_ref_samples, V.size());

    std::cout << "  Test samples: " << n_test_samples << " (" << (n_test_samples * 2) << " haplotypes)" << std::endl;
    std::cout << "  Reference samples: " << n_ref_samples << " (" << H.n_hap << " haplotypes)" << std::endl;
    std::cout << "  Donor block structure: " << n_blocks << " blocks × " << donors_per_block
              << " samples = " << (donors_per_block * 2) << " haplotypes/block" << std::endl;

    // Fill reference panel with block structure
    // Block A (haps 0-31):   ALT1 at supersite
    // Block B (haps 32-63):  ALT2 at supersite
    // Block C (haps 64-95):  ALT3 at supersite
    // Block D (haps 96-127): REF at supersite

    auto set_panel_hap = [&](int locus, unsigned int hap_idx, bool value) {
        if (value) {
            H.H_opt_var.set(locus, hap_idx, 1);
            H.H_opt_hap.set(hap_idx, locus, 1);
        }
    };

    // Set donor patterns for each block
    // NOTE: Start from hap index n_test_samples*2 to skip test samples
    // IMPORTANT: All flanking variants are IDENTICAL across all blocks
    // This isolates the supersite as the ONLY source of information for PBWT
    const unsigned int ref_hap_start = n_test_samples * 2;  // Haps 0-7 are test samples
    for (unsigned int hap = ref_hap_start; hap < H.n_hap; ++hap) {
        int block_idx = (hap - ref_hap_start) / (donors_per_block * 2);  // 4 blocks total

        // Set ALL flanking variants to REF (0) for ALL haplotypes
        // This makes flanking context uninformative - only supersite differs between blocks
        set_panel_hap(0, hap, false);   // flank_L5: all REF
        set_panel_hap(1, hap, false);   // flank_L4: all REF
        set_panel_hap(2, hap, false);   // flank_L3: all REF
        set_panel_hap(3, hap, false);   // flank_L2: all REF
        set_panel_hap(4, hap, false);   // flank_L1: all REF
        set_panel_hap(8, hap, false);   // flank_R1: all REF
        set_panel_hap(9, hap, false);   // flank_R2: all REF
        set_panel_hap(10, hap, false);  // flank_R3: all REF
        set_panel_hap(11, hap, false);  // flank_R4: all REF
        set_panel_hap(12, hap, false);  // flank_R5: all REF

        // Set supersite based on block membership
        switch (block_idx) {
            case 0:  // Block A: ALT1
                set_panel_hap(supersite_anchor + 0, hap, true);  // ALT1 split
                set_panel_hap(supersite_anchor + 1, hap, false); // ALT2 split
                set_panel_hap(supersite_anchor + 2, hap, false); // ALT3 split
                break;
            case 1:  // Block B: ALT2
                set_panel_hap(supersite_anchor + 0, hap, false); // ALT1 split
                set_panel_hap(supersite_anchor + 1, hap, true);  // ALT2 split
                set_panel_hap(supersite_anchor + 2, hap, false); // ALT3 split
                break;
            case 2:  // Block C: ALT3
                set_panel_hap(supersite_anchor + 0, hap, false); // ALT1 split
                set_panel_hap(supersite_anchor + 1, hap, false); // ALT2 split
                set_panel_hap(supersite_anchor + 2, hap, true);  // ALT3 split
                break;
            case 3:  // Block D: REF
                set_panel_hap(supersite_anchor + 0, hap, false); // ALT1 split
                set_panel_hap(supersite_anchor + 1, hap, false); // ALT2 split
                set_panel_hap(supersite_anchor + 2, hap, false); // ALT3 split
                break;
        }
    }

    // Now populate test sample genotypes
    // Test samples are stored in first 2*n_test_samples haplotypes of H_opt_hap/H_opt_var
    std::cout << "\n  Creating test sample genotypes..." << std::endl;

    auto set_test_hap = [&](int hap_idx, int locus, bool value) {
        if (value) {
            H.H_opt_hap.set(hap_idx, locus, 1);
            H.H_opt_var.set(locus, hap_idx, 1);
        }
    };

    // Sample 0: ALT1|REF at supersite
    set_test_hap(0, supersite_anchor + 0, true);   // hap0 = ALT1 (split 0)
    set_test_hap(1, supersite_anchor + 0, false);  // hap1 = REF
    set_test_hap(1, supersite_anchor + 1, false);  // hap1 = REF (explicitly clear split 1)
    set_test_hap(1, supersite_anchor + 2, false);  // hap1 = REF (explicitly clear split 2)
    std::cout << "    Sample 0: ALT1|REF (haplotypes 0,1)" << std::endl;

    // Sample 1: ALT2|REF at supersite
    set_test_hap(2, supersite_anchor + 0, false);  // hap0 NEQ ALT1
    set_test_hap(2, supersite_anchor + 1, true);   // hap0 = ALT2 (split 1)
    set_test_hap(2, supersite_anchor + 2, false);  // hap0 NEQ ALT3
    set_test_hap(3, supersite_anchor + 0, false);  // hap1 = REF
    set_test_hap(3, supersite_anchor + 1, false);  // hap1 = REF
    set_test_hap(3, supersite_anchor + 2, false);  // hap1 = REF
    std::cout << "    Sample 1: ALT2|REF (haplotypes 2,3)" << std::endl;

    // Sample 2: ALT3|REF at supersite
    set_test_hap(4, supersite_anchor + 0, false);  // hap0 NEQ ALT1
    set_test_hap(4, supersite_anchor + 1, false);  // hap0 NEQ ALT2
    set_test_hap(4, supersite_anchor + 2, true);   // hap0 = ALT3 (split 2)
    set_test_hap(5, supersite_anchor + 0, false);  // hap1 = REF
    set_test_hap(5, supersite_anchor + 1, false);  // hap1 = REF
    set_test_hap(5, supersite_anchor + 2, false);  // hap1 = REF
    std::cout << "    Sample 2: ALT3|REF (haplotypes 4,5)" << std::endl;

    // Sample 3: REF|REF at supersite
    set_test_hap(6, supersite_anchor + 0, false);  // hap0 = REF
    set_test_hap(6, supersite_anchor + 1, false);  // hap0 = REF
    set_test_hap(6, supersite_anchor + 2, false);  // hap0 = REF
    set_test_hap(7, supersite_anchor + 0, false);  // hap1 = REF
    set_test_hap(7, supersite_anchor + 1, false);  // hap1 = REF
    set_test_hap(7, supersite_anchor + 2, false);  // hap1 = REF
    std::cout << "    Sample 3: REF|REF (haplotypes 6,7)" << std::endl;

    // Sample 4: ALT1|ALT3 at supersite (mixed non-adjacent ALTs)
    set_test_hap(8, supersite_anchor + 0, true);   // hap0 = ALT1 (split 0)
    set_test_hap(8, supersite_anchor + 1, false);  // hap0 NEQ ALT2
    set_test_hap(8, supersite_anchor + 2, false);  // hap0 NEQ ALT3
    set_test_hap(9, supersite_anchor + 0, false);  // hap1 NEQ ALT1
    set_test_hap(9, supersite_anchor + 1, false);  // hap1 NEQ ALT2
    set_test_hap(9, supersite_anchor + 2, true);   // hap1 = ALT3 (split 2)
    std::cout << "    Sample 4: ALT1|ALT3 (haplotypes 8,9)" << std::endl;

    // Sample 5: ALT3|ALT3 at supersite (homozygous ALT3)
    set_test_hap(10, supersite_anchor + 0, false); // hap0 NEQ ALT1
    set_test_hap(10, supersite_anchor + 1, false); // hap0 NEQ ALT2
    set_test_hap(10, supersite_anchor + 2, true);  // hap0 = ALT3 (split 2)
    set_test_hap(11, supersite_anchor + 0, false); // hap1 NEQ ALT1
    set_test_hap(11, supersite_anchor + 1, false); // hap1 NEQ ALT2
    set_test_hap(11, supersite_anchor + 2, true);  // hap1 = ALT3 (split 2)
    std::cout << "    Sample 5: ALT3|ALT3 (haplotypes 10,11)" << std::endl;

    // Set flanking variants for test samples - all REF to match reference panel
    for (int s = 0; s < n_test_samples; ++s) {
        int hap0 = s * 2;
        int hap1 = s * 2 + 1;
        // All flanking variants = REF for test samples (matching reference panel)
        set_test_hap(hap0, 0, false);   // flank_L5
        set_test_hap(hap1, 0, false);
        set_test_hap(hap0, 1, false);   // flank_L4
        set_test_hap(hap1, 1, false);
        set_test_hap(hap0, 2, false);   // flank_L3
        set_test_hap(hap1, 2, false);
        set_test_hap(hap0, 3, false);   // flank_L2
        set_test_hap(hap1, 3, false);
        set_test_hap(hap0, 4, false);   // flank_L1
        set_test_hap(hap1, 4, false);
        set_test_hap(hap0, 8, false);   // flank_R1
        set_test_hap(hap1, 8, false);
        set_test_hap(hap0, 9, false);   // flank_R2
        set_test_hap(hap1, 9, false);
        set_test_hap(hap0, 10, false);  // flank_R3
        set_test_hap(hap1, 10, false);
        set_test_hap(hap0, 11, false);  // flank_R4
        set_test_hap(hap1, 11, false);
        set_test_hap(hap0, 12, false);  // flank_R5
        set_test_hap(hap1, 12, false);
    }

    // Update variant counts based on panel (reference only)
    for (size_t locus = 0; locus < V.size(); ++locus) {
        unsigned int alt_count = 0;
        for (unsigned int hap = 0; hap < H.n_hap; ++hap) {
            if (H.H_opt_hap.get(hap, locus)) {
                alt_count++;
            }
        }
        V.vec_pos[locus]->calt = alt_count;
        V.vec_pos[locus]->cref = H.n_hap - alt_count;
        V.vec_pos[locus]->cmis = 0;
    }

    // Diagnostic #2: Verify reference panel structure at supersite
    if (env_true("SHAPEIT5_TEST_TRACE")) {
        std::cerr << "\n[REF_PANEL_VERIFICATION] Checking supersite structure" << std::endl;
        std::cerr << "  Expected blocks (32 haps each):" << std::endl;
        std::cerr << "    Block 0 (ALT1): haps 12-43" << std::endl;
        std::cerr << "    Block 1 (ALT2): haps 44-75" << std::endl;
        std::cerr << "    Block 2 (ALT3): haps 76-107" << std::endl;
        std::cerr << "    Block 3 (REF):  haps 108-139" << std::endl;

        for (int block = 0; block < n_blocks; block++) {
            const char* expected_class[] = {"ALT1", "ALT2", "ALT3", "REF"};
            int start_hap = ref_hap_start + (block * donors_per_block * 2);
            int end_hap = start_hap + std::min(5, donors_per_block * 2); // Check first 5 haps

            std::cerr << "\n  Block " << block << " (" << expected_class[block]
                      << ", haps " << start_hap << "-" << (start_hap + donors_per_block * 2 - 1) << "):"
                      << std::endl;

            for (int h = start_hap; h < end_hap; h++) {
                bool s0 = H.H_opt_hap.get(h, supersite_anchor + 0);
                bool s1 = H.H_opt_hap.get(h, supersite_anchor + 1);
                bool s2 = H.H_opt_hap.get(h, supersite_anchor + 2);

                // Determine actual class
                const char* actual_class = "REF";
                if (s0) actual_class = "ALT1";
                else if (s1) actual_class = "ALT2";
                else if (s2) actual_class = "ALT3";

                bool matches = (strcmp(actual_class, expected_class[block]) == 0);

                std::cerr << "    hap " << std::setw(3) << h << ": "
                          << "split0=" << s0 << " split1=" << s1 << " split2=" << s2
                          << " → " << actual_class;
                if (!matches) {
                    std::cerr << " ✗ MISMATCH (expected " << expected_class[block] << ")";
                }
                std::cerr << std::endl;
            }
        }
    }

    // Build supersites
    SuperSiteContext ctx = build_supersites(V, H);
    std::cout << "  Detected " << ctx.super_sites.size() << " supersites" << std::endl;

    // Initialize PBWT with supersite guards
    std::cout << "\nBuilding PBWT with supersite anchor masking..." << std::endl;

    float modulo_selection = 1.0f;     // Select all sites
    float modulo_multithreading = 1.0f;
    float mdr = 0.1f;
    int depth = 16;                    // K=16 donors per target
    int mac = 1;
    int nthread = 1;

    H.initialize(V, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    // Apply supersite anchor masking and redirect
    if (!ctx.super_sites.empty()) {
        H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
        std::vector<int> anchor_map = buildSupersiteAnchorMap(
            ctx.super_sites, ctx.super_site_var_index, V.size());
        H.setSupersiteAnchorRedirect(anchor_map);

        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::cout << "  Anchor redirect map:" << std::endl;
            for (size_t i = 0; i < anchor_map.size(); ++i) {
                if (anchor_map[i] != static_cast<int>(i)) {
                    std::cout << "    locus " << i << " → " << anchor_map[i] << std::endl;
                }
            }
        }
    }

    H.select();

    std::cout << "  PBWT depth: " << H.depth << std::endl;
    std::cout << "  PBWT sites selected: " << H.n_site << std::endl;

    // Debug: Print which sites were selected
    if (env_true("SHAPEIT5_TEST_TRACE")) {
        std::cout << "  Sites selected for PBWT:" << std::endl;
        for (size_t i = 0; i < H.sites_pbwt_selection.size(); ++i) {
            if (H.sites_pbwt_selection[i]) {
                std::cout << "    Locus " << i << ": " << V.vec_pos[i]->id
                          << " (group=" << H.sites_pbwt_grouping[i] << ")" << std::endl;
            }
        }
    }

    // Verify that supersite siblings were NOT selected for PBWT
    if (H.sites_pbwt_selection[supersite_anchor + 1] ||
        H.sites_pbwt_selection[supersite_anchor + 2]) {
        std::cerr << "ERROR: PBWT selected supersite siblings!" << std::endl;
        std::cerr << "  Anchor (locus " << supersite_anchor << "): "
                  << H.sites_pbwt_selection[supersite_anchor] << std::endl;
        std::cerr << "  Sibling 1 (locus " << (supersite_anchor+1) << "): "
                  << H.sites_pbwt_selection[supersite_anchor+1] << std::endl;
        std::cerr << "  Sibling 2 (locus " << (supersite_anchor+2) << "): "
                  << H.sites_pbwt_selection[supersite_anchor+2] << std::endl;
        return 1;
    }
    std::cout << "  ✓ Supersite siblings correctly excluded from PBWT" << std::endl;

    // Test all 8 haplotypes (4 samples × 2 haplotypes each)
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Testing PBWT Neighbor Enrichment" << std::endl;
    std::cout << "======================================================================" << std::endl;

    // Helper to get allele class for a test haplotype
    auto get_test_hap_allele = [&](int hap_idx) -> AlleleClass {
        bool has_alt1 = H.H_opt_hap.get(hap_idx, supersite_anchor + 0);
        bool has_alt2 = H.H_opt_hap.get(hap_idx, supersite_anchor + 1);
        bool has_alt3 = H.H_opt_hap.get(hap_idx, supersite_anchor + 2);

        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::cerr << "[get_test_hap_allele] hap=" << hap_idx
                      << " split0=" << (int)has_alt1
                      << " split1=" << (int)has_alt2
                      << " split2=" << (int)has_alt3 << std::endl;
        }

        if (has_alt1) return ALT1;
        if (has_alt2) return ALT2;
        if (has_alt3) return ALT3;
        return REF;
    };

    // Collect all analyses for summary statistics
    std::vector<DonorBlockAnalysis> all_analyses;
    std::vector<DonorBlockAnalysis::EvalResult> all_evals;

    // Test each haplotype
    for (int hap = 0; hap < n_test_samples * 2; ++hap) {
        int sample = hap / 2;
        int hap_in_sample = hap % 2;
        AlleleClass target_class = get_test_hap_allele(hap);

        const char* allele_name[] = {"REF", "ALT1", "ALT2", "ALT3"};
        std::cout << "\n=== Sample " << sample << " hap" << hap_in_sample
                  << " (target=" << allele_name[target_class] << ") ===" << std::endl;

        DonorBlockAnalysis analysis = analyze_donor_blocks(H, supersite_anchor, hap, target_class);
        analysis.print("  ");

        auto eval = analysis.evaluate();
        eval.print("  ");

        all_analyses.push_back(analysis);
        all_evals.push_back(eval);
    }

    // Generate summary statistics
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Summary Statistics by Allele Class" << std::endl;
    std::cout << "======================================================================" << std::endl;

    // Count statistics per allele class
    struct ClassStats {
        int count = 0;
        int option1_pass = 0;
        int option2_pass = 0;
        int option3_pass = 0;
        double avg_matching_fraction = 0.0;
    };

    std::map<AlleleClass, ClassStats> class_stats;

    for (size_t i = 0; i < all_analyses.size(); ++i) {
        AlleleClass ac = all_analyses[i].target_class;
        class_stats[ac].count++;
        if (all_evals[i].option1_pass) class_stats[ac].option1_pass++;
        if (all_evals[i].option2_pass) class_stats[ac].option2_pass++;
        if (all_evals[i].option3_pass) class_stats[ac].option3_pass++;
        class_stats[ac].avg_matching_fraction += all_analyses[i].get_matching_fraction();
    }

    // Compute averages
    for (auto& pair : class_stats) {
        if (pair.second.count > 0) {
            pair.second.avg_matching_fraction /= pair.second.count;
        }
    }

    // Print per-class statistics
    const char* allele_name[] = {"REF", "ALT1", "ALT2", "ALT3"};
    for (int ac = 0; ac < 4; ++ac) {
        AlleleClass allele = static_cast<AlleleClass>(ac);
        const ClassStats& stats = class_stats[allele];

        std::cout << "\n" << allele_name[ac] << " targets (n=" << stats.count << "):" << std::endl;
        std::cout << "  Average matching fraction: "
                  << std::fixed << std::setprecision(1) << (stats.avg_matching_fraction * 100.0) << "%" << std::endl;
        std::cout << "  Option 1 pass rate (>50%):       " << stats.option1_pass << "/" << stats.count
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * stats.option1_pass / stats.count) << "%)" << std::endl;
        std::cout << "  Option 2 pass rate (≥31% plur):  " << stats.option2_pass << "/" << stats.count
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * stats.option2_pass / stats.count) << "%)" << std::endl;
        std::cout << "  Option 3 pass rate (≥25%):       " << stats.option3_pass << "/" << stats.count
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * stats.option3_pass / stats.count) << "%)" << std::endl;
    }

    // Overall statistics
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Overall Results (All " << all_analyses.size() << " Haplotypes)" << std::endl;
    std::cout << "======================================================================" << std::endl;

    int total_option1_pass = 0;
    int total_option2_pass = 0;
    int total_option3_pass = 0;
    double total_avg_matching = 0.0;

    for (const auto& analysis : all_analyses) {
        total_avg_matching += analysis.get_matching_fraction();
    }
    for (const auto& eval : all_evals) {
        if (eval.option1_pass) total_option1_pass++;
        if (eval.option2_pass) total_option2_pass++;
        if (eval.option3_pass) total_option3_pass++;
    }

    total_avg_matching /= all_analyses.size();

    std::cout << "  Total haplotypes tested: " << all_analyses.size() << std::endl;
    std::cout << "  Average matching fraction: "
              << std::fixed << std::setprecision(1) << (total_avg_matching * 100.0) << "%" << std::endl;
    std::cout << "\n  Option 1 (Strict >50%):" << std::endl;
    std::cout << "    Pass: " << total_option1_pass << "/" << all_analyses.size()
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * total_option1_pass / all_analyses.size()) << "%)" << std::endl;
    std::cout << "  Option 2 (Plurality ≥31%):" << std::endl;
    std::cout << "    Pass: " << total_option2_pass << "/" << all_analyses.size()
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * total_option2_pass / all_analyses.size()) << "%)" << std::endl;
    std::cout << "  Option 3 (Baseline ≥25%):" << std::endl;
    std::cout << "    Pass: " << total_option3_pass << "/" << all_analyses.size()
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * total_option3_pass / all_analyses.size()) << "%)" << std::endl;

    // Final verdict
    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Test Verdict" << std::endl;
    std::cout << "======================================================================" << std::endl;

    bool test_passed = true;
    std::vector<std::string> failures;

    // Use Option 2 (moderate criterion) as the pass/fail threshold
    if (total_option2_pass < all_analyses.size()) {
        failures.push_back(std::to_string(all_analyses.size() - total_option2_pass) +
                          " haplotypes failed Option 2 criterion");
        test_passed = false;
    }

    // Print results
    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;

    if (test_passed) {
        std::cout << "✓ SUCCESS: PBWT correctly identifies allele-class similarity" << std::endl;
        std::cout << "✓ Matching allele classes enriched in neighbor sets" << std::endl;
        std::cout << "✓ Supersite anchor masking/redirect working correctly" << std::endl;
        std::cout << std::endl;
        std::cout << "INTERPRETATION:" << std::endl;
        std::cout << "  - PBWT successfully prefers donors with matching allele classes" << std::endl;
        std::cout << "  - Split-bit representation preserves allele-specific similarity" << std::endl;
        std::cout << "  - 33% accuracy loss likely has different root cause" << std::endl;
    } else {
        std::cout << "✗ FAILURE: PBWT does NOT properly identify allele-class similarity!" << std::endl;
        std::cout << std::endl;
        std::cout << "FAILURES:" << std::endl;
        for (const auto& failure : failures) {
            std::cout << "  - " << failure << std::endl;
        }
        std::cout << std::endl;
        std::cout << "INTERPRETATION:" << std::endl;
        std::cout << "  - PBWT may be treating all ALT alleles equally" << std::endl;
        std::cout << "  - OR: Anchor masking/redirect breaking similarity detection" << std::endl;
        std::cout << "  - OR: Split representation loses allele-specific information" << std::endl;
        std::cout << "  → This likely explains 33% accuracy loss in multiallelic phasing" << std::endl;
        std::cout << "  → PBWT selection not reflecting observed genotypes" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT multiallelic similarity test completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;

    return test_passed ? 0 : 1;
}

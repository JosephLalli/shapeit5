/*******************************************************************************
 * PBWT Conditioning Set Parity Test
 *
 * Tests whether PBWT conditioning set selection produces identical results
 * when processing the same genetic variation through:
 * 1. Biallelic pathway (variants processed separately)
 * 2. Supersite pathway (variants grouped at same position)
 *
 * This test will determine if the hypothesized "PBWT over-conditioning" 
 * actually occurs, or if conditioning sets remain identical between modes.
 *
 * Key Validation:
 * - Compare conditioning set sizes (K values) at each locus
 * - Compare actual donor haplotype IDs selected for conditioning
 * - Ensure identical haplotype similarity calculations
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

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

struct ConditioningComparison {
    int locus;
    int biallelic_k;
    int supersite_k;
    std::set<unsigned int> biallelic_donors;
    std::set<unsigned int> supersite_donors;
    double jaccard_similarity;
    bool sets_identical;
    bool significant_difference;
};

struct PBWTParityAnalysis {
    std::vector<ConditioningComparison> comparisons;
    int total_loci = 0;
    int loci_with_differences = 0;
    int loci_with_k_inflation = 0;
    double avg_k_ratio = 1.0;
    double avg_jaccard_similarity = 1.0;
    double min_jaccard_similarity = 1.0;
    
    void analyze() {
        if (comparisons.empty()) return;
        
        total_loci = comparisons.size();
        loci_with_differences = 0;
        loci_with_k_inflation = 0;
        double sum_k_ratio = 0.0;
        double sum_jaccard = 0.0;
        min_jaccard_similarity = 1.0;
        
        for (const auto& comp : comparisons) {
            if (!comp.sets_identical) {
                loci_with_differences++;
            }
            
            if (comp.supersite_k > comp.biallelic_k) {
                loci_with_k_inflation++;
            }
            
            double k_ratio = (comp.biallelic_k > 0) ? 
                            static_cast<double>(comp.supersite_k) / comp.biallelic_k : 1.0;
            sum_k_ratio += k_ratio;
            sum_jaccard += comp.jaccard_similarity;
            min_jaccard_similarity = std::min(min_jaccard_similarity, comp.jaccard_similarity);
        }
        
        avg_k_ratio = sum_k_ratio / total_loci;
        avg_jaccard_similarity = sum_jaccard / total_loci;
    }
    
    void print_summary() {
        std::cout << "\n=== PBWT CONDITIONING PARITY ANALYSIS ===" << std::endl;
        std::cout << "Total loci compared: " << total_loci << std::endl;
        std::cout << "Loci with different conditioning sets: " << loci_with_differences 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * loci_with_differences / total_loci) << "%)" << std::endl;
        std::cout << "Loci with K inflation: " << loci_with_k_inflation 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * loci_with_k_inflation / total_loci) << "%)" << std::endl;
        std::cout << "Average K ratio (SS/Bial): " << std::fixed << std::setprecision(3) << avg_k_ratio << std::endl;
        std::cout << "Average Jaccard similarity: " << std::fixed << std::setprecision(3) << avg_jaccard_similarity << std::endl;
        std::cout << "Minimum Jaccard similarity: " << std::fixed << std::setprecision(3) << min_jaccard_similarity << std::endl;
        
        if (loci_with_differences == 0) {
            std::cout << "\n✓ PBWT CONDITIONING SETS IDENTICAL" << std::endl;
            std::cout << "  No evidence of over-conditioning bug" << std::endl;
        } else {
            std::cout << "\n⚠️  PBWT CONDITIONING DIFFERENCES DETECTED" << std::endl;
            std::cout << "  Investigate why conditioning sets differ between modes" << std::endl;
        }
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

// Extract conditioning set for a specific locus from PBWT data structures
static std::set<unsigned int> extract_conditioning_set(const conditioning_set& H, int locus) {
    std::set<unsigned int> donors;
    
    if (locus >= static_cast<int>(H.sites_pbwt_selection.size()) || 
        !H.sites_pbwt_selection[locus]) {
        // Locus not selected for conditioning, return empty set
        return donors;
    }
    
    // Extract from indexes_pbwt_neighbour based on PBWT storage layout
    unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
    
    for (int target_hap = 0; target_hap < H.n_hap; ++target_hap) {
        int target_ind = target_hap / 2;
        if (target_ind >= H.n_ind) continue;
        
        unsigned long base_idx = H.sites_pbwt_grouping[locus] * 2UL * H.n_ind + target_hap;
        
        for (int d = 0; d < H.depth && base_idx < H.indexes_pbwt_neighbour.size(); ++d) {
            unsigned long neighbor_idx = d * addr_offset + base_idx;
            if (neighbor_idx < H.indexes_pbwt_neighbour.size()) {
                unsigned int donor_hap = H.indexes_pbwt_neighbour[neighbor_idx];
                donors.insert(donor_hap);
            }
        }
    }
    
    return donors;
}

static double compute_jaccard_similarity(const std::set<unsigned int>& set1, 
                                        const std::set<unsigned int>& set2) {
    if (set1.empty() && set2.empty()) return 1.0;
    
    std::set<unsigned int> intersection, union_set;
    
    std::set_intersection(set1.begin(), set1.end(),
                         set2.begin(), set2.end(),
                         std::inserter(intersection, intersection.begin()));
    
    std::set_union(set1.begin(), set1.end(),
                   set2.begin(), set2.end(),
                   std::inserter(union_set, union_set.begin()));
    
    return static_cast<double>(intersection.size()) / union_set.size();
}

static ConditioningComparison compare_conditioning_sets(
    int locus,
    const conditioning_set& H_bial,
    const conditioning_set& H_ss) {
    
    ConditioningComparison comp;
    comp.locus = locus;
    comp.biallelic_donors = extract_conditioning_set(H_bial, locus);
    comp.supersite_donors = extract_conditioning_set(H_ss, locus);
    comp.biallelic_k = static_cast<int>(comp.biallelic_donors.size());
    comp.supersite_k = static_cast<int>(comp.supersite_donors.size());
    comp.jaccard_similarity = compute_jaccard_similarity(comp.biallelic_donors, comp.supersite_donors);
    comp.sets_identical = (comp.biallelic_donors == comp.supersite_donors);
    comp.significant_difference = (std::abs(comp.supersite_k - comp.biallelic_k) > 1) || 
                                 (comp.jaccard_similarity < 0.95);
    
    return comp;
}

} // namespace

int main() {
    TEST_INIT("test_pbwt_conditioning_parity");
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT Conditioning Set Parity Test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create datasets representing same genetic variation processed differently
    std::cout << "Creating test datasets..." << std::endl;
    
    variant_map V_bial, V_ss;
    
    // Biallelic version: Process triallelic positions as separate biallelic variants
    V_bial.push(make_var("1", 1000, "pos1000_A_T", "A", "T", 0));  // First split at pos 1000
    V_bial.push(make_var("1", 1000, "pos1000_A_G", "A", "G", 1));  // Second split at pos 1000
    V_bial.push(make_var("1", 2000, "pos2000_C_A", "C", "A", 2));  // First split at pos 2000
    V_bial.push(make_var("1", 2000, "pos2000_C_T", "C", "T", 3));  // Second split at pos 2000
    V_bial.push(make_var("1", 3000, "pos3000_G_C", "G", "C", 4));  // Simple biallelic at pos 3000
    V_bial.push(make_var("1", 4000, "pos4000_T_A", "T", "A", 5));  // Simple biallelic at pos 4000

    // Supersite version: Same variants but will be grouped into supersites
    V_ss.push(make_var("1", 1000, "ss1_v0", "A", "T", 0));        // Supersite 1, variant 0
    V_ss.push(make_var("1", 1000, "ss1_v1", "A", "G", 1));        // Supersite 1, variant 1  
    V_ss.push(make_var("1", 2000, "ss2_v0", "C", "A", 2));        // Supersite 2, variant 0
    V_ss.push(make_var("1", 2000, "ss2_v1", "C", "T", 3));        // Supersite 2, variant 1
    V_ss.push(make_var("1", 3000, "v3000", "G", "C", 4));         // Regular variant
    V_ss.push(make_var("1", 4000, "v4000", "T", "A", 5));         // Regular variant

    // Create identical conditioning sets
    conditioning_set H_bial, H_ss;
    const int n_ref_samples = 16; // 16 samples = 32 haplotypes
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

    // Set realistic haplotype patterns that create meaningful PBWT similarity structure
    std::vector<std::vector<int>> realistic_patterns = {
        // Variant 0 (pos1000 A->T): 8 ALT haplotypes
        {1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 1 (pos1000 A->G): 6 ALT haplotypes 
        {0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 2 (pos2000 C->A): 10 ALT haplotypes
        {0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0},
        // Variant 3 (pos2000 C->T): 4 ALT haplotypes  
        {0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        // Variant 4 (pos3000 G->C): 12 ALT haplotypes
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0},
        // Variant 5 (pos4000 T->A): 14 ALT haplotypes
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1}
    };

    // Set identical patterns for both datasets
    for (int i = 0; i < 6; ++i) {
        set_panel(H_bial, i, realistic_patterns[i]);
        set_panel(H_ss, i, realistic_patterns[i]);
    }

    std::cout << "  Built biallelic dataset: " << V_bial.size() << " variants" << std::endl;
    std::cout << "  Built supersite dataset: " << V_ss.size() << " variants" << std::endl;
    
    // Build supersites for supersite dataset
    SuperSiteContext ctx_ss = build_supersites(V_ss, H_ss);
    std::cout << "  Detected " << ctx_ss.super_sites.size() << " supersites" << std::endl;

    // Initialize and build PBWT conditioning sets for both modes
    std::cout << "\nBuilding PBWT conditioning sets..." << std::endl;
    
    float modulo_selection = 1.0f;     // Select all sites for conditioning
    float modulo_multithreading = 1.0f;
    float mdr = 0.1f;
    int depth = 8;                     // Number of conditioning donors per target
    int mac = 1;
    int nthread = 1;

    H_bial.initialize(V_bial, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);
    H_ss.initialize(V_ss, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    H_bial.select();
    H_ss.select();

    std::cout << "  Biallelic PBWT: " << H_bial.n_site << " sites, depth=" << H_bial.depth << std::endl;
    std::cout << "  Supersite PBWT: " << H_ss.n_site << " sites, depth=" << H_ss.depth << std::endl;

    // Compare conditioning sets at each locus
    std::cout << "\nComparing conditioning sets..." << std::endl;
    
    PBWTParityAnalysis analysis;
    
    // Compare corresponding loci between biallelic and supersite modes
    int min_sites = std::min(H_bial.n_site, H_ss.n_site);
    for (int locus = 0; locus < min_sites; ++locus) {
        ConditioningComparison comp = compare_conditioning_sets(locus, H_bial, H_ss);
        analysis.comparisons.push_back(comp);
    }
    
    // Analyze and report results
    analysis.analyze();
    analysis.print_summary();
    // Test results
    bool test_passed = (analysis.loci_with_differences == 0);
    
    if (test_passed) {
        std::cout << std::endl;
        std::cout << "✓ SUCCESS: PBWT conditioning sets identical between modes" << std::endl;
        std::cout << "✓ No evidence of over-conditioning bug" << std::endl;
        std::cout << "✓ K inflation must have different root cause" << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << "✗ FAILURE: PBWT conditioning sets differ between modes!" << std::endl;
        std::cout << "✗ " << analysis.loci_with_differences << "/" << analysis.total_loci 
                  << " loci show conditioning differences" << std::endl;
        std::cout << "✗ This confirms PBWT over-conditioning hypothesis" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT conditioning parity test completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    return test_passed ? 0 : 1;
}

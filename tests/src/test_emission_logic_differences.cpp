/*******************************************************************************
 * Emission Logic Differences Test
 *
 * Investigates how supersite emission logic differs from biallelic emission
 * and whether this causes more donors to match (leading to K inflation).
 *
 * Key Hypothesis: The supersite anchor split semantics in match mask building
 * creates different matching patterns compared to biallelic logic, potentially
 * causing more donors to be considered as "matches", inflating K.
 *
 * Analysis:
 * - Biallelic: Direct bit comparison (donor_alt == expected_bit)  
 * - Supersite: Two-step anchor class mapping then split semantics
 *   - use_split=true uses anchor_class + amb_mask logic
 *   - use_split=false falls back to lane_class (similar to biallelic)
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

#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/models/site_emission_adapter.h"

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
    std::vector<uint8_t> sample_codes_unused;
};

struct EmissionComparison {
    int locus;
    int n_donors;
    int biallelic_matches;
    int supersite_matches;
    int supersite_split_matches;
    double match_ratio;
    bool significant_difference;
    std::vector<uint8_t> biallelic_mask;
    std::vector<uint8_t> supersite_mask;
    std::vector<uint8_t> supersite_split_mask;
};

struct EmissionAnalysis {
    std::vector<EmissionComparison> comparisons;
    int total_loci_compared = 0;
    int loci_with_differences = 0;
    int loci_with_inflation = 0;
    double avg_match_ratio = 1.0;
    double max_match_ratio = 1.0;
    int max_absolute_inflation = 0;
    
    void analyze() {
        if (comparisons.empty()) return;
        
        total_loci_compared = comparisons.size();
        loci_with_differences = 0;
        loci_with_inflation = 0;
        double sum_ratios = 0.0;
        max_match_ratio = 0.0;
        max_absolute_inflation = 0;
        
        for (const auto& comp : comparisons) {
            sum_ratios += comp.match_ratio;
            max_match_ratio = std::max(max_match_ratio, comp.match_ratio);
            
            if (comp.significant_difference) {
                loci_with_differences++;
                
                if (comp.supersite_matches > comp.biallelic_matches ||
                    comp.supersite_split_matches > comp.biallelic_matches) {
                    loci_with_inflation++;
                    
                    int inflation = std::max(comp.supersite_matches - comp.biallelic_matches,
                                           comp.supersite_split_matches - comp.biallelic_matches);
                    max_absolute_inflation = std::max(max_absolute_inflation, inflation);
                }
            }
        }
        
        avg_match_ratio = sum_ratios / total_loci_compared;
    }
    
    void print_summary() {
        std::cout << "\n=== EMISSION LOGIC ANALYSIS ===" << std::endl;
        std::cout << "Total loci compared: " << total_loci_compared << std::endl;
        std::cout << "Loci with differences: " << loci_with_differences 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * loci_with_differences / total_loci_compared) << "%)" << std::endl;
        std::cout << "Loci with match inflation: " << loci_with_inflation 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * loci_with_inflation / total_loci_compared) << "%)" << std::endl;
        std::cout << "Average match ratio: " << std::fixed << std::setprecision(3) << avg_match_ratio << std::endl;
        std::cout << "Maximum match ratio: " << std::fixed << std::setprecision(3) << max_match_ratio << std::endl;
        std::cout << "Maximum absolute inflation: " << max_absolute_inflation << " donors" << std::endl;
        
        if (loci_with_inflation > 0) {
            std::cout << "\n🔴 EMISSION LOGIC DIFFERENCES DETECTED!" << std::endl;
            std::cout << "Supersite emission creates more matches than biallelic at " 
                      << loci_with_inflation << " loci" << std::endl;
        } else {
            std::cout << "\n🟢 No significant emission differences detected" << std::endl;
        }
    }
    
    void print_details() {
        if (!env_true("SHAPEIT5_TEST_TRACE")) return;
        
        std::cout << "\n=== DETAILED EMISSION COMPARISON ===" << std::endl;
        for (const auto& comp : comparisons) {
            if (!comp.significant_difference) continue;
            
            std::cout << "\nLocus " << comp.locus << ":" << std::endl;
            std::cout << "  Donors: " << comp.n_donors << std::endl;
            std::cout << "  Biallelic matches: " << comp.biallelic_matches << std::endl;
            std::cout << "  Supersite matches: " << comp.supersite_matches << std::endl;
            std::cout << "  Supersite split matches: " << comp.supersite_split_matches << std::endl;
            std::cout << "  Match ratio: " << std::fixed << std::setprecision(3) << comp.match_ratio << std::endl;
        }
    }
};

enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static void clear_variant_state(genotype& G, int locus) {
    unsigned char& byte = G.Variants[DIV2(locus)];
    const int shift = (MOD2(locus)) << 2;
    byte &= ~(0x0F << shift);
}

static void set_phase(genotype& G, int locus, PhaseCode code) {
    clear_variant_state(G, locus);
    unsigned char& byte = G.Variants[DIV2(locus)];
    switch (code) {
        case REF_REF: break;
        case ALT_ALT:
            VAR_SET_HAP0(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
        case ALT_REF:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP0(MOD2(locus), byte);
            break;
        case REF_ALT:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
    }
}

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index, ctx.sample_codes_unused);
    return ctx;
}

static EmissionComparison compare_emission_logic(
    int locus,
    const genotype& G_bial, const conditioning_set& H_bial,
    const genotype& G_ss, const conditioning_set& H_ss, 
    const SuperSiteContext& ctx_ss,
    const std::vector<unsigned int>& idxH) {
    
    EmissionComparison comp;
    comp.locus = locus;
    comp.n_donors = static_cast<int>(idxH.size());
    comp.biallelic_matches = 0;
    comp.supersite_matches = 0; 
    comp.supersite_split_matches = 0;
    comp.match_ratio = 1.0;
    comp.significant_difference = false;
    
    // Build biallelic emission view and mask
    BiallelicEmissionAdapter bial_adapter(&G_bial, const_cast<bitmatrix*>(&H_bial.H_opt_var));
    SiteView bial_view;
    MatchMask bial_mask;
    bial_adapter.build_view(locus, 0, bial_view);
    bial_adapter.build_match_mask(bial_view, static_cast<unsigned int>(idxH.size()), locus, bial_mask);
    
    // Count biallelic matches
    for (size_t i = 0; i < bial_mask.by_donor_lane.size(); ++i) {
        if (bial_mask.by_donor_lane[i] == MatchMask::kMatch) {
            comp.biallelic_matches++;
        }
    }
    comp.biallelic_mask.assign(bial_mask.by_donor_lane.begin(), bial_mask.by_donor_lane.end());
    
    // Build supersite emission view and mask (standard mode)
    SupersiteEmissionAdapter ss_adapter(&G_ss, &ctx_ss.super_sites, &ctx_ss.locus_to_super_idx,
                                       &ctx_ss.super_site_var_index,
                                       ctx_ss.packed_codes.empty() ? nullptr : ctx_ss.packed_codes.data(),
                                       &idxH, ctx_ss.packed_codes.size());
    
    SiteView ss_view;
    MatchMask ss_mask;
    bool has_supersite = ss_adapter.build_view(locus, 0, ss_view);
    
    if (has_supersite) {
        // Standard supersite matching (use_anchor_split_semantics = false)
        ss_adapter.build_match_mask(ss_view, static_cast<unsigned int>(idxH.size()), false, ss_mask);
        
        for (size_t i = 0; i < ss_mask.by_donor_lane.size(); ++i) {
            if (ss_mask.by_donor_lane[i] == MatchMask::kMatch) {
                comp.supersite_matches++;
            }
        }
        comp.supersite_mask.assign(ss_mask.by_donor_lane.begin(), ss_mask.by_donor_lane.end());
        
        // Split semantics supersite matching (use_anchor_split_semantics = true)
        MatchMask ss_split_mask;
        ss_adapter.build_match_mask(ss_view, static_cast<unsigned int>(idxH.size()), true, ss_split_mask);
        
        for (size_t i = 0; i < ss_split_mask.by_donor_lane.size(); ++i) {
            if (ss_split_mask.by_donor_lane[i] == MatchMask::kMatch) {
                comp.supersite_split_matches++;
            }
        }
        comp.supersite_split_mask.assign(ss_split_mask.by_donor_lane.begin(), ss_split_mask.by_donor_lane.end());
        
        // Calculate match ratio and detect significant differences
        comp.match_ratio = (comp.biallelic_matches > 0) 
                          ? static_cast<double>(std::max(comp.supersite_matches, comp.supersite_split_matches)) / comp.biallelic_matches
                          : 1.0;
        
        comp.significant_difference = (std::abs(comp.supersite_matches - comp.biallelic_matches) > 0) ||
                                     (std::abs(comp.supersite_split_matches - comp.biallelic_matches) > 0);
        
        if (env_true("SHAPEIT5_TEST_TRACE") && comp.significant_difference) {
            std::cout << "EMISSION_DIFF: locus=" << locus 
                      << " bial=" << comp.biallelic_matches
                      << " ss=" << comp.supersite_matches 
                      << " ss_split=" << comp.supersite_split_matches
                      << " ratio=" << std::fixed << std::setprecision(3) << comp.match_ratio << std::endl;
        }
    }
    
    return comp;
}

} // namespace

int main() {
    std::cout << "======================================================================" << std::endl;
    std::cout << "Emission Logic Differences Investigation" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create identical datasets: biallelic vs supersite versions
    std::cout << "Creating test datasets..." << std::endl;
    
    variant_map V_bial, V_ss;
    
    // Biallelic version: 6 separate variants
    V_bial.push(make_var("1", 1000, "v1000", "A", "T", 0));
    V_bial.push(make_var("1", 2000, "v2000", "A", "C", 1));  
    V_bial.push(make_var("1", 3000, "v3000", "A", "G", 2));
    V_bial.push(make_var("1", 4000, "v4000", "A", "T", 3));
    V_bial.push(make_var("1", 5000, "v5000", "A", "C", 4));
    V_bial.push(make_var("1", 6000, "v6000", "A", "G", 5));

    // Supersite version: 6 variants forming 3 supersites
    V_ss.push(make_var("1", 1000, "ss1_anchor", "A", "T", 0));   
    V_ss.push(make_var("1", 1000, "ss1_sibling", "A", "C", 1));  
    V_ss.push(make_var("1", 3000, "ss2_anchor", "A", "G", 2));   
    V_ss.push(make_var("1", 3000, "ss2_sibling", "A", "T", 3));  
    V_ss.push(make_var("1", 5000, "ss3_anchor", "A", "C", 4));   
    V_ss.push(make_var("1", 5000, "ss3_sibling", "A", "G", 5));  

    // Create conditioning sets with diverse patterns to test emission differences
    conditioning_set H_bial, H_ss;
    H_bial.allocate(0, 12, V_bial.size()); // 12 reference samples => 24 haplotypes
    H_ss.allocate(0, 12, V_ss.size());

    auto set_panel = [](conditioning_set& H, int locus, const std::vector<int>& alt_flags) {
        for (size_t hap = 0; hap < alt_flags.size(); ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    // Set diverse conditioning panel patterns to maximize emission differences
    std::vector<std::vector<int>> patterns = {
        {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0}, // alternating
        {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1}, // reverse alternating
        {1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0}, // pairs
        {1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1}, // triples
        {1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0}, // sixes
        {0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,1,1}  // reverse sixes
    };

    for (int i = 0; i < 6; ++i) {
        set_panel(H_bial, i, patterns[i]);
        set_panel(H_ss, i, patterns[i]);
    }

    // Create target genotypes with HET sites to trigger complex emission logic
    auto create_genotype = [](size_t n_variants) {
        genotype G(0);
        G.n_segments = 1;
        G.n_variants = n_variants;
        G.n_ambiguous = 0;
        G.n_missing = 0;
        G.n_transitions = 0;
        G.n_stored_transitionProbs = 0;
        G.n_storage_events = 0;
        G.double_precision = false;
        G.haploid = false;
        G.Variants.assign((n_variants + 1) / 2, 0);
        G.Ambiguous.clear();
        G.Diplotypes.assign(1, 1ull);
        G.Lengths.assign(1, static_cast<unsigned short>(n_variants));
        return G;
    };

    genotype G_bial = create_genotype(V_bial.size());
    genotype G_ss = create_genotype(V_ss.size());

    // Set genotypes with HET patterns to trigger ambiguous emission logic
    std::vector<PhaseCode> gt_patterns = {ALT_REF, REF_ALT, ALT_REF, REF_ALT, ALT_REF, REF_ALT};
    
    for (size_t i = 0; i < gt_patterns.size(); ++i) {
        set_phase(G_bial, i, gt_patterns[i]);
        set_phase(G_ss, i, gt_patterns[i]);
    }
    
    G_bial.build();
    G_ss.build();

    std::cout << "  Built biallelic dataset: " << V_bial.size() << " variants" << std::endl;
    std::cout << "  Built supersite dataset: " << V_ss.size() << " variants" << std::endl;

    // Build supersites
    SuperSiteContext ctx_ss = build_supersites(V_ss, H_ss);
    std::cout << "  Detected " << ctx_ss.super_sites.size() << " supersites" << std::endl;

    // Compare emission logic at each locus
    std::cout << "\nAnalyzing emission logic differences..." << std::endl;
    
    const std::vector<unsigned int> idxH = {0u,1u,2u,3u,4u,5u,6u,7u,8u,9u,10u,11u,12u,13u,14u,15u,16u,17u,18u,19u,20u,21u,22u,23u};
    
    EmissionAnalysis analysis;
    
    // Compare anchor loci (even indices in supersite dataset)
    for (int anchor_idx = 0; anchor_idx < 3; ++anchor_idx) {
        int bial_locus = anchor_idx * 2; // 0, 2, 4
        int ss_locus = anchor_idx * 2;   // 0, 2, 4 (anchors in supersite dataset)
        
        if (bial_locus < static_cast<int>(V_bial.size()) && ss_locus < static_cast<int>(V_ss.size())) {
            EmissionComparison comp = compare_emission_logic(
                bial_locus, G_bial, H_bial, G_ss, H_ss, ctx_ss, idxH);
            analysis.comparisons.push_back(comp);
        }
    }
    
    // Analyze results
    analysis.analyze();
    analysis.print_summary();
    analysis.print_details();

    // Test results
    bool test_passed = (analysis.loci_with_inflation == 0);
    
    if (test_passed) {
        std::cout << std::endl;
        std::cout << "✓ SUCCESS: No emission logic differences causing match inflation" << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << "✗ FAILURE: Emission logic differences detected!" << std::endl;
        std::cout << "  Supersite emission creates more matches than biallelic" << std::endl;
        std::cout << "  This could be the root cause of K inflation and accuracy issues" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "Emission logic investigation completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    return test_passed ? 0 : 1;
}

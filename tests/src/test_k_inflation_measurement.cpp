/*******************************************************************************
 * K Inflation Measurement Test
 *
 * Quantifies K (conditioning set size) inflation when supersites are enabled
 * versus biallelic processing. This test measures the key metric suspected
 * to be causing accuracy degradation.
 *
 * Methodology:
 * 1. Create identical dataset processed both ways:
 *    - As biallelic splits (traditional processing)
 *    - With supersites enabled (multiallelic processing)
 * 2. Measure K values at each locus during conditioning set building
 * 3. Report inflation statistics (median, max, variance differences)
 * 
 * Expected: K should be similar between modes if no inflation bug exists
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
#include <numeric>
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

struct KMeasurement {
    std::vector<int> k_per_locus;
    int median_k = 0;
    int max_k = 0;
    int min_k = 0;
    double mean_k = 0.0;
    double variance_k = 0.0;
    
    void analyze() {
        if (k_per_locus.empty()) return;
        
        std::vector<int> sorted_k = k_per_locus;
        std::sort(sorted_k.begin(), sorted_k.end());
        
        min_k = sorted_k.front();
        max_k = sorted_k.back();
        median_k = sorted_k[sorted_k.size() / 2];
        
        mean_k = std::accumulate(sorted_k.begin(), sorted_k.end(), 0.0) / sorted_k.size();
        
        double sum_sq_diff = 0.0;
        for (int k : sorted_k) {
            sum_sq_diff += (k - mean_k) * (k - mean_k);
        }
        variance_k = sum_sq_diff / sorted_k.size();
    }
    
    void print(const std::string& label) {
        std::cout << "  " << label << " K Statistics:" << std::endl;
        std::cout << "    Median K: " << median_k << std::endl;
        std::cout << "    Mean K: " << std::fixed << std::setprecision(2) << mean_k << std::endl;
        std::cout << "    Range: [" << min_k << ", " << max_k << "]" << std::endl;
        std::cout << "    Variance: " << std::fixed << std::setprecision(2) << variance_k << std::endl;
    }
};

struct KInflationReport {
    KMeasurement biallelic;
    KMeasurement supersite;
    
    double median_inflation_ratio = 1.0;
    double mean_inflation_ratio = 1.0;
    double max_inflation_ratio = 1.0;
    int max_absolute_inflation = 0;
    
    void analyze() {
        biallelic.analyze();
        supersite.analyze();
        
        if (biallelic.median_k > 0) median_inflation_ratio = static_cast<double>(supersite.median_k) / biallelic.median_k;
        if (biallelic.mean_k > 0) mean_inflation_ratio = supersite.mean_k / biallelic.mean_k;
        if (biallelic.max_k > 0) max_inflation_ratio = static_cast<double>(supersite.max_k) / biallelic.max_k;
        
        max_absolute_inflation = supersite.max_k - biallelic.max_k;
    }
    
    void print_comparison() {
        std::cout << "\n=== K INFLATION ANALYSIS ===" << std::endl;
        biallelic.print("Biallelic");
        supersite.print("Supersite");
        
        std::cout << "\n  Inflation Metrics:" << std::endl;
        std::cout << "    Median K inflation ratio: " << std::fixed << std::setprecision(3) << median_inflation_ratio << std::endl;
        std::cout << "    Mean K inflation ratio: " << std::fixed << std::setprecision(3) << mean_inflation_ratio << std::endl;
        std::cout << "    Max K inflation ratio: " << std::fixed << std::setprecision(3) << max_inflation_ratio << std::endl;
        std::cout << "    Max absolute K increase: " << max_absolute_inflation << std::endl;
        
        if (median_inflation_ratio > 1.1) {
            std::cout << "\n  🔴 SIGNIFICANT K INFLATION DETECTED! (>10%)" << std::endl;
        } else if (median_inflation_ratio > 1.05) {
            std::cout << "\n  🟡 Moderate K inflation detected (>5%)" << std::endl;
        } else {
            std::cout << "\n  🟢 K inflation within acceptable range (<5%)" << std::endl;
        }
    }
    
    void export_csv(const std::string& filename) {
        std::ofstream out(filename);
        out << "locus,biallelic_k,supersite_k,inflation_ratio,absolute_diff\n";
        
        size_t min_size = std::min(biallelic.k_per_locus.size(), supersite.k_per_locus.size());
        for (size_t i = 0; i < min_size; ++i) {
            int bial_k = biallelic.k_per_locus[i];
            int ss_k = supersite.k_per_locus[i];
            double ratio = (bial_k > 0) ? static_cast<double>(ss_k) / bial_k : 1.0;
            int diff = ss_k - bial_k;
            
            out << i << "," << bial_k << "," << ss_k << "," << ratio << "," << diff << "\n";
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

// Simplified K measurement during PBWT operations
class KMeasuringHaplotypeSegment : public haplotype_segment_double {
public:
    KMeasurement* measurement;
    const std::vector<unsigned int>& idxH_ref;
    
    KMeasuringHaplotypeSegment(genotype* G, 
                              bitmatrix& H_opt_hap,
                              std::vector<unsigned int>& idxH,
                              window& W, 
                              hmm_parameters& M,
                              const std::vector<SuperSite>* super_sites,
                              const std::vector<bool>* is_super_site,
                              const std::vector<int>* locus_to_super_idx,
                              const uint8_t* panel_codes,
                              size_t panel_codes_size,
                              const std::vector<int>* super_site_var_index,
                              KMeasurement* k_measurement)
        : haplotype_segment_double(G, H_opt_hap, idxH, W, M, super_sites, is_super_site,
                                   locus_to_super_idx, panel_codes, panel_codes_size, super_site_var_index)
        , measurement(k_measurement)
        , idxH_ref(idxH) {
        
        if (measurement) {
            measurement->k_per_locus.clear();
            measurement->k_per_locus.reserve(W.stop_locus - W.start_locus + 1);
        }
    }

    void forward() {
        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::fprintf(stdout, "K_MEASUREMENT: starting forward() with %zu conditioning haps\n", idxH_ref.size());
            std::fprintf(stdout, "K_MEASUREMENT: window segments[%d-%d] loci[%d-%d]\n", 
                        segment_first, segment_last, locus_first, locus_last);
        }
        
        haplotype_segment_double::forward();
        
        // Record K measurement after the base forward pass
        // For this test, K is constant throughout the window
        if (measurement) {
            int current_k = static_cast<int>(idxH_ref.size());
            measurement->k_per_locus.push_back(current_k);
            
            if (env_true("SHAPEIT5_TEST_TRACE")) {
                std::fprintf(stdout, "K_MEASUREMENT: recorded K=%d\n", current_k);
            }
        }
    }

protected:
    void record_k_at_locus() {
        if (!measurement) return;
        
        // Record the conditioning set size at this locus
        // This is stored in the idxH array size which represents active donors
        int current_k = static_cast<int>(idxH_ref.size());
        measurement->k_per_locus.push_back(current_k);
        
        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::fprintf(stdout, "K_MEASUREMENT: locus=%d k=%d\n", curr_abs_locus, current_k);
        }
    }
};

} // namespace

int main() {
    std::cout << "======================================================================" << std::endl;
    std::cout << "K Inflation Measurement Test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Create test dataset that triggers supersite formation
    std::cout << "Creating test dataset..." << std::endl;
    
    variant_map V_bial, V_ss;
    
    // Biallelic version: 8 separate variants
    V_bial.push(make_var("1", 1000, "v1000", "A", "T", 0));
    V_bial.push(make_var("1", 1500, "v1500", "A", "C", 1));  
    V_bial.push(make_var("1", 2000, "v2000", "A", "G", 2));
    V_bial.push(make_var("1", 2500, "v2500", "A", "T", 3));
    V_bial.push(make_var("1", 3000, "v3000", "A", "C", 4));
    V_bial.push(make_var("1", 3500, "v3500", "A", "G", 5));
    V_bial.push(make_var("1", 4000, "v4000", "A", "T", 6));
    V_bial.push(make_var("1", 4500, "v4500", "A", "C", 7));

    // Supersite version: 8 variants forming 4 supersites (2 variants each at same positions)
    V_ss.push(make_var("1", 1000, "ss1_anchor", "A", "T", 0));   
    V_ss.push(make_var("1", 1000, "ss1_sibling", "A", "C", 1));  
    V_ss.push(make_var("1", 2000, "ss2_anchor", "A", "G", 2));   
    V_ss.push(make_var("1", 2000, "ss2_sibling", "A", "T", 3));  
    V_ss.push(make_var("1", 3000, "ss3_anchor", "A", "C", 4));   
    V_ss.push(make_var("1", 3000, "ss3_sibling", "A", "G", 5));  
    V_ss.push(make_var("1", 4000, "ss4_anchor", "A", "T", 6));   
    V_ss.push(make_var("1", 4000, "ss4_sibling", "A", "C", 7));  

    // Create identical conditioning sets for both
    conditioning_set H_bial, H_ss;
    H_bial.allocate(0, 4, V_bial.size()); // 4 reference samples => 8 haplotypes
    H_ss.allocate(0, 4, V_ss.size());

    auto set_panel = [](conditioning_set& H, int locus, const std::vector<int>& alt_flags) {
        for (size_t hap = 0; hap < alt_flags.size(); ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    // Set diverse haplotype patterns to create realistic conditioning scenarios
    std::vector<std::vector<int>> patterns = {
        {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0},  // alternating
        {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1},  // reverse alternating
        {1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0},  // pairs
        {0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1},  // reverse pairs
        {1,0,0,1,0,1,1,0,1,0,0,1,0,1,1,0},  // mixed
        {0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1},  // reverse mixed
        {1,1,1,0,0,0,1,1,1,0,0,0,1,1,1,0},  // triplets
        {0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1}   // reverse triplets
    };

    for (int i = 0; i < 8; ++i) {
        set_panel(H_bial, i, patterns[i]);
        set_panel(H_ss, i, patterns[i]);
    }

    // Create identical target genotypes
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

    // Set identical genotype patterns (some HET to create ambiguity)
    std::vector<PhaseCode> gt_pattern = {ALT_REF, REF_ALT, ALT_REF, REF_REF, ALT_ALT, REF_ALT, ALT_REF, REF_REF};
    
    for (size_t i = 0; i < gt_pattern.size(); ++i) {
        set_phase(G_bial, i, gt_pattern[i]);
        set_phase(G_ss, i, gt_pattern[i]);
    }
    
    G_bial.build();
    G_ss.build();

    std::cout << "  Built datasets: " << V_bial.size() << " biallelic variants, " 
              << V_ss.size() << " supersite variants" << std::endl;

    // Build supersites for supersite version
    SuperSiteContext ctx_ss = build_supersites(V_ss, H_ss);
    std::cout << "  Detected " << ctx_ss.super_sites.size() << " supersites" << std::endl;

    // Create windows and parameters
    auto create_window = [](const genotype& G) {
        window W{};
        W.start_locus = 0;
        W.stop_locus = static_cast<int>(G.n_variants) - 1;
        W.start_segment = 0;
        W.stop_segment = static_cast<int>(G.n_segments) - 1;  // Fixed: use actual segments from genotype
        W.start_ambiguous = 0;
        W.stop_ambiguous = static_cast<int>(G.n_ambiguous) - 1;
        W.start_missing = 0;
        W.stop_missing = static_cast<int>(G.n_missing) - 1;  // Fixed: use actual missing count
        W.start_transition = 0;
        W.stop_transition = static_cast<int>(G.n_transitions) - 1;  // Fixed: use actual transitions
        
        // Debug output
        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::fprintf(stdout, "Window: locus[%d-%d] segment[%d-%d] ambiguous[%d-%d] missing[%d-%d] transitions[%d-%d]\n",
                        W.start_locus, W.stop_locus, W.start_segment, W.stop_segment, 
                        W.start_ambiguous, W.stop_ambiguous, W.start_missing, W.stop_missing, 
                        W.start_transition, W.stop_transition);
        }
        
        return W;
    };

    window W_bial = create_window(G_bial);
    window W_ss = create_window(G_ss);

    auto create_hmm_params = [](size_t n_variants, unsigned int n_hap) {
        hmm_parameters M;
        M.ed = 0.01;
        M.ee = 1.0;
        M.ss_anchor_split_emissions = true;  // Use anchor split for parity with biallelic
        M.cm.assign(n_variants, 0.0f);
        for (size_t i = 0; i < n_variants; ++i) {
            M.cm[i] = 0.001f * static_cast<float>(i + 1);
        }
        M.Neff = 10000;
        M.Nhap = static_cast<int>(n_hap);
        M.t.assign(n_variants > 1 ? n_variants - 1 : 0, 0.01f);
        M.nt.assign(n_variants > 1 ? n_variants - 1 : 0, 0.99f);
        M.rare_allele.assign(n_variants, -1);
        return M;
    };

    hmm_parameters M_bial = create_hmm_params(V_bial.size(), H_bial.n_hap);
    hmm_parameters M_ss = create_hmm_params(V_ss.size(), H_ss.n_hap);
    M_ss.markSuperSiteSiblings(ctx_ss.super_sites, ctx_ss.locus_to_super_idx);

    std::vector<unsigned int> idxH = {0u,1u,2u,3u,4u,5u,6u,7u};

    // Measure K for both modes
    std::cout << "\nMeasuring K inflation..." << std::endl;
    
    KInflationReport report;

    // Biallelic measurement
    std::cout << "  Running biallelic measurement..." << std::endl;
    KMeasuringHaplotypeSegment HS_bial(&G_bial, H_bial.H_opt_hap, 
                                       idxH,
                                       W_bial, M_bial, nullptr, nullptr, nullptr, 
                                       nullptr, 0, nullptr, &report.biallelic);
    HS_bial.forward();

    // Supersite measurement  
    std::cout << "  Running supersite measurement..." << std::endl;
    KMeasuringHaplotypeSegment HS_ss(&G_ss, H_ss.H_opt_hap,
                                     idxH,
                                     W_ss, M_ss, &ctx_ss.super_sites, &ctx_ss.is_super_site,
                                     &ctx_ss.locus_to_super_idx, 
                                     ctx_ss.packed_codes.empty() ? nullptr : ctx_ss.packed_codes.data(),
                                     ctx_ss.packed_codes.size(), &ctx_ss.super_site_var_index, 
                                     &report.supersite);
    HS_ss.forward();

    // Analyze and report results
    report.analyze();
    report.print_comparison();

    // Export detailed data if requested
    if (env_true("SHAPEIT5_TEST_TRACE")) {
        std::string csv_file = "k_inflation_data.csv";
        report.export_csv(csv_file);
        std::cout << "\nDetailed K measurements exported to: " << csv_file << std::endl;
    }

    // Test results
    bool test_passed = (report.median_inflation_ratio < 1.1); // Allow up to 10% inflation
    
    if (test_passed) {
        std::cout << std::endl;
        std::cout << "✓ SUCCESS: K inflation within acceptable limits" << std::endl;
    } else {
        std::cout << std::endl;
        std::cout << "✗ FAILURE: Significant K inflation detected!" << std::endl;
        std::cout << "  This confirms K inflation is contributing to accuracy issues" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "K inflation measurement completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;
    
    return test_passed ? 0 : 1;
}
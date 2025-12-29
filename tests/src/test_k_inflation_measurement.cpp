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

#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h" 
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/compute_job.h"

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
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
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
        haplotype_segment_double::forward();
        // Record K measurement after the base forward pass
        // For this test, K is constant throughout the window
        if (measurement) {
            int current_k = static_cast<int>(idxH_ref.size());
            measurement->k_per_locus.push_back(current_k);
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
    TEST_INIT("test_k_inflation_measurement");
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
    // Increase reference panel: 16 reference samples => 32 haplotypes (more realistic)
    const int N_REF_SAMPLES = 16;
    // Include one main sample (the target genotype) so PBWT indexes have valid bounds
    H_bial.allocate(1, N_REF_SAMPLES, V_bial.size());
    H_ss.allocate(1, N_REF_SAMPLES, V_ss.size());

    auto set_panel = [](conditioning_set& H, int locus, const std::vector<int>& alt_flags) {
        for (size_t hap = 0; hap < alt_flags.size(); ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    // Programmatically generate diverse haplotype patterns sized to the panel (N_REF_SAMPLES * 2 haps)
    const int N_HAPS = H_bial.n_hap;
    std::vector<std::vector<int>> patterns;
    patterns.reserve(N_REF_SAMPLES);
    for (int s = 0; s < N_REF_SAMPLES; ++s) {
        std::vector<int> pat;
        pat.reserve(N_HAPS);
        // Create pattern with varied periodicities and bit rotations so donors differ
        int shift = 1 + (s % 5);
        for (int h = 0; h < N_HAPS; ++h) {
            // Use a mixture of bit shifts and modular residues to create diversity
            int bit = ((h >> (shift % 6)) & 1) ^ ((s & 1) ? 1 : 0);
            // Add some additional periodic toggles for variety
            if (((h + s) % (2 + (s % 4))) == 0) bit = 1;
            // Keep target sample haplotypes (0/1) reference-only; donors start at hap2
            if (h < 2) bit = 0;
            pat.push_back(bit);
        }
        patterns.push_back(pat);
    }

    for (size_t v = 0; v < V_bial.size(); ++v) {
        const int pat_idx = static_cast<int>(v % patterns.size());
        set_panel(H_bial, static_cast<int>(v), patterns[pat_idx]);
        set_panel(H_ss, static_cast<int>(v), patterns[pat_idx]);
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
        G.Lengths_bio = G.Lengths;
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

    // Ensure variant allele counts are non-zero so PBWT evaluation can run
    for (size_t i = 0; i < V_bial.size(); ++i) {
        V_bial.vec_pos[i]->cref = 8;
        V_bial.vec_pos[i]->calt = 8;
        V_bial.vec_pos[i]->cmis = 0;
        // set a small genetic distance so grouping is non-negative
        V_bial.vec_pos[i]->cm = 0.001 * static_cast<double>(i);
    }
    for (size_t i = 0; i < V_ss.size(); ++i) {
        V_ss.vec_pos[i]->cref = 8;
        V_ss.vec_pos[i]->calt = 8;
        V_ss.vec_pos[i]->cmis = 0;
        V_ss.vec_pos[i]->cm = 0.001 * static_cast<double>(i);
    }

    // Create windows and parameters
    auto create_window = [](const genotype& G) {
        window W{};
        W.start_locus = 0;
        W.stop_locus = static_cast<int>(G.n_variants) - 1;
        W.start_segment = 0;
        W.stop_segment = 0;
        W.start_ambiguous = 0;
        W.stop_ambiguous = static_cast<int>(G.n_ambiguous) - 1;
        W.start_missing = 0;
        W.stop_missing = -1;
        W.start_transition = 0;
        W.stop_transition = -1;
        return W;
    };

    window W_bial = create_window(G_bial);
    window W_ss = create_window(G_ss);

    auto create_hmm_params = [](size_t n_variants, unsigned int n_hap) {
        hmm_parameters M;
        M.ed = 0.01;
        M.ee = 1.0;
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

    // We'll run three full epochs (PBWT selection + compute_job.make + update/transposition)
    std::cout << "\nMeasuring K inflation across 3 epochs..." << std::endl;

    KInflationReport report;

    // Initialize PBWT parameters (simple settings appropriate for tiny test)
    const float modulo_selection = 1.0f;    // coarse grouping
    const float modulo_multithreading = 1.0f;
    // Use a very permissive MDR so synthetic sites are eligible for PBWT evaluation
    const float mdr = 1e6f;
    const int depth = 4;                     // PBWT depth (number of neighbours to propose)
    const int mac = 0;                       // accept all sites for evaluation
    const int nthread = 1;

    H_bial.initialize(V_bial, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);
    H_ss.initialize(V_ss, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    // Debug: print PBWT grouping/evaluation mapping to debug select() bounds
    std::cout << "\n[DEBUG] Biallelic PBWT groups: size=" << H_bial.sites_pbwt_grouping.size()
              << " back=" << (H_bial.sites_pbwt_grouping.empty() ? -1 : H_bial.sites_pbwt_grouping.back())
              << " eval_count=" << std::count(H_bial.sites_pbwt_evaluation.begin(), H_bial.sites_pbwt_evaluation.end(), true)
              << std::endl;
    std::cout << "[DEBUG] Supersite PBWT groups: size=" << H_ss.sites_pbwt_grouping.size()
              << " back=" << (H_ss.sites_pbwt_grouping.empty() ? -1 : H_ss.sites_pbwt_grouping.back())
              << " eval_count=" << std::count(H_ss.sites_pbwt_evaluation.begin(), H_ss.sites_pbwt_evaluation.end(), true)
              << std::endl;

    // Print full grouping vectors for deeper inspection
    std::cout << "[DEBUG] H_bial grouping: ";
    for (size_t i = 0; i < H_bial.sites_pbwt_grouping.size(); ++i) std::cout << H_bial.sites_pbwt_grouping[i] << (i+1<H_bial.sites_pbwt_grouping.size()?",":"\n");
    std::cout << "[DEBUG] H_ss grouping: ";
    for (size_t i = 0; i < H_ss.sites_pbwt_grouping.size(); ++i) std::cout << H_ss.sites_pbwt_grouping[i] << (i+1<H_ss.sites_pbwt_grouping.size()?",":"\n");

    // Prepare genotype_set wrappers expected by compute_job
    genotype_set GS_bial, GS_ss;
    GS_bial.allocate(1, V_bial.size());
    GS_ss.allocate(1, V_ss.size());

    // Copy constructed genotypes into genotype_set objects
    // (allocate created default genotype objects; overwrite their fields)
    GS_bial.vecG[0]->n_segments = G_bial.n_segments;
    GS_bial.vecG[0]->n_variants = G_bial.n_variants;
    GS_bial.vecG[0]->n_ambiguous = G_bial.n_ambiguous;
    GS_bial.vecG[0]->n_missing = G_bial.n_missing;
    GS_bial.vecG[0]->n_transitions = G_bial.n_transitions;
    GS_bial.vecG[0]->n_stored_transitionProbs = G_bial.n_stored_transitionProbs;
    GS_bial.vecG[0]->n_storage_events = G_bial.n_storage_events;
    GS_bial.vecG[0]->double_precision = G_bial.double_precision;
    GS_bial.vecG[0]->haploid = G_bial.haploid;
    GS_bial.vecG[0]->Variants = G_bial.Variants;
    GS_bial.vecG[0]->Ambiguous = G_bial.Ambiguous;
    GS_bial.vecG[0]->Diplotypes = G_bial.Diplotypes;
    GS_bial.vecG[0]->Lengths = G_bial.Lengths;

    GS_ss.vecG[0]->n_segments = G_ss.n_segments;
    GS_ss.vecG[0]->n_variants = G_ss.n_variants;
    GS_ss.vecG[0]->n_ambiguous = G_ss.n_ambiguous;
    GS_ss.vecG[0]->n_missing = G_ss.n_missing;
    GS_ss.vecG[0]->n_transitions = G_ss.n_transitions;
    GS_ss.vecG[0]->n_stored_transitionProbs = G_ss.n_stored_transitionProbs;
    GS_ss.vecG[0]->n_storage_events = G_ss.n_storage_events;
    GS_ss.vecG[0]->double_precision = G_ss.double_precision;
    GS_ss.vecG[0]->haploid = G_ss.haploid;
    GS_ss.vecG[0]->Variants = G_ss.Variants;
    GS_ss.vecG[0]->Ambiguous = G_ss.Ambiguous;
    GS_ss.vecG[0]->Diplotypes = G_ss.Diplotypes;
    GS_ss.vecG[0]->Lengths = G_ss.Lengths;

    // Three epochs: PBWT selection -> compute_job.make() -> update haplotypes -> transpose
    const int N_EPOCHS = 3;
    for (int epoch = 0; epoch < N_EPOCHS; ++epoch) {
        std::cout << "\nEpoch " << epoch+1 << " / " << N_EPOCHS << std::endl;

        // --- Biallelic path ---
        std::cout << "  [Biallelic] PBWT selection..." << std::endl;
        H_bial.select();

        // Run compute_job to build conditioning K-states for sample 0
        compute_job CJ_bial(V_bial, GS_bial, H_bial, 1024, 1024, nullptr, nullptr, nullptr);
        CJ_bial.make(0, 0.0);

        // Record K (number of conditioning haplotypes) for each window
        for (size_t w = 0; w < CJ_bial.Kstates.size(); ++w) {
            int ksize = static_cast<int>(CJ_bial.Kstates[w].size());
            report.biallelic.k_per_locus.push_back(ksize);
            if (env_true("SHAPEIT5_TEST_TRACE")) std::cout << "    K[bial][w=" << w << "]=" << ksize << std::endl;
        }

        // Mimic panel update & transpose (no sampling in this test)
        H_bial.updateHaplotypes(GS_bial, false);
        H_bial.transposeHaplotypes_H2V(true, false);

        // --- Supersite path ---
        std::cout << "  [Supersite] PBWT selection..." << std::endl;
        H_ss.select();

        compute_job CJ_ss(V_ss, GS_ss, H_ss, 1024, 1024, &ctx_ss.super_sites, &ctx_ss.locus_to_super_idx, &ctx_ss.super_site_var_index);
        CJ_ss.make(0, 0.0);

        for (size_t w = 0; w < CJ_ss.Kstates.size(); ++w) {
            int ksize = static_cast<int>(CJ_ss.Kstates[w].size());
            report.supersite.k_per_locus.push_back(ksize);
            if (env_true("SHAPEIT5_TEST_TRACE")) std::cout << "    K[ss][w=" << w << "]=" << ksize << std::endl;
        }

        H_ss.updateHaplotypes(GS_ss, false);
        H_ss.transposeHaplotypes_H2V(true, false);
    }

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

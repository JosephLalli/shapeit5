/*******************************************************************************
 * PBWT Biallelic Similarity Sanity Test
 *
 * Builds a small but polymorphic biallelic panel with simple LD structure,
 * runs the standard PBWT selection (no supersites), and checks that PBWT
 * neighbors are enriched for sharing the target allele at selected loci
 * compared to a global baseline.
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
#include <random>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

#define private public
#define protected public
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#undef private
#undef protected

#include "../../phase_common/src/containers/variant_map.h"

namespace {

struct SimilarityStats {
    double sum_matching_frac = 0.0;
    double sum_global_frac = 0.0;
    int n_pairs = 0;

    int n_pairs_better = 0;
    int n_pairs_much_better = 0;
};

static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

// Extract donors for a given target haplotype at a PBWT-selected locus
static std::vector<unsigned int> extract_donors_for_haplotype(
    const conditioning_set& H,
    int locus,
    int target_hap) {

    std::vector<unsigned int> donors;

    if (locus >= static_cast<int>(H.sites_pbwt_selection.size()) ||
        !H.sites_pbwt_selection[locus]) {
        return donors;
    }

    unsigned long addr_offset = H.sites_pbwt_ngroups * H.n_ind * 2UL;
    unsigned long base_idx = static_cast<unsigned long>(H.sites_pbwt_grouping[locus]) * 2UL * H.n_ind
                           + static_cast<unsigned long>(target_hap);

    for (int d = 0; d < H.depth; ++d) {
        unsigned long neighbor_idx = static_cast<unsigned long>(d) * addr_offset + base_idx;
        if (neighbor_idx >= H.indexes_pbwt_neighbour.size()) continue;
        int donor_raw = H.indexes_pbwt_neighbour[neighbor_idx];
        if (donor_raw >= 0) {
            donors.push_back(static_cast<unsigned int>(donor_raw));
        }
    }

    return donors;
}

} // namespace

int main() {
    TEST_INIT("test_pbwt_multiallelic_similarity_bial");
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT Biallelic Similarity Sanity Test" << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << std::endl;

    // Synthetic biallelic panel with simple LD structure
    const int n_sites = 100;
    const int n_test_samples = 10;
    const int n_ref_samples = 90;
    const int n_clusters = 4;

    std::cout << "Creating biallelic test dataset..." << std::endl;
    std::cout << "  Total sites: " << n_sites << std::endl;
    std::cout << "  Test samples: " << n_test_samples << " (" << n_test_samples * 2 << " haplotypes)" << std::endl;
    std::cout << "  Reference samples: " << n_ref_samples << " (" << (n_test_samples + n_ref_samples) * 2 << " haplotypes)" << std::endl;

    // Build variant map with simple genetic map
    variant_map V;
    for (int i = 0; i < n_sites; ++i) {
        int bp = 1000 + i * 100;
        std::string id = "v" + std::to_string(i);
        std::string ref = "A";
        std::string alt = "C";
        V.push(make_var("1", bp, id, ref, alt, i));
    }
    for (int i = 0; i < n_sites; ++i) {
        V.vec_pos[i]->cm = 0.01 * static_cast<double>(i + 1);
        V.vec_pos[i]->cref = 0;
        V.vec_pos[i]->calt = 0;
        V.vec_pos[i]->cmis = 0;
    }

    // Allocate conditioning set and haplotype matrices
    conditioning_set H;
    H.allocate(n_test_samples, n_ref_samples, V.size());

    // Build cluster-based haplotypes with simple LD structure
    std::mt19937 rng(12345);

    // Design: one dominant LD block where cluster membership strongly predicts allele,
    // remaining sites monomorphic REF. This makes the block (and especially its sites)
    // dominate PBWT similarity.
    const int block_start = 40;
    const int block_end = 60; // inclusive

    std::vector<std::vector<int>> cluster_alleles(n_clusters, std::vector<int>(n_sites, 0));
    for (int c = 0; c < n_clusters; ++c) {
        for (int s = 0; s < n_sites; ++s) {
            if (s >= block_start && s <= block_end) {
                // Clusters 0/1 carry ALT=1 in the block; 2/3 carry REF=0
                cluster_alleles[c][s] = (c < 2) ? 1 : 0;
            } else {
                cluster_alleles[c][s] = 0;
            }
        }
    }

    auto set_hap = [&](int hap_idx, int locus, bool value) {
        if (value) {
            H.H_opt_hap.set(hap_idx, locus, 1);
            H.H_opt_var.set(locus, hap_idx, 1);
        }
    };

    int n_hap = H.n_hap;
    for (int h = 0; h < n_hap; ++h) {
        int cluster = h % n_clusters;
        for (int s = 0; s < n_sites; ++s) {
            int allele = cluster_alleles[cluster][s];
            if (allele) set_hap(h, s, true);
        }
    }

    // Update variant MAC/MAF fields from H_opt_hap
    for (int locus = 0; locus < n_sites; ++locus) {
        unsigned int alt_count = 0;
        for (int h = 0; h < n_hap; ++h) {
            if (H.H_opt_hap.get(h, locus)) alt_count++;
        }
        V.vec_pos[locus]->calt = alt_count;
        V.vec_pos[locus]->cref = static_cast<unsigned int>(n_hap) - alt_count;
        V.vec_pos[locus]->cmis = 0;
    }

    // PBWT initialization: standard parameters
    std::cout << "\nBuilding PBWT (biallelic, no supersites)..." << std::endl;

    float modulo_selection = 0.05f;
    float modulo_multithreading = 1.0f;
    float mdr = 0.0f;
    int depth = 16;
    int mac = 1;
    int nthread = 1;

    H.initialize(V, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);
    H.select();

    std::cout << "  PBWT depth: " << H.depth << std::endl;
    std::cout << "  PBWT sites selected: " << H.n_site << std::endl;

    // Collect similarity statistics across selected loci and test haplotypes
    SimilarityStats stats;

    for (int locus = 0; locus < H.n_site; ++locus) {
        if (!H.sites_pbwt_selection[locus]) continue;

        unsigned int mac_at_locus = V.vec_pos[locus]->getMAC();
        if (mac_at_locus == 0) continue;

        for (int hap = 0; hap < n_test_samples * 2; ++hap) {
            bool target_allele = H.H_opt_hap.get(hap, locus);

            std::vector<unsigned int> donors = extract_donors_for_haplotype(H, locus, hap);
            if (donors.empty()) continue;

            int match_count = 0;
            for (unsigned int donor_hap : donors) {
                bool donor_allele = H.H_opt_hap.get(donor_hap, locus);
                if (donor_allele == target_allele) match_count++;
            }

            double matching_frac = static_cast<double>(match_count) / donors.size();

            int global_match = 0;
            for (int h = 0; h < n_hap; ++h) {
                bool allele = H.H_opt_hap.get(h, locus);
                if (allele == target_allele) global_match++;
            }
            double global_frac = static_cast<double>(global_match) / n_hap;

            stats.sum_matching_frac += matching_frac;
            stats.sum_global_frac += global_frac;
            stats.n_pairs++;
            if (matching_frac > global_frac) stats.n_pairs_better++;
            if (matching_frac >= global_frac + 0.15) stats.n_pairs_much_better++;
        }
    }

    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Similarity Statistics (PBWT neighbors vs global)" << std::endl;
    std::cout << "======================================================================" << std::endl;

    if (stats.n_pairs == 0) {
        std::cout << "No (haplotype, locus) pairs with PBWT neighbors found; test inconclusive." << std::endl;
        return 1;
    }

    double avg_matching = stats.sum_matching_frac / stats.n_pairs;
    double avg_global = stats.sum_global_frac / stats.n_pairs;

    std::cout << "  Pairs evaluated: " << stats.n_pairs << std::endl;
    std::cout << "  Average neighbor matching fraction: "
              << std::fixed << std::setprecision(3) << avg_matching << std::endl;
    std::cout << "  Average global matching fraction:   "
              << std::fixed << std::setprecision(3) << avg_global << std::endl;
    std::cout << "  Pairs with neighbor > global:       "
              << stats.n_pairs_better << " / " << stats.n_pairs
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * stats.n_pairs_better / stats.n_pairs) << "%)" << std::endl;
    std::cout << "  Pairs with neighbor ≥ global+0.15:  "
              << stats.n_pairs_much_better << " / " << stats.n_pairs
              << " (" << std::fixed << std::setprecision(1)
              << (100.0 * stats.n_pairs_much_better / stats.n_pairs) << "%)" << std::endl;

    std::cout << "\n======================================================================" << std::endl;
    std::cout << "Test Verdict" << std::endl;
    std::cout << "======================================================================" << std::endl;

    bool pass = true;
    std::vector<std::string> failures;

    const double min_avg_delta = 0.05;
    const double min_frac_better = 0.55;

    if (avg_matching < avg_global + min_avg_delta) {
        pass = false;
        failures.push_back("Average neighbor matching fraction not ≥ global + 0.05");
    }
    double frac_better = static_cast<double>(stats.n_pairs_better) / stats.n_pairs;
    if (frac_better < min_frac_better) {
        pass = false;
        failures.push_back("Fewer than 55% of pairs have neighbor > global");
    }

    if (pass) {
        std::cout << "✓ SUCCESS: PBWT neighbors enriched for matching alleles (bial path)" << std::endl;
    } else {
        std::cout << "✗ FAILURE: PBWT neighbors not clearly enriched over baseline (bial path)" << std::endl;
        std::cout << std::endl;
        std::cout << "FAILURES:" << std::endl;
        for (const auto& f : failures) {
            std::cout << "  - " << f << std::endl;
        }
    }

    std::cout << std::endl;
    std::cout << "======================================================================" << std::endl;
    std::cout << "PBWT biallelic similarity sanity test completed!" << std::endl;
    std::cout << "======================================================================" << std::endl;

    return pass ? 0 : 1;
}

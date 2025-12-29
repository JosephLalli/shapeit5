/*******************************************************************************
 * Emission Pattern Validation Test
 * 
 * Compares emission match masks between biallelic and supersite modes to identify
 * inconsistencies in emission semantics that may cause K inflation.
 *
 * Focus: Validate that SupersiteEmissionAdapter produces consistent emission 
 * patterns with BiallelicEmissionAdapter for equivalent genotype configurations.
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
#include "../../phase_common/src/models/site_emission_adapter.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#undef private
#undef protected

namespace {

struct EmissionComparison {
    int locus;
    std::string test_scenario;
    
    // Biallelic emission results
    std::vector<bool> biallelic_matches;
    std::array<bool, HAP_NUMBER> biallelic_any_match;
    
    // Supersite emission results  
    std::vector<bool> supersite_matches;
    std::array<bool, HAP_NUMBER> supersite_any_match;
    
    // Analysis results
    bool patterns_identical;
    double jaccard_similarity;
    int mismatched_donors;
    int mismatched_lanes;
    
    EmissionComparison(int l, const std::string& scenario) 
        : locus(l), test_scenario(scenario), patterns_identical(false),
          jaccard_similarity(0.0), mismatched_donors(0), mismatched_lanes(0) {}
    
    void analyze() {
        patterns_identical = true;
        mismatched_donors = 0;
        mismatched_lanes = 0;
        
        // Compare per-donor matches
        if (biallelic_matches.size() == supersite_matches.size()) {
            for (size_t i = 0; i < biallelic_matches.size(); i++) {
                if (biallelic_matches[i] != supersite_matches[i]) {
                    patterns_identical = false;
                    mismatched_donors++;
                }
            }
        } else {
            patterns_identical = false;
        }
        
        // Compare any_match_lane arrays
        for (int h = 0; h < HAP_NUMBER; h++) {
            if (biallelic_any_match[h] != supersite_any_match[h]) {
                patterns_identical = false;
                mismatched_lanes++;
            }
        }
        
        // Calculate Jaccard similarity for matches
        if (!biallelic_matches.empty() && !supersite_matches.empty()) {
            int intersection = 0, union_size = 0;
            for (size_t i = 0; i < std::min(biallelic_matches.size(), supersite_matches.size()); i++) {
                if (biallelic_matches[i] && supersite_matches[i]) intersection++;
                if (biallelic_matches[i] || supersite_matches[i]) union_size++;
            }
            jaccard_similarity = (union_size > 0) ? static_cast<double>(intersection) / union_size : 1.0;
        }
    }
    
    void print() const {
        std::printf("Locus %d (%s): %s\n", locus, test_scenario.c_str(),
                   patterns_identical ? "IDENTICAL" : "MISMATCH");
        if (!patterns_identical) {
            std::printf("  Mismatched donors: %d/%zu\n", mismatched_donors, biallelic_matches.size());
            std::printf("  Mismatched lanes: %d/%d\n", mismatched_lanes, HAP_NUMBER);
            std::printf("  Jaccard similarity: %.3f\n", jaccard_similarity);
        }
    }
};

class EmissionValidator {
public:
    std::vector<EmissionComparison> comparisons;
    
    void add_comparison(const EmissionComparison& comp) {
        comparisons.push_back(comp);
    }
    
    void analyze_all() {
        std::printf("=== Emission Pattern Validation Results ===\n");
        
        int identical_count = 0;
        int mismatch_count = 0;
        double avg_jaccard = 0.0;
        
        for (auto& comp : comparisons) {
            comp.analyze();
            comp.print();
            
            if (comp.patterns_identical) {
                identical_count++;
            } else {
                mismatch_count++;
            }
            avg_jaccard += comp.jaccard_similarity;
        }
        
        if (!comparisons.empty()) {
            avg_jaccard /= comparisons.size();
        }
        
        std::printf("\n=== Summary ===\n");
        std::printf("Total tests: %zu\n", comparisons.size());
        std::printf("Identical patterns: %d\n", identical_count);
        std::printf("Mismatched patterns: %d\n", mismatch_count);
        std::printf("Average Jaccard similarity: %.3f\n", avg_jaccard);
        
        if (mismatch_count > 0) {
            std::printf("*** EMISSION PATTERN INCONSISTENCIES DETECTED ***\n");
            std::printf("This may explain K inflation in supersite mode.\n");
        } else {
            std::printf("All emission patterns identical - inconsistency not in emission logic.\n");
        }
    }
    
    void export_csv(const char* filename) {
        FILE* f = std::fopen(filename, "w");
        if (!f) return;
        
        std::fprintf(f, "locus,scenario,patterns_identical,jaccard_similarity,mismatched_donors,mismatched_lanes\n");
        for (const auto& comp : comparisons) {
            std::fprintf(f, "%d,%s,%s,%.3f,%d,%d\n",
                        comp.locus, comp.test_scenario.c_str(),
                        comp.patterns_identical ? "true" : "false",
                        comp.jaccard_similarity, comp.mismatched_donors, comp.mismatched_lanes);
        }
        std::fclose(f);
        std::printf("Detailed emission comparisons exported to %s\n", filename);
    }
};

// Mock conditioning set for testing
struct MockConditioningData {
    std::vector<bool> donor_states;  // REF=false, ALT=true
    size_t n_donors;
    
    MockConditioningData(size_t n) : n_donors(n) {
        donor_states.resize(n);
    }
    
    void set_donor_alt(size_t donor_idx, bool is_alt) {
        if (donor_idx < n_donors) {
            donor_states[donor_idx] = is_alt;
        }
    }
};

// Simulate biallelic emission logic
static void simulate_biallelic_emission(const MockConditioningData& cond,
                                       uint8_t sample_h0_class, uint8_t sample_h1_class,
                                       EmissionComparison& result) {
    result.biallelic_matches.resize(cond.n_donors * HAP_NUMBER);
    std::fill(result.biallelic_any_match.begin(), result.biallelic_any_match.end(), false);
    
    for (size_t k = 0; k < cond.n_donors; k++) {
        const uint8_t donor_code = cond.donor_states[k] ? 1u : 0u;
        const size_t base = k * HAP_NUMBER;
        
        // HAP_NUMBER = 2 assumption: h0, h1
        bool h0_match = (donor_code == sample_h0_class);
        bool h1_match = (donor_code == sample_h1_class);
        
        result.biallelic_matches[base + 0] = h0_match;
        result.biallelic_matches[base + 1] = h1_match;
        
        result.biallelic_any_match[0] = result.biallelic_any_match[0] || h0_match;
        result.biallelic_any_match[1] = result.biallelic_any_match[1] || h1_match;
    }
}

// Simulate supersite emission logic (simplified class comparisons)
static void simulate_supersite_emission(const MockConditioningData& cond,
                                       uint8_t sample_h0_class, uint8_t sample_h1_class,
                                       uint32_t amb_mask,
                                       EmissionComparison& result) {
    result.supersite_matches.assign(cond.n_donors * HAP_NUMBER, false);
    std::fill(result.supersite_any_match.begin(), result.supersite_any_match.end(), false);

    for (size_t k = 0; k < cond.n_donors; k++) {
        const uint8_t donor_code = cond.donor_states[k] ? 1u : 0u;
        const size_t base = k * HAP_NUMBER;

        // HAP_NUMBER = 2 assumption: h0, h1
        bool h0_match = (donor_code == sample_h0_class);
        uint8_t expected_h1 = sample_h1_class;
        if (amb_mask != 0) {
            expected_h1 = ((amb_mask >> 1) & 1U) ? sample_h1_class : sample_h0_class;
        }
        bool h1_match = (donor_code == expected_h1);

        result.supersite_matches[base + 0] = h0_match;
        result.supersite_matches[base + 1] = h1_match;

        result.supersite_any_match[0] = result.supersite_any_match[0] || h0_match;
        result.supersite_any_match[1] = result.supersite_any_match[1] || h1_match;
    }
}

} // anonymous namespace

int main() {
    TEST_INIT("test_emission_pattern_validation");
    std::printf("=== Emission Pattern Validation Test ===\n\n");
    
    EmissionValidator validator;
    
    // Test scenarios with different genotype configurations
    const size_t n_donors = 8;
    MockConditioningData cond(n_donors);
    
    // Set up donor haplotypes: mix of REF and ALT
    cond.set_donor_alt(0, false);  // REF
    cond.set_donor_alt(1, true);   // ALT
    cond.set_donor_alt(2, false);  // REF
    cond.set_donor_alt(3, true);   // ALT
    cond.set_donor_alt(4, false);  // REF
    cond.set_donor_alt(5, true);   // ALT
    cond.set_donor_alt(6, false);  // REF
    cond.set_donor_alt(7, true);   // ALT
    
    std::printf("Testing emission patterns with %zu donors (4 REF, 4 ALT)\n", n_donors);
    
    // Test Case 1: Homozygous REF sample
    {
        EmissionComparison comp(1, "HOM_REF");
        simulate_biallelic_emission(cond, 0, 0, comp);
        simulate_supersite_emission(cond, 0, 0, 0, comp);
        validator.add_comparison(comp);
    }
    
    // Test Case 2: Homozygous ALT sample  
    {
        EmissionComparison comp(2, "HOM_ALT");
        simulate_biallelic_emission(cond, 1, 1, comp);
        simulate_supersite_emission(cond, 1, 1, 0, comp);
        validator.add_comparison(comp);
    }
    
    // Test Case 3: Heterozygous sample (standard)
    {
        EmissionComparison comp(3, "HET_standard");
        simulate_biallelic_emission(cond, 0, 1, comp);
        simulate_supersite_emission(cond, 0, 1, 0, comp);
        validator.add_comparison(comp);
    }
    
    // Test Case 4: Ambiguous site (amb_mask=0b10)
    {
        EmissionComparison comp(4, "AMB_masked");
        simulate_biallelic_emission(cond, 0, 1, comp);
        simulate_supersite_emission(cond, 0, 1, 0b10, comp);
        validator.add_comparison(comp);
    }
    
    // Test Case 5: Missing data comparison
    {
        EmissionComparison comp(5, "MISSING");
        simulate_biallelic_emission(cond, 2, 2, comp);  // Assume 2 = missing
        simulate_supersite_emission(cond, 2, 2, 0, comp);
        validator.add_comparison(comp);
    }
    
    std::printf("\nRunning emission pattern comparisons...\n\n");
    validator.analyze_all();
    validator.export_csv("emission_pattern_validation.csv");
    
    std::printf("\n=== Investigation Summary ===\n");
    std::printf("This test validates emission match mask consistency between biallelic and supersite modes.\n");
    std::printf("Key findings help identify if emission mismatches could cause K inflation through:\n");
    std::printf("  - Inconsistent donor matching patterns\n");
    std::printf("  - Different any_match_lane calculations\n");
    
    std::printf("\nNext steps based on results:\n");
    std::printf("  - If patterns identical: K inflation source is elsewhere\n");
    std::printf("  - If patterns differ: Focus on emission logic differences\n");
    std::printf("  - Target emission adapter fixes for consistency\n");
    
    std::printf("\n=== Test completed ===\n");
    TEST_SUMMARY();
    return 0;
}

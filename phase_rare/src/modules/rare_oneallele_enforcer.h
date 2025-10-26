/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#ifndef _RARE_ONEALLELE_ENFORCER_H
#define _RARE_ONEALLELE_ENFORCER_H

#include <utils/otools.h>
#include <objects/sparse_genotype.h>
#include <containers/variant_map.h>
#include <containers/genotype_set/genotype_set_header.h>
#include <phaser/phaser_header.h>

namespace shapeit5 {
namespace modules {

// Unify stats type with phaser: alias the nested type
using OneAlleleRareStats = ::phaser::OneAlleleRareStats;

/**
 * @brief Enhanced rare variant multiallelic constraint enforcement
 * 
 * Provides sophisticated algorithms for resolving multiallelic constraint violations
 * in sparse rare variant data, with support for PBWT donor integration and 
 * Li-Stephens emission modeling.
 */
class RareOneAlleleEnforcer {
public:
    
    // Enforcement mode options
    enum class Mode {
        PP_BASIC,           // Basic PP-based enforcement (existing algorithm)
        PP_ENHANCED,        // PP-based with Li-Stephens enhancement
        SPARSE_TRANSITION,  // Sparse transition scoring with donor context
        SPARSE_MICRO        // Sparse enumeration for complex cases
    };
    
    // Structure for tracking sparse multiallelic groups
    struct SparseMultiallelicGroup {
        std::vector<int> rare_variant_indices;                    // Indices into V.vec_rare
        std::vector<std::vector<sparse_genotype*>> sample_buckets; // Per-sample sparse entries
        long long position_key;                                   // (chr_id << 32) | bp
        int total_samples_with_violations = 0;
    };
    
    // Change analysis for tracking genotype vs phase modifications
    struct SparseChangeAnalysis {
        bool has_genotype_changes = false;     // REF<->ALT changes
        bool has_phase_only_changes = false;   // Only haplotype assignment changes
        std::vector<std::string> change_descriptions;  // For debugging/logging
        
        // Specific to sparse representation
        bool used_donor_context = false;       // PBWT donors were available and used
        bool used_li_stephens = false;         // Li-Stephens emission modeling applied
        int complexity_level = 0;              // 0=simple, 1=moderate, 2=complex
    };

public:
    RareOneAlleleEnforcer();
    ~RareOneAlleleEnforcer();
    
    // Main enforcement interface
    void enforce_multiallelic_constraints(
        genotype_set& G,
        variant_map& V,
        Mode enforcement_mode,
        OneAlleleRareStats& stats
    );
    
    // Mode-specific enforcement algorithms
    bool enforce_pp_basic(
        const SparseMultiallelicGroup& group,
        genotype_set& G,
        OneAlleleRareStats& stats
    );
    
    bool enforce_pp_enhanced(
        const SparseMultiallelicGroup& group,
        genotype_set& G,
        variant_map& V,
        OneAlleleRareStats& stats
    );
    
    bool enforce_sparse_transition(
        const SparseMultiallelicGroup& group,
        genotype_set& G,
        variant_map& V,
        OneAlleleRareStats& stats
    );
    
    bool enforce_sparse_micro(
        const SparseMultiallelicGroup& group,
        genotype_set& G,
        variant_map& V,
        OneAlleleRareStats& stats
    );

private:
    // Group building and analysis
    std::vector<SparseMultiallelicGroup> build_multiallelic_groups(
        const genotype_set& G,
        variant_map& V
    );
    
    bool analyze_violations_in_group(
        const SparseMultiallelicGroup& group,
        genotype_set& G
    );
    
    // Li-Stephens emission modeling for sparse data
    double compute_sparse_li_stephens_score(
        const std::vector<sparse_genotype*>& sample_genotypes,
        const std::vector<sparse_genotype*>& proposed_assignment
    );
    
    // PBWT donor context integration for sparse transition scoring
    double compute_sparse_donor_transition_score(
        const std::vector<sparse_genotype*>& sample_genotypes,
        const std::vector<sparse_genotype*>& proposed_assignment,
        int sample_id,
        const genotype_set& G,
        variant_map& V
    );
    
    // Sparse micro mode enumeration algorithms
    bool enumerate_optimal_sparse_assignment(
        const std::vector<sparse_genotype*>& het_genotypes,
        int sample_id,
        const genotype_set& G,
        variant_map& V,
        OneAlleleRareStats& stats
    );
    
    bool optimize_large_sparse_group(
        const std::vector<sparse_genotype*>& het_genotypes,
        int sample_id,
        const genotype_set& G,
        variant_map& V,
        OneAlleleRareStats& stats
    );
    
    // Change analysis and tracking
    SparseChangeAnalysis analyze_sparse_changes(
        const std::vector<sparse_genotype*>& before,
        const std::vector<sparse_genotype*>& after,
        Mode enforcement_mode
    );
    
    void update_statistics(
        const SparseChangeAnalysis& analysis,
        Mode enforcement_mode,
        OneAlleleRareStats& stats
    );
    
    // Utility functions
    double clamp_probability(double prob, double min_val = 1e-6, double max_val = 1.0 - 1e-6);
    double logit(double prob);
    double sigmoid(double logit_val);
    
    // Constants
    static constexpr double DEFAULT_ERROR_RATE = 0.01;
    static constexpr int MAX_ITERATIONS_PER_GROUP = 4;
    static constexpr int COMPLEX_GROUP_THRESHOLD = 3; // >3 ALTs = complex
};

}  // namespace modules
}  // namespace shapeit5

#endif /* _RARE_ONEALLELE_ENFORCER_H */

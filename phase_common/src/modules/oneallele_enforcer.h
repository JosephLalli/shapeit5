#pragma once

#include <cstdint>
#include <vector>
#include <string>

#include "multiallelic_position_map.h"

class genotype;
class genotype_set;
class variant_map;

namespace shapeit5 {
namespace modules {

enum class OneAlleleMode {
  TRANSITION,
  MICRO,
  MICRO_DONOR
};

struct OneAlleleStats {
  std::uint64_t positions_checked = 0;
  std::uint64_t violations_found = 0;
  std::uint64_t flips_applied = 0;
  
  // Micro-donor mode specific statistics
  std::uint64_t emission_dominated_decisions = 0;    // Emission score > transition score
  std::uint64_t transition_dominated_decisions = 0;  // Transition score > emission score
  std::uint64_t donor_weighted_changes = 0;          // Different outcome vs simple emission
  std::uint64_t genotype_changes = 0;                // REF<->ALT changes (not just phase)
  std::uint64_t phase_only_changes = 0;              // Only haplotype assignment changes
  std::uint64_t fallback_to_simple_emission = 0;     // No PBWT donors available
};

struct OneAlleleEpochStats {
  std::uint64_t violations_found = 0;
  std::uint64_t flips_applied = 0;
  
  // Micro-donor mode specific epoch statistics
  std::uint64_t emission_dominated_decisions = 0;
  std::uint64_t transition_dominated_decisions = 0;
  std::uint64_t donor_weighted_changes = 0;
  std::uint64_t genotype_changes = 0;
  std::uint64_t phase_only_changes = 0;
  std::uint64_t fallback_to_simple_emission = 0;
};

class OneAlleleEnforcer {
 public:
  OneAlleleEnforcer();
  void set_enabled(bool enabled);
  bool enabled() const;
  void set_mode(OneAlleleMode mode);
  OneAlleleMode mode() const;
  void reset_stats();
  const OneAlleleStats& stats() const;
  void reset_epoch_stats();
  const OneAlleleEpochStats& epoch_stats() const;

  void set_conditioning_size(int m);
  void set_min_distance_cm(double d);
  void set_debug_output_file(const std::string& filename);
  void set_debug_target(const std::string& sample_name, const std::string& chrom, int pos);

  void enforce(const MultiallelicPositionMap& map,
               genotype_set& G,
               const variant_map& V,
               const std::string& iteration_context = "unknown");

  void enforce_sample(const MultiallelicPositionMap& map,
                      genotype& g,
                      const variant_map& V,
                      const std::vector<std::vector<unsigned int>>& Kstates,  // Fix: unsigned int is correct
                      const std::string& iteration_context = "unknown",
                      int sample_index = -1);

  // Get per-sample statistics (for worker threads)
  const OneAlleleEpochStats& sample_epoch_stats() const;

  // Accumulate sample stats into global epoch stats (mutex-protected by caller)
  void accumulate_sample_stats(const OneAlleleEpochStats& sample_stats);

  // Reset per-sample statistics before enforce_sample() call
  void reset_sample_epoch_stats();

  // Check and log state of debug target from genotype_set 
  void check_debug_target_state(const std::string& stage,
                                 const MultiallelicPositionMap& map,
                                 genotype_set& G,
                                 const variant_map& V) const;

 private:
  bool enabled_ = false;
  OneAlleleMode mode_ = OneAlleleMode::TRANSITION;
  OneAlleleStats stats_;
  OneAlleleEpochStats epoch_stats_;
  OneAlleleEpochStats sample_epoch_stats_;
  int default_conditioning_size_ = 16;
  double min_distance_cm_ = 1e-8;
  std::string debug_output_file_;
  std::string current_iteration_context_;
  int current_sample_index_;
  std::string current_sample_name_;

  // Target debugging
  std::string debug_target_sample_;
  std::string debug_target_chrom_;
  int debug_target_pos_;
  bool debug_target_enabled_;

  void debug_log_multiallelic_site(const std::string& message) const;
  bool is_debug_target(const std::string& sample_name, const std::string& chrom, int pos) const;
  void log_debug_target_state(const std::string& stage,
                               const std::string& sample_name,
                               const std::string& chrom,
                               int pos,
                               const std::vector<int>& variant_indices,
                               const std::vector<uint8_t>& hap0_bits,
                               const std::vector<uint8_t>& hap1_bits,
                               const variant_map& V) const;

  int find_left_neighbor(const genotype& g,
                         const std::vector<int>& variant_to_segment,
                         int idx,
                         const variant_map& V) const;

  int find_right_neighbor(const genotype& g,
                          const std::vector<int>& variant_to_segment,
                          int idx,
                          const variant_map& V) const;

  void flip_alt_to_other_hap(genotype& g,
                             int variant_idx,
                             bool from_hap0,
                             std::vector<uint8_t>& hap0_bits,
                             std::vector<uint8_t>& hap1_bits) const;

  bool resolve_violation_for_hap(genotype& g,
                                 const variant_map& V,
                                 const std::vector<int>& variant_to_segment,
                                 std::vector<uint8_t>& target_hap,
                                 std::vector<uint8_t>& other_hap,
                                 std::vector<int>& alt_indices,
                                 bool target_is_hap0);

  // Micro re-decode functions
  struct MicroCandidate {
    std::vector<uint8_t> hap0_assignment;  // 0=REF, 1=ALT for each variant in group
    std::vector<uint8_t> hap1_assignment;
    double score;
    bool is_valid;
    
    // For tracking decision components (micro-donor mode)
    double transition_score = 0.0;
    double emission_score = 0.0;
  };
  
  struct GenotypeChangeAnalysis {
    bool has_genotype_changes = false;     // REF<->ALT changes
    bool has_phase_only_changes = false;   // Only haplotype assignment changes
    std::vector<std::string> change_descriptions;  // For verbose logging
  };

  // Enforcement mode implementations
  // Per-sample enforcement helpers (used by enforce_sample)
  bool enforce_group_transition(const PositionGroup& position_group,
                                genotype& g,
                                const variant_map& V,
                                const std::vector<int>& variant_to_segment,
                                std::vector<uint8_t>& hap0_bits,
                                std::vector<uint8_t>& hap1_bits,
                                const std::vector<int>& alt_indices_h0,
                                const std::vector<int>& alt_indices_h1);

  // Batch mode version (used by enforce()) - donor-agnostic
  bool enforce_group_micro(genotype& g,
                          const variant_map& V,
                          const std::vector<int>& variant_to_segment,
                          const std::vector<int>& position_group_indices,
                          std::vector<uint8_t>& hap0_bits,
                          std::vector<uint8_t>& hap1_bits);

  // Per-sample mode version (used by enforce_sample()) - with donor context
  bool enforce_group_micro(const PositionGroup& position_group,
                          genotype& g,
                          const variant_map& V,
                          const std::vector<std::vector<unsigned int>>& Kstates,
                          const std::vector<int>& variant_to_segment,
                          std::vector<uint8_t>& hap0_bits,
                          std::vector<uint8_t>& hap1_bits,
                          const std::vector<int>& alt_indices_h0,
                          const std::vector<int>& alt_indices_h1,
                          bool use_donors);
                          
// Candidate enumeration and scoring
  std::vector<MicroCandidate> enumerate_micro_candidates(
      const std::vector<int>& position_group_indices,
      const std::vector<uint8_t>& current_hap0,
      const std::vector<uint8_t>& current_hap1);

  double evaluate_candidate_micro(const MicroCandidate& candidate,
                                 genotype& g,
                                 const variant_map& V,
                                 const std::vector<int>& variant_to_segment,
                                 const std::vector<int>& position_group_indices);

  double evaluate_candidate_raw_score(const MicroCandidate& candidate,
                                     genotype& g,
                                     const variant_map& V,
                                     const std::vector<int>& variant_to_segment,
                                     const std::vector<int>& position_group_indices);

  double evaluate_candidate_micro_donor(const MicroCandidate& candidate,
                                       genotype& g,
                                       const variant_map& V,
                                       const std::vector<int>& variant_to_segment,
                                       const std::vector<int>& position_group_indices,
                                       const std::vector<std::vector<unsigned int>>& Kstates);

  double compute_transition_score(genotype& g,
                                 const variant_map& V,
                                 const std::vector<int>& variant_to_segment,
                                 int left_anchor,
                                 int right_anchor,
                                 const std::vector<int>& group_indices,
                                 const std::vector<uint8_t>& hap0_assignment,
                                 const std::vector<uint8_t>& hap1_assignment);

  double compute_emission_score(const std::vector<int>& group_indices,
                               const std::vector<uint8_t>& hap0_assignment,
                               const std::vector<uint8_t>& hap1_assignment);

  double compute_donor_weighted_emission_score(const std::vector<int>& group_indices,
                                              const std::vector<uint8_t>& hap0_assignment,
                                              const std::vector<uint8_t>& hap1_assignment,
                                              const variant_map& V,
                                              const std::vector<std::vector<unsigned int>>& Kstates);

  void apply_micro_candidate(genotype& g,
                            const MicroCandidate& candidate,
                            const std::vector<int>& position_group_indices,
                            std::vector<uint8_t>& hap0_bits,
                            std::vector<uint8_t>& hap1_bits);

  double compute_chain_score_with_donors(genotype& g,
                                         const variant_map& V,
                                         const std::vector<int>& variant_to_segment,
                                         int left_anchor,
                                         int right_anchor,
                                         const std::vector<int>& group_indices,
                                         const std::vector<uint8_t>& hap0_assignment,
                                         const std::vector<uint8_t>& hap1_assignment,
                                         const std::vector<std::vector<unsigned int>>& Kstates);

  double compute_flip_score(genotype& g,
                            const variant_map& V,
                            int variant_idx,
                            int left_anchor,
                            int right_anchor,
                            bool is_hap0) const;

  // Enhanced tracking and analysis methods
  GenotypeChangeAnalysis analyze_genotype_changes(const MicroCandidate& candidate,
                                                  const std::vector<int>& position_group_indices,
                                                  const std::vector<uint8_t>& current_hap0,
                                                  const std::vector<uint8_t>& current_hap1,
                                                  const variant_map& V) const;

  void log_genotype_changes(const GenotypeChangeAnalysis& analysis,
                           const std::string& sample_name,
                           const std::vector<int>& position_group_indices,
                           const variant_map& V) const;

  void update_micro_donor_stats(const MicroCandidate& best_candidate,
                               const GenotypeChangeAnalysis& analysis,
                               bool used_donor_weighting);
};

}  // namespace modules
}  // namespace shapeit5

#pragma once

#include <cstdint>
#include <vector>

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
};

struct OneAlleleEpochStats {
  std::uint64_t violations_found = 0;
  std::uint64_t flips_applied = 0;
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

  void enforce(const MultiallelicPositionMap& map,
               genotype_set& G,
               const variant_map& V);

  void enforce_sample(const MultiallelicPositionMap& map,
                      genotype& g,
                      const variant_map& V,
                      const std::vector<std::vector<unsigned int>>& Kstates);

 private:
  bool enabled_ = false;
  OneAlleleMode mode_ = OneAlleleMode::TRANSITION;
  OneAlleleStats stats_;
  OneAlleleEpochStats epoch_stats_;
  int default_conditioning_size_ = 16;
  double min_distance_cm_ = 1e-8;

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
  };

  bool enforce_group_micro(genotype& g,
                          const variant_map& V,
                          const std::vector<int>& variant_to_segment,
                          const std::vector<int>& position_group_indices,
                          std::vector<uint8_t>& hap0_bits,
                          std::vector<uint8_t>& hap1_bits);

  bool enforce_group_micro_donor(genotype& g,
                                 const variant_map& V,
                                 const std::vector<int>& variant_to_segment,
                                 const std::vector<int>& position_group_indices,
                                 std::vector<uint8_t>& hap0_bits,
                                 std::vector<uint8_t>& hap1_bits,
                                 const std::vector<std::vector<unsigned int>>& Kstates);

  std::vector<MicroCandidate> enumerate_micro_candidates(
      const std::vector<int>& position_group_indices,
      const std::vector<uint8_t>& current_hap0,
      const std::vector<uint8_t>& current_hap1);

  double evaluate_candidate_micro(const MicroCandidate& candidate,
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
};

}  // namespace modules
}  // namespace shapeit5

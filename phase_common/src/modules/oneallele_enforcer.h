#pragma once

#include <cstdint>
#include <vector>

#include "multiallelic_position_map.h"

class genotype;
class genotype_set;
class variant_map;

namespace shapeit5 {
namespace modules {

struct OneAlleleStats {
  std::uint64_t positions_checked = 0;
  std::uint64_t violations_found = 0;
  std::uint64_t flips_applied = 0;
};

class OneAlleleEnforcer {
 public:
  OneAlleleEnforcer();
  void set_enabled(bool enabled);
  bool enabled() const;
  void reset_stats();
  const OneAlleleStats& stats() const;

  void set_conditioning_size(int m);
  void set_min_distance_cm(double d);

  void enforce(const MultiallelicPositionMap& map,
               genotype_set& G,
               const variant_map& V);

 private:
  bool enabled_ = false;
  OneAlleleStats stats_;
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
};

}  // namespace modules
}  // namespace shapeit5

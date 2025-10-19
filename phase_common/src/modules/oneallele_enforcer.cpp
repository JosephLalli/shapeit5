#include "oneallele_enforcer.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include <containers/genotype_set.h>
#include <containers/variant_map.h>
#include <objects/genotype/genotype_header.h>

#include "transition_scorer.h"

namespace shapeit5 {
namespace modules {

OneAlleleEnforcer::OneAlleleEnforcer() = default;

void OneAlleleEnforcer::set_enabled(bool enabled) { enabled_ = enabled; }

bool OneAlleleEnforcer::enabled() const { return enabled_; }

void OneAlleleEnforcer::reset_stats() { stats_ = OneAlleleStats{}; }

const OneAlleleStats& OneAlleleEnforcer::stats() const { return stats_; }

void OneAlleleEnforcer::set_conditioning_size(int m) {
  default_conditioning_size_ = std::max(m, 2);
}

void OneAlleleEnforcer::set_min_distance_cm(double d) {
  if (d > 0.0) {
    min_distance_cm_ = d;
  }
}

void OneAlleleEnforcer::enforce(const MultiallelicPositionMap& map,
                                genotype_set& G,
                                const variant_map& V) {
  if (!enabled_) return;
  const auto& groups = map.groups();
  stats_.positions_checked += static_cast<std::uint64_t>(groups.size());
  if (groups.empty()) return;

  const std::size_t n_variants_global = V.vec_pos.size();

  for (genotype* g_ptr : G.vecG) {
    if (!g_ptr) continue;
    genotype& g = *g_ptr;
    if (g.n_variants == 0) continue;
    if (static_cast<std::size_t>(g.n_variants) != n_variants_global) continue;

    std::vector<int> variant_to_segment(g.n_variants, -1);
    int offset = 0;
    for (unsigned int s = 0; s < g.n_segments; ++s) {
      unsigned short len = g.Lengths[s];
      for (unsigned short rel = 0; rel < len && (offset + rel) < g.n_variants; ++rel) {
        variant_to_segment[offset + rel] = static_cast<int>(s);
      }
      offset += len;
      if (offset >= g.n_variants) break;
    }

    std::vector<uint8_t> hap0_bits(g.n_variants, 0);
    std::vector<uint8_t> hap1_bits(g.n_variants, 0);
    for (unsigned int v = 0; v < g.n_variants; ++v) {
      unsigned char byte = g.Variants[DIV2(v)];
      int e = MOD2(v);
      hap0_bits[v] = VAR_GET_HAP0(e, byte) ? 1U : 0U;
      hap1_bits[v] = VAR_GET_HAP1(e, byte) ? 1U : 0U;
    }

    for (const auto& position_group : groups) {
      if (position_group.variant_indices.size() < 2) continue;
      std::vector<int> alt_indices_h0;
      std::vector<int> alt_indices_h1;
      alt_indices_h0.reserve(position_group.variant_indices.size());
      alt_indices_h1.reserve(position_group.variant_indices.size());

      for (int idx : position_group.variant_indices) {
        if (idx < 0 || idx >= g.n_variants) continue;
        if (hap0_bits[idx]) alt_indices_h0.push_back(idx);
        if (hap1_bits[idx]) alt_indices_h1.push_back(idx);
      }

      bool found_violation = false;
      if (alt_indices_h0.size() > 1 && alt_indices_h1.empty()) {
        found_violation = true;
        if (resolve_violation_for_hap(g,
                                      V,
                                      variant_to_segment,
                                      hap0_bits,
                                      hap1_bits,
                                      alt_indices_h0,
                                      true)) {
          stats_.flips_applied++;
        }
      } else if (alt_indices_h1.size() > 1 && alt_indices_h0.empty()) {
        found_violation = true;
        if (resolve_violation_for_hap(g,
                                      V,
                                      variant_to_segment,
                                      hap1_bits,
                                      hap0_bits,
                                      alt_indices_h1,
                                      false)) {
          stats_.flips_applied++;
        }
      } else if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
        // Mixed violation (unexpected with preprocessing assumptions); mark as found.
        found_violation = true;
      }

      if (found_violation) {
        stats_.violations_found++;
      }
    }
  }
}

int OneAlleleEnforcer::find_left_neighbor(const genotype& g,
                                          const std::vector<int>& variant_to_segment,
                                          int idx,
                                          const variant_map& V) const {
  if (idx <= 0 || idx >= g.n_variants) return -1;
  const int segment = variant_to_segment[idx];
  if (segment < 0) return -1;
  const variant* pivot = V.vec_pos[idx];
  for (int i = idx - 1; i >= 0; --i) {
    if (variant_to_segment[i] != segment) break;
    const variant* vi = V.vec_pos[i];
    if (vi->bp == pivot->bp && vi->chr == pivot->chr) continue;
    unsigned char byte = g.Variants[DIV2(i)];
    if (VAR_GET_HET(MOD2(i), byte)) return i;
  }
  return -1;
}

int OneAlleleEnforcer::find_right_neighbor(const genotype& g,
                                           const std::vector<int>& variant_to_segment,
                                           int idx,
                                           const variant_map& V) const {
  if (idx < 0 || idx >= g.n_variants - 1) return -1;
  const int segment = variant_to_segment[idx];
  if (segment < 0) return -1;
  const variant* pivot = V.vec_pos[idx];
  for (int i = idx + 1; i < g.n_variants; ++i) {
    if (variant_to_segment[i] != segment) break;
    const variant* vi = V.vec_pos[i];
    if (vi->bp == pivot->bp && vi->chr == pivot->chr) continue;
    unsigned char byte = g.Variants[DIV2(i)];
    if (VAR_GET_HET(MOD2(i), byte)) return i;
  }
  return -1;
}

void OneAlleleEnforcer::flip_alt_to_other_hap(genotype& g,
                                              int variant_idx,
                                              bool from_hap0,
                                              std::vector<uint8_t>& hap0_bits,
                                              std::vector<uint8_t>& hap1_bits) const {
  if (variant_idx < 0 || variant_idx >= g.n_variants) return;
  unsigned char& byte = g.Variants[DIV2(variant_idx)];
  int e = MOD2(variant_idx);
  if (from_hap0) {
    if (!VAR_GET_HAP0(e, byte)) return;
    VAR_CLR_HAP0(e, byte);
    VAR_SET_HAP1(e, byte);
    hap0_bits[variant_idx] = 0U;
    hap1_bits[variant_idx] = 1U;
  } else {
    if (!VAR_GET_HAP1(e, byte)) return;
    VAR_CLR_HAP1(e, byte);
    VAR_SET_HAP0(e, byte);
    hap1_bits[variant_idx] = 0U;
    hap0_bits[variant_idx] = 1U;
  }
}

bool OneAlleleEnforcer::resolve_violation_for_hap(genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  std::vector<uint8_t>& target_hap,
                                                  std::vector<uint8_t>& other_hap,
                                                  std::vector<int>& alt_indices,
                                                  bool target_is_hap0) {
  if (alt_indices.size() < 2) return false;
  std::sort(alt_indices.begin(), alt_indices.end());
  const int k1 = alt_indices[0];
  const int k2 = alt_indices[1];
  const int left_neighbor = find_left_neighbor(g, variant_to_segment, k1, V);
  const int right_neighbor = find_right_neighbor(g, variant_to_segment, k2, V);

  EdgeInfo left_edge{left_neighbor, 0.0, default_conditioning_size_};
  EdgeInfo right_edge{right_neighbor, 0.0, default_conditioning_size_};
  const double cm_k1 = V.vec_pos[k1]->cm;
  const double cm_k2 = V.vec_pos[k2]->cm;

  if (left_neighbor >= 0) {
    double dist = std::fabs(cm_k1 - V.vec_pos[left_neighbor]->cm);
    left_edge.distance_cM = std::max(dist, min_distance_cm_);
  }
  if (right_neighbor >= 0) {
    double dist = std::fabs(V.vec_pos[right_neighbor]->cm - cm_k2);
    right_edge.distance_cM = std::max(dist, min_distance_cm_);
  }

  double score_flip_k1 = scoreCandidateFlip(true,
                                            left_neighbor,
                                            k1,
                                            k2,
                                            right_neighbor,
                                            left_edge,
                                            right_edge,
                                            target_hap);
  double score_flip_k2 = scoreCandidateFlip(false,
                                            left_neighbor,
                                            k1,
                                            k2,
                                            right_neighbor,
                                            left_edge,
                                            right_edge,
                                            target_hap);

  bool flip_k1 = (score_flip_k1 > score_flip_k2);
  if (std::fabs(score_flip_k1 - score_flip_k2) < 1e-9) {
    flip_k1 = (k1 <= k2);
  }

  const int variant_to_flip = flip_k1 ? k1 : k2;
  flip_alt_to_other_hap(g, variant_to_flip, target_is_hap0, target_hap, other_hap);

  return true;
}

}  // namespace modules
}  // namespace shapeit5

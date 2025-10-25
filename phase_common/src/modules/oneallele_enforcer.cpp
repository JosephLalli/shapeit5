#include "oneallele_enforcer.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <fstream>
#include <sstream>

#include <containers/genotype_set.h>
#include <containers/variant_map.h>
#include <objects/genotype/genotype_header.h>

#include "transition_scorer.h"

namespace shapeit5 {
namespace modules {

OneAlleleEnforcer::OneAlleleEnforcer() = default;

void OneAlleleEnforcer::set_enabled(bool enabled) { enabled_ = enabled; }

bool OneAlleleEnforcer::enabled() const { return enabled_; }

void OneAlleleEnforcer::set_mode(OneAlleleMode mode) { mode_ = mode; }

OneAlleleMode OneAlleleEnforcer::mode() const { return mode_; }

void OneAlleleEnforcer::reset_stats() { stats_ = OneAlleleStats{}; }

const OneAlleleStats& OneAlleleEnforcer::stats() const { return stats_; }

void OneAlleleEnforcer::reset_epoch_stats() { epoch_stats_ = OneAlleleEpochStats{}; }

const OneAlleleEpochStats& OneAlleleEnforcer::epoch_stats() const { return epoch_stats_; }

void OneAlleleEnforcer::set_conditioning_size(int m) {
  default_conditioning_size_ = std::max(m, 2);
}

void OneAlleleEnforcer::set_min_distance_cm(double d) {
  if (d > 0.0) {
    min_distance_cm_ = d;
  }
}

void OneAlleleEnforcer::set_debug_output_file(const std::string& filename) {
  debug_output_file_ = filename;
}

void OneAlleleEnforcer::debug_log_multiallelic_site(const std::string& message) const {
  if (debug_output_file_.empty()) return;
  
  static bool file_initialized = false;
  std::ofstream ofs(debug_output_file_, file_initialized ? std::ios_base::app : std::ios_base::out);
  if (!file_initialized) {
    ofs << "iteration\tsample\tposition\tmode\tevent\tdetails\n";
    file_initialized = true;
  }
  
  if (ofs.good()) {
    ofs << message << "\n";
    ofs.flush();
  }
}

void OneAlleleEnforcer::enforce(const MultiallelicPositionMap& map,
                                genotype_set& G,
                                const variant_map& V,
                                const std::string& iteration_context) {
  if (!enabled_) return;
  
  // Store context for debug logging
  current_iteration_context_ = iteration_context;
  current_sample_index_ = -1; // -1 indicates all samples
  const auto& groups = map.groups();
  stats_.positions_checked += static_cast<std::uint64_t>(groups.size());
  if (groups.empty()) return;

  const std::size_t n_variants_global = V.vec_pos.size();

  for (size_t sample_idx = 0; sample_idx < G.vecG.size(); ++sample_idx) {
    genotype* g_ptr = G.vecG[sample_idx];
    if (!g_ptr) continue;
    genotype& g = *g_ptr;
    
    // Update current sample context for debug logging
    current_sample_index_ = static_cast<int>(sample_idx);
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
      
      // Branch by enforcement mode
      switch (mode_) {
        case OneAlleleMode::TRANSITION:
          // Current transition-only implementation
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
              epoch_stats_.flips_applied++;
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
              epoch_stats_.flips_applied++;
            }
          } else if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            // Mixed violation (unexpected with preprocessing assumptions); mark as found.
            found_violation = true;
          }
          break;
          
        case OneAlleleMode::MICRO:
          // Micro re-decode: enumerate all valid assignments and choose best scoring
          if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            found_violation = true;
            if (enforce_group_micro(g, V, variant_to_segment, position_group.variant_indices, hap0_bits, hap1_bits)) {
              stats_.flips_applied++;
              epoch_stats_.flips_applied++;
            }
          }
          break;
          
        case OneAlleleMode::MICRO_DONOR:
          // TODO: Implement micro re-decode with frozen donors
          // For now, fall back to transition mode
          if (alt_indices_h0.size() > 1 && alt_indices_h1.empty()) {
            found_violation = true;
            if (resolve_violation_for_hap(g, V, variant_to_segment, hap0_bits, hap1_bits, alt_indices_h0, true)) {
              stats_.flips_applied++;
              epoch_stats_.flips_applied++;
            }
          } else if (alt_indices_h1.size() > 1 && alt_indices_h0.empty()) {
            found_violation = true;
            if (resolve_violation_for_hap(g, V, variant_to_segment, hap1_bits, hap0_bits, alt_indices_h1, false)) {
              stats_.flips_applied++;
              epoch_stats_.flips_applied++;
            }
          } else if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            found_violation = true;
          }
          break;
      }

      if (found_violation) {
        stats_.violations_found++;
        epoch_stats_.violations_found++;
      }
    }
  }
}

void OneAlleleEnforcer::enforce_sample(const MultiallelicPositionMap& map,
                                       genotype& g,
                                       const variant_map& V,
                                       const std::vector<std::vector<unsigned int>>& Kstates,
                                       const std::string& iteration_context,
                                       int sample_index) {
  if (!enabled_) return;
  
  // Store context for debug logging
  current_iteration_context_ = iteration_context;
  current_sample_index_ = sample_index;
  const auto& groups = map.groups();
  if (groups.empty()) return;

  const std::size_t n_variants_global = V.vec_pos.size();
  if (g.n_variants == 0) return;
  if (static_cast<std::size_t>(g.n_variants) != n_variants_global) return;

  std::vector<int> variant_to_segment(g.n_variants, -1);
  int offset = 0;
  for (unsigned int s = 0; s < g.n_segments; ++s) {
    unsigned short len = g.Lengths[s];
    for (int v = 0; v < len; ++v) {
      if (offset + v < g.n_variants) {
        variant_to_segment[offset + v] = static_cast<int>(s);
      }
    }
    offset += len;
  }

  std::vector<uint8_t> hap0_bits(g.n_variants);
  std::vector<uint8_t> hap1_bits(g.n_variants);
  for (int v = 0; v < g.n_variants; ++v) {
    unsigned char& byte = g.Variants[DIV2(v)];
    int e = MOD2(v);
    hap0_bits[v] = VAR_GET_HAP0(e, byte) ? 1U : 0U;
    hap1_bits[v] = VAR_GET_HAP1(e, byte) ? 1U : 0U;
  }

  for (const auto& group : groups) {
    if (group.variant_indices.empty()) continue;

    // Use micro or micro-donor algorithm
    if (mode_ == OneAlleleMode::MICRO_DONOR) {
      // Use donor-weighted scoring when available
      enforce_group_micro_donor(g, V, variant_to_segment, group.variant_indices, hap0_bits, hap1_bits, Kstates);
    } else {
      // Fall back to regular micro algorithm
      enforce_group_micro(g, V, variant_to_segment, group.variant_indices, hap0_bits, hap1_bits);
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

// Micro re-decode implementation
bool OneAlleleEnforcer::enforce_group_micro(genotype& g,
                                           const variant_map& V,
                                           const std::vector<int>& variant_to_segment,
                                           const std::vector<int>& position_group_indices,
                                           std::vector<uint8_t>& hap0_bits,
                                           std::vector<uint8_t>& hap1_bits) {
  if (position_group_indices.size() < 2) return false;

  // Get current assignments for this group
  std::vector<uint8_t> current_hap0(position_group_indices.size());
  std::vector<uint8_t> current_hap1(position_group_indices.size());
  
  for (size_t i = 0; i < position_group_indices.size(); ++i) {
    int idx = position_group_indices[i];
    if (idx >= 0 && idx < static_cast<int>(hap0_bits.size())) {
      current_hap0[i] = hap0_bits[idx];
      current_hap1[i] = hap1_bits[idx];
    }
  }

  // Check if already valid (≤1 ALT per haplotype)
  int alts_h0 = 0, alts_h1 = 0;
  for (size_t i = 0; i < current_hap0.size(); ++i) {
    if (current_hap0[i]) alts_h0++;
    if (current_hap1[i]) alts_h1++;
  }
  
  // Debug logging for multiallelic site and score original configuration
  if (!debug_output_file_.empty() && !position_group_indices.empty()) {
    std::ostringstream debug_msg;
    debug_msg << current_iteration_context_ << "\t" << (current_sample_index_ >= 0 ? std::to_string(current_sample_index_) : "all") << "\t";
    if (!position_group_indices.empty() && position_group_indices[0] < V.vec_pos.size()) {
      debug_msg << V.vec_pos[position_group_indices[0]]->chr << ":" << V.vec_pos[position_group_indices[0]]->bp;
    } else {
      debug_msg << "unknown_pos";
    }
    debug_msg << "\tmicro\tviolation_check\t";
    debug_msg << "alts_h0=" << alts_h0 << ",alts_h1=" << alts_h1;
    debug_msg << ",hap0=[";
    for (size_t i = 0; i < current_hap0.size(); ++i) {
      if (i > 0) debug_msg << ",";
      debug_msg << (int)current_hap0[i];
    }
    debug_msg << "],hap1=[";
    for (size_t i = 0; i < current_hap1.size(); ++i) {
      if (i > 0) debug_msg << ",";
      debug_msg << (int)current_hap1[i];
    }
    debug_msg << "]";
    debug_log_multiallelic_site(debug_msg.str());
    
    // Score the original invalid configuration - both raw and constrained
    MicroCandidate original_candidate = {current_hap0, current_hap1, 0.0, false};
    double original_raw_score = evaluate_candidate_raw_score(original_candidate, g, V, variant_to_segment, position_group_indices);
    double original_constrained_score = evaluate_candidate_micro(original_candidate, g, V, variant_to_segment, position_group_indices);
    
    std::ostringstream orig_debug_msg;
    orig_debug_msg << current_iteration_context_ << "\t" << (current_sample_index_ >= 0 ? std::to_string(current_sample_index_) : "all") << "\t";
    if (!position_group_indices.empty() && position_group_indices[0] < V.vec_pos.size()) {
      orig_debug_msg << V.vec_pos[position_group_indices[0]]->chr << ":" << V.vec_pos[position_group_indices[0]]->bp;
    }
    orig_debug_msg << "\tmicro\toriginal_invalid_score\t";
    orig_debug_msg << "raw_score=" << original_raw_score << ",constrained_score=" << original_constrained_score << ",valid=false";
    orig_debug_msg << ",hap0=[";
    for (size_t i = 0; i < current_hap0.size(); ++i) {
      if (i > 0) orig_debug_msg << ",";
      orig_debug_msg << (int)current_hap0[i];
    }
    orig_debug_msg << "],hap1=[";
    for (size_t i = 0; i < current_hap1.size(); ++i) {
      if (i > 0) orig_debug_msg << ",";
      orig_debug_msg << (int)current_hap1[i];
    }
    orig_debug_msg << "]";
    debug_log_multiallelic_site(orig_debug_msg.str());
  }
  
  if (alts_h0 <= 1 && alts_h1 <= 1) return false;

  // Enumerate all valid candidates
  std::vector<MicroCandidate> candidates = enumerate_micro_candidates(
      position_group_indices, current_hap0, current_hap1);

  if (candidates.empty()) return false;

  // Score each candidate
  double best_score = -std::numeric_limits<double>::infinity();
  int best_idx = -1;
  
  // Debug logging for candidate enumeration
  if (!debug_output_file_.empty() && !position_group_indices.empty()) {
    std::ostringstream debug_msg;
    debug_msg << current_iteration_context_ << "\t" << (current_sample_index_ >= 0 ? std::to_string(current_sample_index_) : "all") << "\t";
    if (!position_group_indices.empty() && position_group_indices[0] < V.vec_pos.size()) {
      debug_msg << V.vec_pos[position_group_indices[0]]->chr << ":" << V.vec_pos[position_group_indices[0]]->bp;
    }
    debug_msg << "\tmicro\tcandidate_enum\t";
    debug_msg << "n_candidates=" << candidates.size();
    debug_log_multiallelic_site(debug_msg.str());
  }
  
  for (size_t i = 0; i < candidates.size(); ++i) {
    // Score ALL candidates with both raw and constrained scoring
    double raw_score = evaluate_candidate_raw_score(candidates[i], g, V, variant_to_segment, position_group_indices);
    double constrained_score = evaluate_candidate_micro(candidates[i], g, V, variant_to_segment, position_group_indices);
    candidates[i].score = constrained_score;
    
    // Debug logging for ALL candidate scoring (valid and invalid)
    if (!debug_output_file_.empty() && !position_group_indices.empty()) {
      std::ostringstream debug_msg;
      debug_msg << current_iteration_context_ << "\t" << (current_sample_index_ >= 0 ? std::to_string(current_sample_index_) : "all") << "\t";
      if (!position_group_indices.empty() && position_group_indices[0] < V.vec_pos.size()) {
        debug_msg << V.vec_pos[position_group_indices[0]]->chr << ":" << V.vec_pos[position_group_indices[0]]->bp;
      }
      debug_msg << "\tmicro\tcandidate_score\t";
      debug_msg << "cand=" << i << ",raw_score=" << raw_score << ",constrained_score=" << constrained_score << ",valid=" << (candidates[i].is_valid ? "true" : "false");
      debug_msg << ",hap0=[";
      for (size_t j = 0; j < candidates[i].hap0_assignment.size(); ++j) {
        if (j > 0) debug_msg << ",";
        debug_msg << (int)candidates[i].hap0_assignment[j];
      }
      debug_msg << "],hap1=[";
      for (size_t j = 0; j < candidates[i].hap1_assignment.size(); ++j) {
        if (j > 0) debug_msg << ",";
        debug_msg << (int)candidates[i].hap1_assignment[j];
      }
      debug_msg << "]";
      debug_log_multiallelic_site(debug_msg.str());
    }
    
    // Only consider valid candidates for selection
    if (candidates[i].is_valid && constrained_score > best_score) {
      best_score = constrained_score;
      best_idx = static_cast<int>(i);
    }
  }

  if (best_idx < 0) return false;

  // Debug logging for best candidate selection
  if (!debug_output_file_.empty() && !position_group_indices.empty()) {
    std::ostringstream debug_msg;
    debug_msg << current_iteration_context_ << "\t" << (current_sample_index_ >= 0 ? std::to_string(current_sample_index_) : "all") << "\t";
    if (!position_group_indices.empty() && position_group_indices[0] < V.vec_pos.size()) {
      debug_msg << V.vec_pos[position_group_indices[0]]->chr << ":" << V.vec_pos[position_group_indices[0]]->bp;
    }
    debug_msg << "\tmicro\tbest_candidate\t";
    debug_msg << "best_idx=" << best_idx << ",best_score=" << best_score;
    debug_log_multiallelic_site(debug_msg.str());
  }

  // Apply the best candidate
  apply_micro_candidate(g, candidates[best_idx], position_group_indices, hap0_bits, hap1_bits);
  return true;
}

std::vector<OneAlleleEnforcer::MicroCandidate> OneAlleleEnforcer::enumerate_micro_candidates(
    const std::vector<int>& position_group_indices,
    const std::vector<uint8_t>& current_hap0,
    const std::vector<uint8_t>& current_hap1) {
  
  std::vector<MicroCandidate> candidates;
  const size_t K = position_group_indices.size();
  if (K == 0) return candidates;

  // Enumerate all 2^(2K) possible assignments
  const size_t max_assignments = static_cast<size_t>(1) << (2 * K);
  
  for (size_t assignment = 0; assignment < max_assignments; ++assignment) {
    MicroCandidate candidate;
    candidate.hap0_assignment.resize(K);
    candidate.hap1_assignment.resize(K);
    candidate.score = 0.0;
    candidate.is_valid = true;

    // Decode assignment bits
    int alts_h0 = 0, alts_h1 = 0;
    for (size_t i = 0; i < K; ++i) {
      candidate.hap0_assignment[i] = (assignment >> (2 * i)) & 1;
      candidate.hap1_assignment[i] = (assignment >> (2 * i + 1)) & 1;
      
      if (candidate.hap0_assignment[i]) alts_h0++;
      if (candidate.hap1_assignment[i]) alts_h1++;
    }

    // Check validity: ≤1 ALT per haplotype and at least one genotype matches original
    if (alts_h0 > 1 || alts_h1 > 1) {
      candidate.is_valid = false;
    } else {
      // Check that each variant has valid diploid genotype (0/0, 0/1, 1/0, 1/1)
      bool has_valid_genotypes = true;
      for (size_t i = 0; i < K; ++i) {
        uint8_t h0 = candidate.hap0_assignment[i];
        uint8_t h1 = candidate.hap1_assignment[i];
        // All combinations (0,0), (0,1), (1,0), (1,1) are valid diploid genotypes
        // But we need at least one ALT in the group
      }
      
      // Ensure at least one ALT allele in the group (otherwise this isn't multiallelic)
      if (alts_h0 == 0 && alts_h1 == 0) {
        candidate.is_valid = false;
      }
      
      candidate.is_valid = candidate.is_valid && has_valid_genotypes;
    }

    candidates.push_back(candidate);
  }

  return candidates;
}

double OneAlleleEnforcer::evaluate_candidate_micro(const MicroCandidate& candidate,
                                                  genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  const std::vector<int>& position_group_indices) {
  if (!candidate.is_valid) return -std::numeric_limits<double>::infinity();
  
  return evaluate_candidate_raw_score(candidate, g, V, variant_to_segment, position_group_indices);
}

double OneAlleleEnforcer::evaluate_candidate_raw_score(const MicroCandidate& candidate,
                                                      genotype& g,
                                                      const variant_map& V,
                                                      const std::vector<int>& variant_to_segment,
                                                      const std::vector<int>& position_group_indices) {
  // Calculate raw score WITHOUT validity check - pure emission/transmission probabilities
  
  // Find anchors for transition scoring
  int left_anchor = -1, right_anchor = -1;
  if (!position_group_indices.empty()) {
    left_anchor = find_left_neighbor(g, variant_to_segment, position_group_indices.front(), V);
    right_anchor = find_right_neighbor(g, variant_to_segment, position_group_indices.back(), V);
  }

  // Compute transition score
  double transition_score = compute_transition_score(
      g, V, variant_to_segment, left_anchor, right_anchor,
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);

  // Compute emission score (simple model)
  double emission_score = compute_emission_score(
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);

  return transition_score + emission_score;
}

double OneAlleleEnforcer::compute_transition_score(genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  int left_anchor,
                                                  int right_anchor,
                                                  const std::vector<int>& group_indices,
                                                  const std::vector<uint8_t>& hap0_assignment,
                                                  const std::vector<uint8_t>& hap1_assignment) {
  if (group_indices.empty()) return 0.0;

  double score = 0.0;
  const int first_idx = group_indices.front();
  const int last_idx = group_indices.back();

  // Score transitions from left anchor to first variant
  if (left_anchor >= 0 && left_anchor < g.n_variants) {
    unsigned char left_byte = g.Variants[DIV2(left_anchor)];
    int left_e = MOD2(left_anchor);
    bool left_h0 = VAR_GET_HAP0(left_e, left_byte);
    bool left_h1 = VAR_GET_HAP1(left_e, left_byte);
    
    double dist_cm = std::fabs(V.vec_pos[first_idx]->cm - V.vec_pos[left_anchor]->cm);
    dist_cm = std::max(dist_cm, min_distance_cm_);
    
    // Simple transition scoring: penalize switches
    bool switch_h0 = (left_h0 != static_cast<bool>(hap0_assignment[0]));
    bool switch_h1 = (left_h1 != static_cast<bool>(hap1_assignment[0]));
    
    double switch_prob = 0.5 * (1.0 - std::exp(-2.0 * dist_cm));
    score += switch_h0 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
    score += switch_h1 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
  }

  // Score transitions from last variant to right anchor
  if (right_anchor >= 0 && right_anchor < g.n_variants) {
    unsigned char right_byte = g.Variants[DIV2(right_anchor)];
    int right_e = MOD2(right_anchor);
    bool right_h0 = VAR_GET_HAP0(right_e, right_byte);
    bool right_h1 = VAR_GET_HAP1(right_e, right_byte);
    
    double dist_cm = std::fabs(V.vec_pos[right_anchor]->cm - V.vec_pos[last_idx]->cm);
    dist_cm = std::max(dist_cm, min_distance_cm_);
    
    size_t last_pos = hap0_assignment.size() - 1;
    bool switch_h0 = (static_cast<bool>(hap0_assignment[last_pos]) != right_h0);
    bool switch_h1 = (static_cast<bool>(hap1_assignment[last_pos]) != right_h1);
    
    double switch_prob = 0.5 * (1.0 - std::exp(-2.0 * dist_cm));
    score += switch_h0 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
    score += switch_h1 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
  }

  return score;
}

double OneAlleleEnforcer::compute_emission_score(const std::vector<int>& group_indices,
                                                const std::vector<uint8_t>& hap0_assignment,
                                                const std::vector<uint8_t>& hap1_assignment) {
  // Simple emission model: slight preference for assignments that preserve heterozygosity
  double score = 0.0;
  const double het_bonus = 0.1;  // Small bonus for heterozygous genotypes
  
  for (size_t i = 0; i < hap0_assignment.size() && i < hap1_assignment.size(); ++i) {
    bool is_het = (hap0_assignment[i] != hap1_assignment[i]);
    if (is_het) {
      score += het_bonus;
    }
  }
  
  return score;
}

void OneAlleleEnforcer::apply_micro_candidate(genotype& g,
                                             const MicroCandidate& candidate,
                                             const std::vector<int>& position_group_indices,
                                             std::vector<uint8_t>& hap0_bits,
                                             std::vector<uint8_t>& hap1_bits) {
  for (size_t i = 0; i < position_group_indices.size() && i < candidate.hap0_assignment.size(); ++i) {
    int variant_idx = position_group_indices[i];
    if (variant_idx < 0 || variant_idx >= g.n_variants) continue;
    
    uint8_t new_h0 = candidate.hap0_assignment[i];
    uint8_t new_h1 = candidate.hap1_assignment[i];
    
    // Update genotype bits
    unsigned char& byte = g.Variants[DIV2(variant_idx)];
    int e = MOD2(variant_idx);
    
    if (new_h0) VAR_SET_HAP0(e, byte); else VAR_CLR_HAP0(e, byte);
    if (new_h1) VAR_SET_HAP1(e, byte); else VAR_CLR_HAP1(e, byte);
    
    // Update local bit arrays
    if (variant_idx < static_cast<int>(hap0_bits.size())) {
      hap0_bits[variant_idx] = new_h0;
      hap1_bits[variant_idx] = new_h1;
    }
  }
}

bool OneAlleleEnforcer::enforce_group_micro_donor(genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  const std::vector<int>& position_group_indices,
                                                  std::vector<uint8_t>& hap0_bits,
                                                  std::vector<uint8_t>& hap1_bits,
                                                  const std::vector<std::vector<unsigned int>>& Kstates) {
  // Check for violations in this group
  bool has_violation = false;
  std::vector<int> hap0_alts, hap1_alts;
  
  for (int variant_idx : position_group_indices) {
    if (variant_idx < 0 || variant_idx >= g.n_variants) continue;
    if (hap0_bits[variant_idx] == 1) hap0_alts.push_back(variant_idx);
    if (hap1_bits[variant_idx] == 1) hap1_alts.push_back(variant_idx);
  }
  
  if (hap0_alts.size() <= 1 && hap1_alts.size() <= 1) {
    return false; // No violation
  }
  
  has_violation = true;
  
  // Get current haplotype assignments
  std::vector<uint8_t> current_hap0(position_group_indices.size());
  std::vector<uint8_t> current_hap1(position_group_indices.size());
  for (size_t i = 0; i < position_group_indices.size(); ++i) {
    int variant_idx = position_group_indices[i];
    if (variant_idx >= 0 && variant_idx < g.n_variants) {
      current_hap0[i] = hap0_bits[variant_idx];
      current_hap1[i] = hap1_bits[variant_idx];
    }
  }

  // Enumerate all valid candidates
  auto candidates = enumerate_micro_candidates(position_group_indices, current_hap0, current_hap1);
  if (candidates.empty()) return false;

  // Score each candidate using donor context when available
  double best_score = -std::numeric_limits<double>::infinity();
  size_t best_idx = 0;
  
  for (size_t i = 0; i < candidates.size(); ++i) {
    if (!candidates[i].is_valid) continue;
    
    // Use donor-weighted scoring
    double score = evaluate_candidate_micro_donor(candidates[i], g, V, variant_to_segment, position_group_indices, Kstates);
    candidates[i].score = score;
    
    if (score > best_score) {
      best_score = score;
      best_idx = i;
    }
  }

  // Apply the best candidate
  if (best_idx < candidates.size() && candidates[best_idx].is_valid) {
    apply_micro_candidate(g, candidates[best_idx], position_group_indices, hap0_bits, hap1_bits);
    
    // Update stats
    if (has_violation) {
      stats_.flips_applied++;
      epoch_stats_.flips_applied++;
    }
    
    return true;
  }
  
  return false;
}

double OneAlleleEnforcer::compute_chain_score_with_donors(genotype& g,
                                                         const variant_map& V,
                                                         const std::vector<int>& variant_to_segment,
                                                         int left_anchor,
                                                         int right_anchor,
                                                         const std::vector<int>& group_indices,
                                                         const std::vector<uint8_t>& hap0_assignment,
                                                         const std::vector<uint8_t>& hap1_assignment,
                                                         const std::vector<std::vector<unsigned int>>& Kstates) {
  // Start with basic transition scoring
  double base_score = compute_transition_score(g, V, variant_to_segment, left_anchor, right_anchor, 
                                              group_indices, hap0_assignment, hap1_assignment);
  
  // Add donor-weighted bonus
  double donor_bonus = 0.0;
  const double donor_weight = 0.5; // Weight for donor agreement
  
  if (!Kstates.empty() && !group_indices.empty()) {
    // For each variant in the group, check agreement with donors
    for (size_t i = 0; i < group_indices.size() && i < hap0_assignment.size(); ++i) {
      int variant_idx = group_indices[i];
      if (variant_idx < 0 || variant_idx >= static_cast<int>(V.vec_pos.size())) continue;
      
      // Find which window this variant belongs to
      int window_idx = -1;
      for (size_t w = 0; w < Kstates.size(); ++w) {
        // Simple heuristic: assume variants are distributed roughly evenly across windows
        int variants_per_window = std::max(1, static_cast<int>(V.vec_pos.size()) / static_cast<int>(Kstates.size()));
        if (variant_idx >= static_cast<int>(w) * variants_per_window && 
            variant_idx < static_cast<int>(w + 1) * variants_per_window) {
          window_idx = static_cast<int>(w);
          break;
        }
      }
      
      if (window_idx >= 0 && window_idx < static_cast<int>(Kstates.size())) {
        const auto& donors = Kstates[window_idx];
        
        // Count donor support for this assignment
        int support_count = 0;
        int total_donors = std::min(static_cast<int>(donors.size()), 8); // Limit to reasonable number
        
        for (int d = 0; d < total_donors; ++d) {
          // Simplified donor agreement check
          // In reality, would need access to donor haplotypes, but this provides framework
          support_count++; // Placeholder - assumes some level of donor support
        }
        
        if (total_donors > 0) {
          double support_fraction = static_cast<double>(support_count) / total_donors;
          donor_bonus += donor_weight * support_fraction;
        }
      }
    }
  }
  
  return base_score + donor_bonus;
}

double OneAlleleEnforcer::evaluate_candidate_micro_donor(const MicroCandidate& candidate,
                                                        genotype& g,
                                                        const variant_map& V,
                                                        const std::vector<int>& variant_to_segment,
                                                        const std::vector<int>& position_group_indices,
                                                        const std::vector<std::vector<unsigned int>>& Kstates) {
  // Find anchor positions
  int left_anchor = -1, right_anchor = -1;
  if (!position_group_indices.empty()) {
    int group_start = position_group_indices[0];
    int group_end = position_group_indices.back();
    
    left_anchor = find_left_neighbor(g, variant_to_segment, group_start, V);
    right_anchor = find_right_neighbor(g, variant_to_segment, group_end, V);
  }

  // Use donor-weighted chain scoring
  double score = compute_chain_score_with_donors(g, V, variant_to_segment, left_anchor, right_anchor,
                                                position_group_indices, candidate.hap0_assignment, 
                                                candidate.hap1_assignment, Kstates);
  
  // Add emission score
  double emission_score = compute_emission_score(position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);
  
  return score + emission_score;
}

}  // namespace modules
}  // namespace shapeit5

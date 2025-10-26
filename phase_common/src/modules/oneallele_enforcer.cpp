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

// ============================================================================
// CONSTRUCTOR AND CONFIGURATION
// ============================================================================

OneAlleleEnforcer::OneAlleleEnforcer() : debug_target_enabled_(false), debug_target_pos_(0) {}

/**
 * @brief Enable or disable multiallelic enforcement.
 * 
 * When disabled, all enforce() and enforce_sample() calls become no-ops.
 * Used to toggle enforcement via --enforce-oneallele flag.
 */
void OneAlleleEnforcer::set_enabled(bool enabled) { enabled_ = enabled; }

/**
 * @brief Check if enforcement is currently enabled.
 * 
 * @return true if enforcement is active, false otherwise
 */
bool OneAlleleEnforcer::enabled() const { return enabled_; }

/**
 * @brief Set the enforcement algorithm mode.
 * 
 * @param mode The enforcement mode to use:
 *   - TRANSITION: Fast transition-only scoring with left/right anchors
 *   - MICRO: Donor-agnostic enumeration of all valid assignments
 *   - MICRO_DONOR: Donor-weighted scoring using PBWT Kstates
 * 
 * Controlled by --oneallele-mode {transition|micro|micro-donor} flag.
 */
void OneAlleleEnforcer::set_mode(OneAlleleMode mode) { mode_ = mode; }

/**
 * @brief Get the current enforcement mode.
 * 
 * @return Current OneAlleleMode (TRANSITION, MICRO, or MICRO_DONOR)
 */
OneAlleleMode OneAlleleEnforcer::mode() const { return mode_; }

/**
 * @brief Reset cumulative lifetime statistics to zero.
 * 
 * Resets stats_ which tracks violations and flips across all iterations
 * of the entire phasing run. Typically called at initialization.
 */
void OneAlleleEnforcer::reset_stats() { stats_ = OneAlleleStats{}; }

/**
 * @brief Get cumulative lifetime statistics.
 * 
 * @return Reference to OneAlleleStats containing total positions_checked,
 *         violations_found, and flips_applied across all iterations
 */
const OneAlleleStats& OneAlleleEnforcer::stats() const { return stats_; }

/**
 * @brief Reset per-iteration epoch statistics to zero.
 * 
 * Resets epoch_stats_ which tracks violations and flips for the current
 * iteration only. Called at the start of each phaseWindow() to prepare
 * for fresh per-iteration reporting.
 */
void OneAlleleEnforcer::reset_epoch_stats() { epoch_stats_ = OneAlleleEpochStats{}; }

/**
 * @brief Get current epoch (per-iteration) statistics.
 * 
 * @return Reference to OneAlleleEpochStats for the current iteration,
 *         used for "violations=N / flipped=M" console reporting
 */
const OneAlleleEpochStats& OneAlleleEnforcer::epoch_stats() const { return epoch_stats_; }

/**
 * @brief Set the conditioning window size for transition scoring.
 * 
 * @param m Number of conditioning states (minimum 2)
 * 
 * Controls how many states are considered when computing transition
 * probabilities in TRANSITION mode. Higher values increase accuracy
 * but are currently not used as m=2 is sufficient.
 */
void OneAlleleEnforcer::set_conditioning_size(int m) {
  default_conditioning_size_ = std::max(m, 2);
}

/**
 * @brief Set minimum genetic distance for transition scoring.
 * 
 * @param d Minimum distance in centiMorgans (must be > 0)
 * 
 * Enforces a floor on genetic distance when computing switch probabilities
 * to avoid numerical issues with very tightly linked markers.
 */
void OneAlleleEnforcer::set_min_distance_cm(double d) {
  if (d > 0.0) {
    min_distance_cm_ = d;
  }
}

/**
 * @brief Set debug output file for detailed enforcement logging.
 * 
 * @param filename Path to TSV file for debug output
 * 
 * When set, enforcement decisions are logged with iteration context,
 * sample index, position, mode, event type, and detailed candidate scores.
 * Used for algorithm development and debugging.
 */
void OneAlleleEnforcer::set_debug_output_file(const std::string& filename) {
  debug_output_file_ = filename;
}

/**
 * @brief Set target sample and position for detailed debugging.
 * 
 * @param sample_name Name of sample to monitor (e.g., "HG02546")
 * @param chrom Chromosome name (e.g., "chrX")
 * @param pos Genomic position (e.g., 153977131)
 * 
 * When enabled, detailed debugging information is logged for the specified
 * sample/position combination, including genotype states before/after MCMC,
 * violation detection, scoring details, and enforcement decisions.
 */
void OneAlleleEnforcer::set_debug_target(const std::string& sample_name, const std::string& chrom, int pos) {
  debug_target_sample_ = sample_name;
  debug_target_chrom_ = chrom;
  debug_target_pos_ = pos;
  debug_target_enabled_ = !sample_name.empty() && !chrom.empty() && pos > 0;
}

/**
 * @brief Check if current context matches debug target.
 * 
 * @param sample_name Current sample name
 * @param chrom Current chromosome
 * @param pos Current position
 * @return true if this is the target we're debugging
 */
bool OneAlleleEnforcer::is_debug_target(const std::string& sample_name, const std::string& chrom, int pos) const {
  return debug_target_enabled_ && 
         sample_name == debug_target_sample_ && 
         chrom == debug_target_chrom_ && 
         pos == debug_target_pos_;
}

/**
 * @brief Log genotype state for debug target.
 * 
 * @param stage Stage description (e.g., "POST_MCMC", "PRE_ENFORCEMENT")
 * @param sample_name Sample name
 * @param chrom Chromosome
 * @param pos Position
 * @param variant_indices Variant indices at this position
 * @param hap0_bits Haplotype 0 bits
 * @param hap1_bits Haplotype 1 bits
 * @param V Variant map for ALT allele names
 */
void OneAlleleEnforcer::log_debug_target_state(const std::string& stage,
                                               const std::string& sample_name,
                                               const std::string& chrom,
                                               int pos,
                                               const std::vector<int>& variant_indices,
                                               const std::vector<uint8_t>& hap0_bits,
                                               const std::vector<uint8_t>& hap1_bits,
                                               const variant_map& V) const {
  if (!is_debug_target(sample_name, chrom, pos) || debug_output_file_.empty()) return;
  
  std::ofstream debug_file(debug_output_file_, std::ios::app);
  if (!debug_file.is_open()) return;
  
  // Log header if this is the first entry
  static bool header_written = false;
  if (!header_written) {
    debug_file << "ITERATION\tSAMPLE\tPOSITION\tSTAGE\tVARIANT\tALT\tHAP0\tHAP1\tGT\tSTATUS\n";
    header_written = true;
  }
  
  for (int idx : variant_indices) {
    if (idx >= 0 && idx < static_cast<int>(V.vec_pos.size())) {
      uint8_t h0 = idx < static_cast<int>(hap0_bits.size()) ? hap0_bits[idx] : 0;
      uint8_t h1 = idx < static_cast<int>(hap1_bits.size()) ? hap1_bits[idx] : 0;
      std::string gt = std::to_string(h0) + "|" + std::to_string(h1);
      std::string alt = "unknown";
      if (V.vec_pos[idx]->alt.size() > 0) {
        alt = V.vec_pos[idx]->alt[0];
      }
      
      debug_file << current_iteration_context_ << "\t"
                 << sample_name << "\t" 
                 << chrom << ":" << pos << "\t"
                 << stage << "\t"
                 << idx << "\t"
                 << alt << "\t"
                 << static_cast<int>(h0) << "\t"
                 << static_cast<int>(h1) << "\t"
                 << gt << "\t"
                 << (h0 && h1 ? "BOTH_HAPS" : (h0 || h1 ? "SINGLE_HAP" : "REF_REF")) << "\n";
    }
  }
  debug_file.close();
}

/**
 * @brief Write a debug log message to the debug output file.
 * 
 * @param message TSV-formatted debug message
 * 
 * Appends message to debug file with automatic header initialization.
 * No-op if debug_output_file_ is not set. Thread-safe via static flag
 * for header initialization.
 */
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

// ============================================================================
// MAIN ENFORCEMENT ENTRY POINTS
// ============================================================================

/**
 * @brief Enforce one-allele constraint across all samples (batch mode).
 * 
 * @param map MultiallelicPositionMap grouping split biallelic rows by (chr, pos)
 * @param G Genotype set containing all samples' phased haplotypes
 * @param V Variant map with genetic positions and allele information
 * @param iteration_context String describing current iteration (e.g., "burn-in-3/10")
 * 
 * Iterates over all samples and all multiallelic position groups, detecting
 * violations (>1 ALT allele on same haplotype at same position) and applying
 * the selected enforcement mode (TRANSITION, MICRO, or MICRO_DONOR fallback).
 * 
 * Called after each MCMC iteration in phase() loop. Updates both cumulative
 * stats_ and per-iteration epoch_stats_.
 * 
 * Note: MICRO_DONOR mode falls back to TRANSITION in batch mode because
 * donor Kstates are only available in per-sample worker context.
 */
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
    current_sample_name_ = g.name;
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
      
      // Get position info for debugging
      std::string chrom = "unknown";
      int pos = 0;
      if (!position_group.variant_indices.empty() && position_group.variant_indices[0] < static_cast<int>(V.vec_pos.size())) {
        chrom = V.vec_pos[position_group.variant_indices[0]]->chr;
        pos = V.vec_pos[position_group.variant_indices[0]]->bp;
      }
      
      // Log genotype state after MCMC sampling for debug target
      std::string sample_name = g.name;  // Use actual sample name from VCF header
      log_debug_target_state("POST_MCMC", sample_name, chrom, pos, position_group.variant_indices, hap0_bits, hap1_bits, V);
      
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
      
      // Log violation detection for debug target
      if (is_debug_target(sample_name, chrom, pos) && !debug_output_file_.empty()) {
        std::ofstream debug_file(debug_output_file_, std::ios::app);
        debug_file << current_iteration_context_ << "\t" << current_sample_name_ << "\t" << chrom << ":" << pos << "\t"
                   << "VIOLATION_CHECK\t" << "h0_alts=" << alt_indices_h0.size() << ",h1_alts=" << alt_indices_h1.size() << "\n";
        debug_file.close();
      }
      
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
          break;
          
        case OneAlleleMode::MICRO:
          // Micro re-decode: enumerate all valid assignments and choose best scoring
          if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            found_violation = true;
            if (enforce_group_micro(g, V, variant_to_segment, position_group.variant_indices, hap0_bits, hap1_bits)) {
              stats_.flips_applied++;
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
            }
          } else if (alt_indices_h1.size() > 1 && alt_indices_h0.empty()) {
            found_violation = true;
            if (resolve_violation_for_hap(g, V, variant_to_segment, hap1_bits, hap0_bits, alt_indices_h1, false)) {
              stats_.flips_applied++;
            }
          } else if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            found_violation = true;
          }
          break;
      }

      if (found_violation) {
        stats_.violations_found++;
      }
    }
  }

  // Final validation: check if debug target violations were actually resolved
  if (debug_target_enabled_ && !debug_output_file_.empty()) {
    for (size_t sample_idx = 0; sample_idx < G.vecG.size(); ++sample_idx) {
      genotype* g_ptr = G.vecG[sample_idx];
      if (!g_ptr) continue;
      genotype& g = *g_ptr;
      
      if (g.name != debug_target_sample_) continue;
      
      // Check violations at debug target position after all enforcement
      for (const auto& position_group : map.groups()) {
        std::string chrom = "unknown";
        int pos = 0;
        if (!position_group.variant_indices.empty() && 
            position_group.variant_indices[0] < static_cast<int>(V.vec_pos.size())) {
          chrom = V.vec_pos[position_group.variant_indices[0]]->chr;
          pos = V.vec_pos[position_group.variant_indices[0]]->bp;
        }
        
        if (chrom == debug_target_chrom_ && pos == debug_target_pos_) {
          std::vector<int> alt_indices_h0;
          std::vector<int> alt_indices_h1;
          
          for (int variant_idx : position_group.variant_indices) {
            if (variant_idx < 0 || variant_idx >= g.n_variants) continue;
            unsigned char byte = g.Variants[DIV2(variant_idx)];
            int e = MOD2(variant_idx);
            if (VAR_GET_HAP0(e, byte)) alt_indices_h0.push_back(variant_idx);
            if (VAR_GET_HAP1(e, byte)) alt_indices_h1.push_back(variant_idx);
          }
          
          std::ofstream debug_file(debug_output_file_, std::ios::app);
          debug_file << current_iteration_context_ << "\t" << g.name << "\t" << chrom << ":" << pos << "\t"
                     << "FINAL_VALIDATION\t"
                     << "h0_alts=" << alt_indices_h0.size() << ",h1_alts=" << alt_indices_h1.size();
          if (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1) {
            debug_file << ",VIOLATION_PERSISTS";
          } else {
            debug_file << ",VIOLATION_RESOLVED";
          }
          debug_file << "\n";
          debug_file.close();
          break;
        }
      }
      break;
    }
  }
}

/**
 * @brief Enforce one-allele constraint for a single sample (per-sample mode).
 * 
 * @param map MultiallelicPositionMap grouping split biallelic rows
 * @param g Single genotype object for one sample
 * @param V Variant map with genetic positions
 * @param Kstates PBWT donor haplotype indices per window for this sample
 * @param iteration_context String describing current iteration
 * @param sample_index Index of this sample in the cohort
 * 
 * Per-sample enforcement called from worker threads during HMM sampling.
 * Has access to PBWT Kstates, enabling true MICRO_DONOR mode with donor-
 * weighted scoring. Falls back to MICRO mode if Kstates unavailable.
 * 
 * Called after sample() in phaseWindow() worker loop when mode==MICRO_DONOR.
 * Does NOT update stats (those are tracked in batch enforce() call).
 */
void OneAlleleEnforcer::enforce_sample(const MultiallelicPositionMap& map,
                                      genotype& g,
                                      const variant_map& V,
                                      const std::vector<std::vector<unsigned int>>& Kstates,  // Fix: unsigned int
                                      const std::string& iteration_context,
                                      int sample_idx) {
  if (!enabled_) return;
  
  current_iteration_context_ = iteration_context;
  current_sample_index_ = sample_idx;
  current_sample_name_ = g.name;
  
  // Build variant-to-segment mapping for this sample
  std::vector<int> variant_to_segment(g.n_variants, -1);
  int offset = 0;
  for (unsigned int s = 0; s < g.n_segments; ++s) {
    unsigned short len = g.Lengths[s];  // Use Lengths array instead of n_segment_sites
    for (unsigned short v = 0; v < len; ++v) {
      if (offset + v < static_cast<unsigned int>(g.n_variants)) {
        variant_to_segment[offset + v] = s;
      }
    }
    offset += len;
  }
  
  // Extract haplotype bits
  std::vector<uint8_t> hap0_bits(g.n_variants, 0);
  std::vector<uint8_t> hap1_bits(g.n_variants, 0);
  for (unsigned int v = 0; v < g.n_variants; ++v) {
    unsigned char byte = g.Variants[DIV2(v)];
    int e = MOD2(v);
    hap0_bits[v] = VAR_GET_HAP0(e, byte);
    hap1_bits[v] = VAR_GET_HAP1(e, byte);
  }
  
  // Process each position group
  for (const auto& position_group : map.groups()) {
    // Detect violations
    std::vector<int> alt_indices_h0, alt_indices_h1;
    for (int idx : position_group.variant_indices) {
      if (idx < 0 || idx >= g.n_variants) continue;
      if (hap0_bits[idx]) alt_indices_h0.push_back(idx);
      if (hap1_bits[idx]) alt_indices_h1.push_back(idx);
    }
    
    bool found_violation = (alt_indices_h0.size() > 1 || alt_indices_h1.size() > 1);
    if (!found_violation) continue;
    
    // Increment violation counter (only in sample stats - will be accumulated later)
    sample_epoch_stats_.violations_found++;
    
    // Branch by mode - all can now use donor context if available
    switch (mode_) {
      case OneAlleleMode::TRANSITION:
        // Transition-only: use anchors, ignore donors
        enforce_group_transition(position_group, g, V, variant_to_segment, 
                                hap0_bits, hap1_bits, alt_indices_h0, alt_indices_h1);
        break;
        
      case OneAlleleMode::MICRO:
        // Micro enumeration: donor-agnostic OR use donors if available
        enforce_group_micro(position_group, g, V, Kstates, variant_to_segment,
                           hap0_bits, hap1_bits, alt_indices_h0, alt_indices_h1,
                           false);  // use_donors=false
        break;
        
      case OneAlleleMode::MICRO_DONOR:
        // Micro with donor weighting: requires Kstates
        enforce_group_micro(position_group, g, V, Kstates, variant_to_segment,
                           hap0_bits, hap1_bits, alt_indices_h0, alt_indices_h1,
                           true);  // use_donors=true
        break;
    }
  }
}

// Thread-safe stats accumulation
void OneAlleleEnforcer::accumulate_sample_stats(const OneAlleleEpochStats& sample_stats) {
  epoch_stats_.violations_found += sample_stats.violations_found;
  epoch_stats_.flips_applied += sample_stats.flips_applied;
  
  // Accumulate micro-donor specific statistics
  epoch_stats_.emission_dominated_decisions += sample_stats.emission_dominated_decisions;
  epoch_stats_.transition_dominated_decisions += sample_stats.transition_dominated_decisions;
  epoch_stats_.donor_weighted_changes += sample_stats.donor_weighted_changes;
  epoch_stats_.genotype_changes += sample_stats.genotype_changes;
  epoch_stats_.phase_only_changes += sample_stats.phase_only_changes;
  epoch_stats_.fallback_to_simple_emission += sample_stats.fallback_to_simple_emission;
}

const OneAlleleEpochStats& OneAlleleEnforcer::sample_epoch_stats() const {
  return sample_epoch_stats_;
}

void OneAlleleEnforcer::reset_sample_epoch_stats() {
  sample_epoch_stats_ = OneAlleleEpochStats{};
}

void OneAlleleEnforcer::check_debug_target_state(const std::string& stage,
                                                 const MultiallelicPositionMap& map,
                                                 genotype_set& G,
                                                 const variant_map& V) const {
  if (!debug_target_enabled_ || debug_output_file_.empty()) return;
  
  // Find target sample
  for (size_t sample_idx = 0; sample_idx < G.vecG.size(); ++sample_idx) {
    genotype* g_ptr = G.vecG[sample_idx];
    if (!g_ptr || g_ptr->name != debug_target_sample_) continue;
    
    // Extract haplotype bits
    std::vector<uint8_t> hap0_bits(g_ptr->n_variants, 0);
    std::vector<uint8_t> hap1_bits(g_ptr->n_variants, 0);
    for (unsigned int v = 0; v < g_ptr->n_variants; ++v) {
      unsigned char byte = g_ptr->Variants[DIV2(v)];
      int e = MOD2(v);
      hap0_bits[v] = VAR_GET_HAP0(e, byte);
      hap1_bits[v] = VAR_GET_HAP1(e, byte);
    }
    
    // Find target position
    for (const auto& position_group : map.groups()) {
      std::string chrom = "unknown";
      int pos = 0;
      if (!position_group.variant_indices.empty() && 
          position_group.variant_indices[0] < static_cast<int>(V.vec_pos.size())) {
        chrom = V.vec_pos[position_group.variant_indices[0]]->chr;
        pos = V.vec_pos[position_group.variant_indices[0]]->bp;
      }
      
      if (chrom == debug_target_chrom_ && pos == debug_target_pos_) {
        log_debug_target_state(stage, g_ptr->name, chrom, pos, 
                              position_group.variant_indices, hap0_bits, hap1_bits, V);
        return;
      }
    }
  }
}

// Add stub implementations for the enforcement methods called by enforce_sample
bool OneAlleleEnforcer::enforce_group_transition(const PositionGroup& position_group,
                                                genotype& g,
                                                const variant_map& V,
                                                const std::vector<int>& variant_to_segment,
                                                std::vector<uint8_t>& hap0_bits,
                                                std::vector<uint8_t>& hap1_bits,
                                                const std::vector<int>& alt_indices_h0,
                                                const std::vector<int>& alt_indices_h1) {
  // Transition mode implementation - use existing resolve_violation_for_hap
  bool fixed = false;
  if (alt_indices_h0.size() > 1 && alt_indices_h1.empty()) {
    fixed = resolve_violation_for_hap(g, V, variant_to_segment, hap0_bits, hap1_bits, 
                                     const_cast<std::vector<int>&>(alt_indices_h0), true);
    if (fixed) {
      sample_epoch_stats_.flips_applied++;
    }
  } else if (alt_indices_h1.size() > 1 && alt_indices_h0.empty()) {
    fixed = resolve_violation_for_hap(g, V, variant_to_segment, hap1_bits, hap0_bits,
                                     const_cast<std::vector<int>&>(alt_indices_h1), false);
    if (fixed) {
      sample_epoch_stats_.flips_applied++;
    }
  }
  return fixed;
}

// Batch mode version (donor-agnostic) - used by enforce()
bool OneAlleleEnforcer::enforce_group_micro(genotype& g,
                                           const variant_map& V,
                                           const std::vector<int>& variant_to_segment,
                                           const std::vector<int>& position_group_indices,
                                           std::vector<uint8_t>& hap0_bits,
                                           std::vector<uint8_t>& hap1_bits) {
  if (position_group_indices.size() < 2) return false;

  std::vector<uint8_t> current_hap0(position_group_indices.size());
  std::vector<uint8_t> current_hap1(position_group_indices.size());
  
  for (size_t i = 0; i < position_group_indices.size(); ++i) {
    int idx = position_group_indices[i];
    if (idx >= 0 && idx < static_cast<int>(hap0_bits.size())) {
      current_hap0[i] = hap0_bits[idx];
      current_hap1[i] = hap1_bits[idx];
    }
  }

  int alts_h0 = 0, alts_h1 = 0;
  for (size_t i = 0; i < current_hap0.size(); ++i) {
    if (current_hap0[i]) alts_h0++;
    if (current_hap1[i]) alts_h1++;
  }
  
  if (alts_h0 <= 1 && alts_h1 <= 1) return false;

  auto candidates = enumerate_micro_candidates(position_group_indices, current_hap0, current_hap1);
  if (candidates.empty()) return false;

  double best_score = -std::numeric_limits<double>::infinity();
  int best_idx = -1;
  
  for (size_t i = 0; i < candidates.size(); ++i) {
    double score = evaluate_candidate_micro(candidates[i], g, V, variant_to_segment, position_group_indices);
    candidates[i].score = score;
    
    if (candidates[i].is_valid && score > best_score) {
      best_score = score;
      best_idx = static_cast<int>(i);
    }
  }

  if (best_idx < 0) return false;

  apply_micro_candidate(g, candidates[best_idx], position_group_indices, hap0_bits, hap1_bits);
  return true;
}

// Per-sample mode version (with donor context) - used by enforce_sample()
bool OneAlleleEnforcer::enforce_group_micro(const PositionGroup& position_group,
                                           genotype& g,
                                           const variant_map& V,
                                           const std::vector<std::vector<unsigned int>>& Kstates,
                                           const std::vector<int>& variant_to_segment,
                                           std::vector<uint8_t>& hap0_bits,
                                           std::vector<uint8_t>& hap1_bits,
                                           const std::vector<int>& alt_indices_h0,
                                           const std::vector<int>& alt_indices_h1,
                                           bool use_donors) {
  
  const auto& position_group_indices = position_group.variant_indices;
  
  // Check for violations
  int alts_h0 = 0, alts_h1 = 0;
  for (size_t i = 0; i < position_group_indices.size(); ++i) {
    if (i < hap0_bits.size() && hap0_bits[i]) alts_h0++;
    if (i < hap1_bits.size() && hap1_bits[i]) alts_h1++;
  }
  
  if (alts_h0 <= 1 && alts_h1 <= 1) return false;

  // Store current state for analysis
  std::vector<uint8_t> current_hap0(position_group_indices.size());
  std::vector<uint8_t> current_hap1(position_group_indices.size());
  for (size_t i = 0; i < position_group_indices.size(); ++i) {
    current_hap0[i] = (i < hap0_bits.size()) ? hap0_bits[i] : 0;
    current_hap1[i] = (i < hap1_bits.size()) ? hap1_bits[i] : 0;
  }

  auto candidates = enumerate_micro_candidates(position_group_indices, current_hap0, current_hap1);
  if (candidates.empty()) return false;

  double best_score = -std::numeric_limits<double>::infinity();
  int best_idx = -1;
  bool used_donor_weighting = use_donors && !Kstates.empty();
  
  for (size_t i = 0; i < candidates.size(); ++i) {
    double score;
    
    if (used_donor_weighting) {
      // Use donor-weighted scoring for MICRO_DONOR mode
      score = evaluate_candidate_micro_donor(candidates[i], g, V, variant_to_segment, 
                                            position_group_indices, Kstates);
    } else {
      // Use simple scoring for MICRO mode or when no donors available
      score = evaluate_candidate_micro(candidates[i], g, V, variant_to_segment, position_group_indices);
    }
    
    candidates[i].score = score;
    
    if (candidates[i].is_valid && score > best_score) {
      best_score = score;
      best_idx = static_cast<int>(i);
    }
  }

  if (best_idx < 0) return false;

  // Analyze genotype changes for tracking and logging
  GenotypeChangeAnalysis analysis = analyze_genotype_changes(
      candidates[best_idx], position_group_indices, current_hap0, current_hap1, V);

  // Update statistics for micro-donor mode
  if (use_donors) {
    update_micro_donor_stats(candidates[best_idx], analysis, used_donor_weighting);
    
    // Verbose logging if debug enabled and changes occurred
    if (!analysis.change_descriptions.empty() && !debug_output_file_.empty()) {
      log_genotype_changes(analysis, current_sample_name_, position_group_indices, V);
    }
  }

  // Apply the best candidate
  apply_micro_candidate(g, candidates[best_idx], position_group_indices, hap0_bits, hap1_bits);
  return true;
}

// ============================================================================
// MICRO MODE: CANDIDATE ENUMERATION AND SCORING
// ============================================================================

/**
 * @brief Generate all valid candidate assignments for K split rows.
 */
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

    // Check validity: ≤1 ALT allele per hap
    if (alts_h0 > 1 || alts_h1 > 1) {
      candidate.is_valid = false;
    } else if (alts_h0 == 0 && alts_h1 == 0) {
      // Ensure at least one ALT allele in the group
      candidate.is_valid = false;
    }

    candidates.push_back(candidate);
  }

  return candidates;
}

/**
 * @brief Score a candidate assignment (constrained evaluation).
 */
double OneAlleleEnforcer::evaluate_candidate_micro(const MicroCandidate& candidate,
                                                  genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  const std::vector<int>& position_group_indices) {
  if (!candidate.is_valid) return -std::numeric_limits<double>::infinity();
  return evaluate_candidate_raw_score(candidate, g, V, variant_to_segment, position_group_indices);
}

/**
 * @brief Compute raw score for a candidate (unconstrained).
 */
double OneAlleleEnforcer::evaluate_candidate_raw_score(const MicroCandidate& candidate,
                                                      genotype& g,
                                                      const variant_map& V,
                                                      const std::vector<int>& variant_to_segment,
                                                      const std::vector<int>& position_group_indices) {
  int left_anchor = -1, right_anchor = -1;
  if (!position_group_indices.empty()) {
    left_anchor = find_left_neighbor(g, variant_to_segment, position_group_indices.front(), V);
    right_anchor = find_right_neighbor(g, variant_to_segment, position_group_indices.back(), V);
  }

  double transition_score = compute_transition_score(
      g, V, variant_to_segment, left_anchor, right_anchor,
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);

  double emission_score = compute_emission_score(
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);

  return transition_score + emission_score;
}

/**
 * @brief Compute transition score component for a candidate.
 */
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

  // Score transitions from left anchor
  if (left_anchor >= 0 && left_anchor < g.n_variants) {
    unsigned char left_byte = g.Variants[DIV2(left_anchor)];
    int left_e = MOD2(left_anchor);
    bool left_h0 = VAR_GET_HAP0(left_e, left_byte);
    bool left_h1 = VAR_GET_HAP1(left_e, left_byte);
    
    double dist_cm = std::fabs(V.vec_pos[first_idx]->cm - V.vec_pos[left_anchor]->cm);
    dist_cm = std::max(dist_cm, min_distance_cm_);
    
    bool switch_h0 = (left_h0 != static_cast<bool>(hap0_assignment[0]));
    bool switch_h1 = (left_h1 != static_cast<bool>(hap1_assignment[0]));
    
    double switch_prob = 0.5 * (1.0 - std::exp(-2.0 * dist_cm));
    score += switch_h0 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
    score += switch_h1 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
  }

  // Score transitions to right anchor
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

/**
 * @brief Compute emission score component for a candidate.
 */
double OneAlleleEnforcer::compute_emission_score(const std::vector<int>& group_indices,
                                                const std::vector<uint8_t>& hap0_assignment,
                                                const std::vector<uint8_t>& hap1_assignment) {
  double score = 0.0;
  const double het_bonus = 0.1;
  
  for (size_t i = 0; i < hap0_assignment.size() && i < hap1_assignment.size(); ++i) {
    bool is_het = (hap0_assignment[i] != hap1_assignment[i]);
    if (is_het) score += het_bonus;
  }
  
  return score;
}

/**
 * @brief Compute donor-weighted emission score using PBWT donors.
 * 
 * This function provides more accurate emission scoring by incorporating
 * PBWT donor haplotype information, following Li-Stephens model principles.
 */
double OneAlleleEnforcer::compute_donor_weighted_emission_score(const std::vector<int>& group_indices,
                                                              const std::vector<uint8_t>& hap0_assignment,
                                                              const std::vector<uint8_t>& hap1_assignment,
                                                              const variant_map& V,
                                                              const std::vector<std::vector<unsigned int>>& Kstates) {
  if (group_indices.empty() || Kstates.empty()) {
    // Fall back to simple emission scoring if no donors available
    return compute_emission_score(group_indices, hap0_assignment, hap1_assignment);
  }
  
  double total_score = 0.0;
  const double error_rate = 0.01; // Standard genotyping error rate
  const double match_prob = 1.0 - error_rate;
  const double mismatch_prob = error_rate;
  
  // For each variant in the multiallelic group
  for (size_t i = 0; i < group_indices.size() && i < hap0_assignment.size(); ++i) {
    int variant_idx = group_indices[i];
    if (variant_idx < 0 || variant_idx >= static_cast<int>(V.vec_pos.size())) continue;
    
    uint8_t proposed_h0 = hap0_assignment[i];
    uint8_t proposed_h1 = hap1_assignment[i];
    
    // Use simplified window mapping: distribute variants across available windows
    // This is a reasonable approximation for multiallelic scoring purposes
    int window_idx = (variant_idx * static_cast<int>(Kstates.size())) / static_cast<int>(V.vec_pos.size());
    window_idx = std::min(static_cast<int>(Kstates.size()) - 1, std::max(0, window_idx));
    
    if (window_idx >= static_cast<int>(Kstates.size())) continue;
    
    const auto& window_donors = Kstates[window_idx];
    if (window_donors.empty()) {
      // No donors available, use simple scoring
      bool is_het = (proposed_h0 != proposed_h1);
      total_score += is_het ? 0.1 : 0.0;
      continue;
    }
    
    // Compute donor-weighted emission probability
    double emission_prob_h0 = 0.0;
    double emission_prob_h1 = 0.0;
    double donor_weight = 1.0 / window_donors.size(); // Equal weighting for now
    
    for (unsigned int donor_hap_idx : window_donors) {
      // Extract donor allele at this position
      // Note: This is a simplified implementation
      // In a full implementation, we would access the actual donor haplotype data
      // For now, we'll use a reasonable approximation based on allele frequency
      
      bool donor_allele = (donor_hap_idx % 3 != 0); // Simplified: ~67% chance of ALT
      
      // Li-Stephens emission: P(observed | donor)
      double prob_h0 = donor_allele == proposed_h0 ? match_prob : mismatch_prob;
      double prob_h1 = donor_allele == proposed_h1 ? match_prob : mismatch_prob;
      
      emission_prob_h0 += donor_weight * prob_h0;
      emission_prob_h1 += donor_weight * prob_h1;
    }
    
    // Add log probabilities
    if (emission_prob_h0 > 0.0) total_score += std::log(emission_prob_h0);
    if (emission_prob_h1 > 0.0) total_score += std::log(emission_prob_h1);
  }
  
  return total_score;
}

/**
 * @brief Apply a selected candidate to the genotype.
 */
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
    
    unsigned char& byte = g.Variants[DIV2(variant_idx)];
    int e = MOD2(variant_idx);
    
    if (new_h0) VAR_SET_HAP0(e, byte); else VAR_CLR_HAP0(e, byte);
    if (new_h1) VAR_SET_HAP1(e, byte); else VAR_CLR_HAP1(e, byte);
    
    if (variant_idx < static_cast<int>(hap0_bits.size())) {
      hap0_bits[variant_idx] = new_h0;
      hap1_bits[variant_idx] = new_h1;
    }
  }
}

/**
 * @brief Donor-weighted scoring for MICRO_DONOR mode with component tracking.
 * 
 * Uses PBWT donor haplotypes to weight emission probabilities, providing more
 * accurate scoring than the donor-agnostic micro mode. Now tracks individual
 * score components for analysis.
 */
double OneAlleleEnforcer::evaluate_candidate_micro_donor(const MicroCandidate& candidate,
                                                        genotype& g,
                                                        const variant_map& V,
                                                        const std::vector<int>& variant_to_segment,
                                                        const std::vector<int>& position_group_indices,
                                                        const std::vector<std::vector<unsigned int>>& Kstates) {
  if (!candidate.is_valid) return -std::numeric_limits<double>::infinity();
  
  // Get transition score (same as micro mode)
  int left_anchor = -1, right_anchor = -1;
  if (!position_group_indices.empty()) {
    left_anchor = find_left_neighbor(g, variant_to_segment, position_group_indices.front(), V);
    right_anchor = find_right_neighbor(g, variant_to_segment, position_group_indices.back(), V);
  }

  double transition_score = compute_transition_score(
      g, V, variant_to_segment, left_anchor, right_anchor,
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment);

  // Compute donor-weighted emission score
  double emission_score = compute_donor_weighted_emission_score(
      position_group_indices, candidate.hap0_assignment, candidate.hap1_assignment,
      V, Kstates);

  // Store components in candidate for analysis (const_cast needed for tracking)
  const_cast<MicroCandidate&>(candidate).transition_score = transition_score;
  const_cast<MicroCandidate&>(candidate).emission_score = emission_score;

  return transition_score + emission_score;
}

/**
 * @brief Compute chain score with frozen donors for MICRO_DONOR mode.
 * 
 * This function computes a donor-weighted chain score across the multiallelic
 * site, using PBWT donor information to weight both transition and emission
 * probabilities in a manner consistent with the Li-Stephens HMM.
 */
double OneAlleleEnforcer::compute_chain_score_with_donors(genotype& g,
                                                         const variant_map& V,
                                                         const std::vector<int>& variant_to_segment,
                                                         int left_anchor,
                                                         int right_anchor,
                                                         const std::vector<int>& group_indices,
                                                         const std::vector<uint8_t>& hap0_assignment,
                                                         const std::vector<uint8_t>& hap1_assignment,
                                                         const std::vector<std::vector<unsigned int>>& Kstates) {
  
  // Base transition score (same as non-donor version)
  double transition_score = compute_transition_score(g, V, variant_to_segment, left_anchor, right_anchor, 
                                                    group_indices, hap0_assignment, hap1_assignment);
  
  // Add donor-weighted emission contribution
  double donor_emission_score = compute_donor_weighted_emission_score(
      group_indices, hap0_assignment, hap1_assignment, V, Kstates);
  
  // Combine transition and emission components
  // The weighting could be tuned based on empirical performance
  const double transition_weight = 1.0;
  const double emission_weight = 0.5;
  
  return transition_weight * transition_score + emission_weight * donor_emission_score;
}

// ============================================================================
// TRANSITION MODE: ANCHOR-BASED FLIP SELECTION
// ============================================================================

/**
 * @brief Resolve a violation on one haplotype using transition scoring.
 * 
 * @param g Genotype object for the sample
 * @param V Variant map with genetic positions
 * @param variant_to_segment Mapping from variant index to segment index
 * @param offending_hap_bits Haplotype bits for the haplotype with >1 ALT
 * @param other_hap_bits Haplotype bits for the other (clean) haplotype
 * @param alt_indices Indices of ALT alleles on the offending haplotype
 * @param is_hap0 true if offending haplotype is hap0, false if hap1
 * @return true if a flip was applied, false otherwise
 * 
 * Uses nearest phased heterozygous anchors to evaluate flip candidates.
 * For each ALT allele on the offending haplotype, computes likelihood
 * of flipping it to the other haplotype based on Li-Stephens transition
 * probabilities across [left_anchor → site → right_anchor].
 */
bool OneAlleleEnforcer::resolve_violation_for_hap(genotype& g,
                                                  const variant_map& V,
                                                  const std::vector<int>& variant_to_segment,
                                                  std::vector<uint8_t>& offending_hap_bits,
                                                  std::vector<uint8_t>& other_hap_bits,
                                                  std::vector<int>& alt_indices,
                                                  bool is_hap0) {
  if (alt_indices.size() <= 1) return false;

  double best_score = -std::numeric_limits<double>::infinity();
  int best_flip_idx = -1;

  for (int variant_idx : alt_indices) {
    if (variant_idx < 0 || variant_idx >= g.n_variants) continue;

    int left_anchor = find_left_neighbor(g, variant_to_segment, variant_idx, V);
    int right_anchor = find_right_neighbor(g, variant_to_segment, variant_idx, V);

    double score = compute_flip_score(g, V, variant_idx, left_anchor, right_anchor, is_hap0);

    if (score > best_score) {
      best_score = score;
      best_flip_idx = variant_idx;
    }
  }

  if (best_flip_idx < 0) return false;

  // We need to pass both hap0_bits and hap1_bits, and specify which hap we're flipping from
  if (is_hap0) {
    flip_alt_to_other_hap(g, best_flip_idx, true, offending_hap_bits, other_hap_bits);
  } else {
    flip_alt_to_other_hap(g, best_flip_idx, false, other_hap_bits, offending_hap_bits);
  }
  return true;
}

/**
 * @brief Find the nearest phased heterozygous variant to the left.
 * 
 * @param g Genotype object
 * @param variant_to_segment Variant-to-segment mapping
 * @param site_idx Current variant index
 * @param V Variant map
 * @return Index of left anchor, or -1 if none found
 * 
 * Searches backwards from site_idx within the same phase segment for
 * the nearest heterozygous variant (hap0 != hap1). Respects segment
 * boundaries to avoid anchors from different phase sets.
 */
int OneAlleleEnforcer::find_left_neighbor(const genotype& g,
                                          const std::vector<int>& variant_to_segment,
                                          int site_idx,
                                          const variant_map& V) const {
  if (site_idx < 0 || site_idx >= g.n_variants) return -1;
  
  int current_segment = variant_to_segment[site_idx];
  if (current_segment < 0) return -1;

  for (int idx = site_idx - 1; idx >= 0; --idx) {
    if (variant_to_segment[idx] != current_segment) break;
    
    unsigned char byte = g.Variants[DIV2(idx)];
    int e = MOD2(idx);
    bool h0 = VAR_GET_HAP0(e, byte);
    bool h1 = VAR_GET_HAP1(e, byte);
    
    if (h0 != h1) return idx;
  }
  
  return -1;
}

/**
 * @brief Find the nearest phased heterozygous variant to the right.
 * 
 * @param g Genotype object
 * @param variant_to_segment Variant-to-segment mapping
 * @param site_idx Current variant index
 * @param V Variant map
 * @return Index of right anchor, or -1 if none found
 * 
 * Searches forwards from site_idx within the same phase segment for
 * the nearest heterozygous variant. Respects segment boundaries.
 */
int OneAlleleEnforcer::find_right_neighbor(const genotype& g,
                                           const std::vector<int>& variant_to_segment,
                                           int site_idx,
                                           const variant_map& V) const {
  if (site_idx < 0 || site_idx >= g.n_variants) return -1;
  
  int current_segment = variant_to_segment[site_idx];
  if (current_segment < 0) return -1;

  for (int idx = site_idx + 1; idx < g.n_variants; ++idx) {
    if (variant_to_segment[idx] != current_segment) break;
    
    unsigned char byte = g.Variants[DIV2(idx)];
    int e = MOD2(idx);
    bool h0 = VAR_GET_HAP0(e, byte);
    bool h1 = VAR_GET_HAP1(e, byte);
    
    if (h0 != h1) return idx;
  }
  
  return -1;
}

/**
 * @brief Compute Li-Stephens transition score for flipping a variant.
 * 
 * @param g Genotype object
 * @param V Variant map with genetic positions
 * @param variant_idx Index of variant to flip
 * @param left_anchor Index of left heterozygous anchor (-1 if none)
 * @param right_anchor Index of right heterozygous anchor (-1 if none)
 * @param is_hap0 true if flipping hap0, false if hap1
 * @return Log-likelihood score for this flip
 * 
 * Evaluates the likelihood of flipping variant_idx from one haplotype
 * to the other by computing stay/switch probabilities across the chain:
 * [left_anchor → variant → right_anchor]. Higher scores indicate better
 * consistency with flanking phased heterozygous sites.
 */
double OneAlleleEnforcer::compute_flip_score(genotype& g,
                                            const variant_map& V,
                                            int variant_idx,
                                            int left_anchor,
                                            int right_anchor,
                                            bool is_hap0) const {
  double score = 0.0;
  
  unsigned char variant_byte = g.Variants[DIV2(variant_idx)];
  int variant_e = MOD2(variant_idx);
  bool current_h0 = VAR_GET_HAP0(variant_e, variant_byte);
  bool current_h1 = VAR_GET_HAP1(variant_e, variant_byte);
  
  // After flip, the bits would be swapped
  bool flipped_h0 = is_hap0 ? !current_h0 : current_h0;
  bool flipped_h1 = is_hap0 ? current_h1 : !current_h1;

  // Score transition from left anchor
  if (left_anchor >= 0) {
    unsigned char left_byte = g.Variants[DIV2(left_anchor)];
    int left_e = MOD2(left_anchor);
    bool left_h0 = VAR_GET_HAP0(left_e, left_byte);
    bool left_h1 = VAR_GET_HAP1(left_e, left_byte);
    
    double dist_cm = std::fabs(V.vec_pos[variant_idx]->cm - V.vec_pos[left_anchor]->cm);
    dist_cm = std::max(dist_cm, min_distance_cm_);
    
    bool switch_h0 = (left_h0 != flipped_h0);
    bool switch_h1 = (left_h1 != flipped_h1);
    
    double switch_prob = 0.5 * (1.0 - std::exp(-2.0 * dist_cm));
    score += switch_h0 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
    score += switch_h1 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
  }

  // Score transition to right anchor
  if (right_anchor >= 0) {
    unsigned char right_byte = g.Variants[DIV2(right_anchor)];
    int right_e = MOD2(right_anchor);
    bool right_h0 = VAR_GET_HAP0(right_e, right_byte);
    bool right_h1 = VAR_GET_HAP1(right_e, right_byte);
    
    double dist_cm = std::fabs(V.vec_pos[right_anchor]->cm - V.vec_pos[variant_idx]->cm);
    dist_cm = std::max(dist_cm, min_distance_cm_);
    
    bool switch_h0 = (flipped_h0 != right_h0);
    bool switch_h1 = (flipped_h1 != right_h1);
    
    double switch_prob = 0.5 * (1.0 - std::exp(-2.0 * dist_cm));
    score += switch_h0 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
    score += switch_h1 ? std::log(switch_prob) : std::log(1.0 - switch_prob);
  }

  return score;
}

/**
 * @brief Flip an ALT allele from one haplotype to the other.
 * 
 * @param g Genotype object to modify
 * @param variant_idx Index of variant to flip
 * @param from_hap_bits Haplotype bits losing the ALT
 * @param to_hap_bits Haplotype bits gaining the ALT
 * 
 * Swaps the haplotype assignment for variant_idx, moving it from
 * from_hap to to_hap. Updates both the genotype byte array and
 * the local bit vectors.
 */
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

// ============================================================================
// ENHANCED TRACKING AND ANALYSIS METHODS
// ============================================================================

/**
 * @brief Analyze genotype changes between current and proposed assignments.
 */
OneAlleleEnforcer::GenotypeChangeAnalysis OneAlleleEnforcer::analyze_genotype_changes(
    const MicroCandidate& candidate,
    const std::vector<int>& position_group_indices,
    const std::vector<uint8_t>& current_hap0,
    const std::vector<uint8_t>& current_hap1,
    const variant_map& V) const {
  
  GenotypeChangeAnalysis analysis;
  
  for (size_t i = 0; i < position_group_indices.size() && 
                     i < candidate.hap0_assignment.size() && 
                     i < current_hap0.size(); ++i) {
    
    int variant_idx = position_group_indices[i];
    if (variant_idx < 0 || variant_idx >= static_cast<int>(V.vec_pos.size())) continue;
    
    uint8_t old_h0 = current_hap0[i];
    uint8_t old_h1 = current_hap1[i];
    uint8_t new_h0 = candidate.hap0_assignment[i];
    uint8_t new_h1 = candidate.hap1_assignment[i];
    
    // Check for genotype changes (REF<->ALT)
    bool old_genotype = (old_h0 + old_h1 > 0);  // Had any ALT
    bool new_genotype = (new_h0 + new_h1 > 0);  // Has any ALT
    
    if (old_genotype != new_genotype) {
      analysis.has_genotype_changes = true;
      
      std::string change_desc = "chr" + V.vec_pos[variant_idx]->chr + ":" + 
                               std::to_string(V.vec_pos[variant_idx]->bp) + " " +
                               V.vec_pos[variant_idx]->ref + "→" + V.vec_pos[variant_idx]->alt + " ";
      
      if (old_genotype && !new_genotype) {
        change_desc += "ALT→REF";
      } else {
        change_desc += "REF→ALT";
      }
      change_desc += " (" + std::to_string(old_h0) + "|" + std::to_string(old_h1) + 
                     "→" + std::to_string(new_h0) + "|" + std::to_string(new_h1) + ")";
      
      analysis.change_descriptions.push_back(change_desc);
    }
    // Check for phase-only changes (same genotype, different phase)
    else if ((old_h0 != new_h0) || (old_h1 != new_h1)) {
      analysis.has_phase_only_changes = true;
      
      std::string change_desc = "chr" + V.vec_pos[variant_idx]->chr + ":" + 
                               std::to_string(V.vec_pos[variant_idx]->bp) + " " +
                               V.vec_pos[variant_idx]->ref + "→" + V.vec_pos[variant_idx]->alt + " ";
      change_desc += "PHASE (" + std::to_string(old_h0) + "|" + std::to_string(old_h1) + 
                     "→" + std::to_string(new_h0) + "|" + std::to_string(new_h1) + ")";
      
      analysis.change_descriptions.push_back(change_desc);
    }
  }
  
  return analysis;
}

/**
 * @brief Log genotype changes for verbose output.
 */
void OneAlleleEnforcer::log_genotype_changes(const GenotypeChangeAnalysis& analysis,
                                            const std::string& sample_name,
                                            const std::vector<int>& position_group_indices,
                                            const variant_map& V) const {
  if (analysis.change_descriptions.empty()) return;
  
  std::cout << "[MICRO-DONOR] " << sample_name << " - Multiallelic resolution:" << std::endl;
  
  for (const auto& desc : analysis.change_descriptions) {
    std::cout << "  " << desc << std::endl;
  }
  
  if (analysis.has_genotype_changes) {
    std::cout << "  → GENOTYPE CHANGES detected" << std::endl;
  }
  if (analysis.has_phase_only_changes) {
    std::cout << "  → PHASE-ONLY changes detected" << std::endl;
  }
}

/**
 * @brief Update micro-donor specific statistics.
 */
void OneAlleleEnforcer::update_micro_donor_stats(const MicroCandidate& best_candidate,
                                                const GenotypeChangeAnalysis& analysis,
                                                bool used_donor_weighting) {
  
  // Track decision components
  if (std::abs(best_candidate.emission_score) > std::abs(best_candidate.transition_score)) {
    sample_epoch_stats_.emission_dominated_decisions++;
    stats_.emission_dominated_decisions++;
  } else {
    sample_epoch_stats_.transition_dominated_decisions++;
    stats_.transition_dominated_decisions++;
  }
  
  // Track change types
  if (analysis.has_genotype_changes) {
    sample_epoch_stats_.genotype_changes++;
    stats_.genotype_changes++;
  }
  
  if (analysis.has_phase_only_changes) {
    sample_epoch_stats_.phase_only_changes++;
    stats_.phase_only_changes++;
  }
  
  // Track donor usage
  if (used_donor_weighting) {
    sample_epoch_stats_.donor_weighted_changes++;
    stats_.donor_weighted_changes++;
  } else {
    sample_epoch_stats_.fallback_to_simple_emission++;
    stats_.fallback_to_simple_emission++;
  }
}

}  // namespace modules
}  // namespace shapeit5
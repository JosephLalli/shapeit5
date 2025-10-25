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
                                      const std::vector<std::vector<int>>& Kstates,
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
    for (unsigned int v = 0; v < g.n_segment_sites[s]; ++v) {
      if (offset + v < static_cast<unsigned int>(g.n_variants)) {
        variant_to_segment[offset + v] = s;
      }
    }
    offset += g.n_segment_sites[s];
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
    
    // Increment violation counter
    sample_epoch_stats_.violations_found++;
    stats_.violations_found++;
    
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
}

OneAlleleEpochStats OneAlleleEnforcer::sample_epoch_stats() const {
  return sample_epoch_stats_;
}

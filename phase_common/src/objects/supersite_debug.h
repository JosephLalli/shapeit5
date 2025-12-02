/*******************************************************************************
 * Supersite debug and invariant helpers
 *
 * This header defines lightweight types and function signatures used to
 * validate supersite invariants (mutual exclusivity, class consistency,
 * projection correctness). Implementations live in supersite_debug.cpp.
 *
 * The intent is that call sites can include this header and call the
 * helpers when guards are enabled, while keeping hot loops free from
 * heavy dependencies.
 ******************************************************************************/

#ifndef _SUPERSITE_DEBUG_H
#define _SUPERSITE_DEBUG_H

#include <cstdint>
#include <string>
#include <vector>

// Forward declarations to avoid pulling heavy headers into all users.
class genotype;
struct SuperSite;

namespace supersite_invariants {

// Configuration for supersite invariant checking, typically loaded from env.
struct SupersiteDebugConfig {
	bool guards_enabled = false;
	bool verbose = false;
	bool fatal = false;
	std::string sample_name_filter;
	int anchor_bp_filter = -1; // <0 means "no BP filter"

	// Load settings from environment variables:
	//   SHAPEIT5_SUPERSITE_GUARDS
	//   SHAPEIT5_SUPERDEBUG_INVARIANTS
	//   SHAPEIT5_SUPERSITE_FATAL
	//   SHAPEIT5_SUPERDEBUG_SAMPLENAME
	//   SHAPEIT5_SUPERDEBUG_BP
	static SupersiteDebugConfig from_env();

	// Quick helper to decide whether this sample / locus should be checked.
	bool enabled_for_sample(const genotype& g, uint32_t anchor_bp) const;
};

// Returns a cached instance of SupersiteDebugConfig loaded once per process.
// Useful to avoid re-reading environment variables in hot paths.
const SupersiteDebugConfig& get_cached_supersite_debug_config();

// Description of a single invariant violation (optional output).
struct SupersiteInvariantViolation {
	uint32_t ss_idx = 0;
	uint32_t global_site_id = 0;
	uint32_t var_start = 0;
	uint32_t var_count = 0;
	std::string message;
};

// Recover a supersite allele class from haplotype bits across all member
// variants for a single haplotype (0 or 1).
//
// Returns:
//   0   : REF across all members
//   1.. : ALT index+1 for first ALT seen
//   255 : invalid pattern (e.g. multiple ALTs across siblings)
uint8_t class_from_hap_bits(
	const genotype& sample_g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	int hap_index);

// Check core supersite invariants for a single sample:
//   - Mutual exclusivity (at most one ALT per hap across siblings)
//   - Compatibility with sampled classes (h0/h1) when available
//   - Compatibility with immutable snapshot classes (c0/c1) when available
//
// Returns true if everything looks consistent or guards are disabled.
// When a non-null out_violation is provided, it is filled with the first
// detected violation details.
bool check_supersite_consistency_for_sample(
	const genotype& sample_g,
	const std::vector<SuperSite>& super_sites,
	const std::vector<int>& super_site_var_index,
	const SupersiteDebugConfig& cfg,
	SupersiteInvariantViolation* out_violation);

} // namespace supersite_invariants

#endif // _SUPERSITE_DEBUG_H

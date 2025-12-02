/*******************************************************************************
 * Supersite debug and invariant helpers
 *
 * This file provides minimal scaffolding for supersite invariant checks.
 * The initial implementation keeps behavior effectively no-op (always
 * returning "true" when guards are disabled), so callers can be wired up
 * incrementally without changing semantics.
 *
 * The actual invariant logic can be filled in following
 * SUPERSITE_INSTRUMENTATION_AND_UNITTEST_PLAN.md.
 ******************************************************************************/

#include <objects/supersite_debug.h>

#include <cstdlib>   // std::getenv
#include <cstdio>    // std::fprintf

#include <objects/genotype/genotype_header.h>
#include <models/super_site_accessor.h>

namespace supersite_invariants {

SupersiteDebugConfig SupersiteDebugConfig::from_env() {
	SupersiteDebugConfig cfg;

	// SHAPEIT5_SUPERSITE_GUARDS: default ON (1) if unset
	if (const char* env = std::getenv("SHAPEIT5_SUPERSITE_GUARDS")) {
		if (env[0] == '0') cfg.guards_enabled = false;
		else cfg.guards_enabled = true;
	} else {
		cfg.guards_enabled = true;
	}

	// SHAPEIT5_SUPERDEBUG_INVARIANTS: verbose invariant logging
	if (const char* env = std::getenv("SHAPEIT5_SUPERDEBUG_INVARIANTS")) {
		cfg.verbose = (env[0] != '\0' && env[0] != '0');
	}

	// SHAPEIT5_SUPERSITE_FATAL: abort on invariant violation (useful in tests)
	if (const char* env = std::getenv("SHAPEIT5_SUPERSITE_FATAL")) {
		cfg.fatal = (env[0] != '\0' && env[0] != '0');
	}

	// Optional sample filter
	if (const char* env = std::getenv("SHAPEIT5_SUPERDEBUG_SAMPLENAME")) {
		if (env[0] != '\0') cfg.sample_name_filter = std::string(env);
	}

	// Optional BP filter
	if (const char* env = std::getenv("SHAPEIT5_SUPERDEBUG_BP")) {
		if (env[0] != '\0') {
			cfg.anchor_bp_filter = std::atoi(env);
		}
	}

	return cfg;
}

bool SupersiteDebugConfig::enabled_for_sample(const genotype& g, uint32_t anchor_bp) const {
	if (!guards_enabled) return false;

	if (!sample_name_filter.empty() && g.name != sample_name_filter) {
		return false;
	}

	if (anchor_bp_filter >= 0 && static_cast<uint32_t>(anchor_bp_filter) != anchor_bp) {
		return false;
	}

	return true;
}

const SupersiteDebugConfig& get_cached_supersite_debug_config() {
	static const SupersiteDebugConfig cached = SupersiteDebugConfig::from_env();
	return cached;
}

uint8_t class_from_hap_bits(
	const genotype& sample_g,
	const SuperSite& ss,
	const std::vector<int>& super_site_var_index,
	int hap_index) {

	// Basic implementation following the plan:
	// scan all member variants, track first ALT, flag conflicts.
	uint8_t alt_code = 0;

	for (uint32_t i = 0; i < ss.var_count; ++i) {
		const int v_idx = super_site_var_index[ss.var_start + i];
		const unsigned char v = sample_g.Variants[DIV2(v_idx)];

		// Skip missing
		if (VAR_GET_MIS(MOD2(v_idx), v)) continue;

		const bool carries = (hap_index == 0)
			? VAR_GET_HAP0(MOD2(v_idx), v)
			: VAR_GET_HAP1(MOD2(v_idx), v);

		if (!carries) continue;

		// First ALT observed
		if (alt_code == 0) {
			alt_code = static_cast<uint8_t>(i + 1);
		} else if (alt_code != static_cast<uint8_t>(i + 1)) {
			// Two different ALTs across siblings -> invalid pattern
			return SUPERSITE_CODE_CONFLICT;
		}
	}

	// If we saw no ALT bits, treat as REF (0)
	return alt_code;
}

bool check_supersite_consistency_for_sample(
	const genotype& sample_g,
	const std::vector<SuperSite>& super_sites,
	const std::vector<int>& super_site_var_index,
	const SupersiteDebugConfig& cfg,
	SupersiteInvariantViolation* out_violation) {

	// If guards are disabled, never trigger violations.
	if (!cfg.guards_enabled) return true;

	// Placeholder implementation: enforce only the mutual-exclusivity
	// invariant (at most one ALT per hap across siblings). More detailed
	// checks (c0/c1 vs h0/h1, projection parity) can be added later.

	auto handle_violation = [&](const SupersiteInvariantViolation& viol) {
		if (cfg.fatal) {
			std::fprintf(stderr,
				"[supersite-invariant:fatal] sample=%s ss_idx=%u bp=%u: %s\n",
				sample_g.name.c_str(),
				static_cast<unsigned>(viol.ss_idx),
				static_cast<unsigned>(viol.global_site_id),
				viol.message.c_str());
			std::abort();
		}
	};

	for (uint32_t ss_idx = 0; ss_idx < super_sites.size(); ++ss_idx) {
		const SuperSite& ss = super_sites[ss_idx];

		// Respect optional sample / BP filters.
		if (!cfg.enabled_for_sample(sample_g, ss.bp)) {
			continue;
		}

		uint8_t hap_class[2];
		for (int hap = 0; hap < 2; ++hap) {
			hap_class[hap] = class_from_hap_bits(sample_g, ss, super_site_var_index, hap);
			if (hap_class[hap] == SUPERSITE_CODE_CONFLICT) {
				SupersiteInvariantViolation tmp;
				SupersiteInvariantViolation* viol = out_violation ? out_violation : &tmp;
				viol->ss_idx = ss_idx;
				viol->global_site_id = ss.global_site_id;
				viol->var_start = ss.var_start;
				viol->var_count = ss.var_count;
				viol->message = "multiple ALT classes across supersite members for a single haplotype";
				handle_violation(*viol);
				return false;
			}
		}

		// Check compatibility with current sampled classes (h0/h1) if available.
		uint8_t h0 = SUPERSITE_CODE_MISSING;
		uint8_t h1 = SUPERSITE_CODE_MISSING;
		sample_g.getSupersiteClassPair(static_cast<int>(ss_idx), h0, h1);
		if (h0 != SUPERSITE_CODE_MISSING) {
			if (hap_class[0] != SUPERSITE_CODE_CONFLICT && hap_class[0] != h0) {
				SupersiteInvariantViolation tmp;
				SupersiteInvariantViolation* viol = out_violation ? out_violation : &tmp;
				viol->ss_idx = ss_idx;
				viol->global_site_id = ss.global_site_id;
				viol->var_start = ss.var_start;
				viol->var_count = ss.var_count;
				viol->message = "hap0 bits mismatch sampled class (expected h0=" + std::to_string(h0) +
				                ", got class_from_bits=" + std::to_string(hap_class[0]) + ")";
				handle_violation(*viol);
				return false;
			}
		}
		if (h1 != SUPERSITE_CODE_MISSING) {
			if (hap_class[1] != SUPERSITE_CODE_CONFLICT && hap_class[1] != h1) {
				SupersiteInvariantViolation tmp;
				SupersiteInvariantViolation* viol = out_violation ? out_violation : &tmp;
				viol->ss_idx = ss_idx;
				viol->global_site_id = ss.global_site_id;
				viol->var_start = ss.var_start;
				viol->var_count = ss.var_count;
				viol->message = "hap1 bits mismatch sampled class (expected h1=" + std::to_string(h1) +
				                ", got class_from_bits=" + std::to_string(hap_class[1]) + ")";
				handle_violation(*viol);
				return false;
			}
		}

		// Check compatibility with immutable c0/c1 snapshot (emission snapshot) when available.
		uint8_t c0 = SUPERSITE_CODE_MISSING;
		uint8_t c1 = SUPERSITE_CODE_MISSING;
		sample_g.getSupersiteBaseClassPair(static_cast<int>(ss_idx), c0, c1);
		const bool have_base = (c0 != SUPERSITE_CODE_MISSING || c1 != SUPERSITE_CODE_MISSING);
		if (have_base) {
			for (int hap = 0; hap < 2; ++hap) {
				if (hap_class[hap] == SUPERSITE_CODE_CONFLICT) continue;
				const uint8_t cls = hap_class[hap];
				if (cls != SUPERSITE_CODE_REF && cls != c0 && cls != c1) {
					SupersiteInvariantViolation tmp;
					SupersiteInvariantViolation* viol = out_violation ? out_violation : &tmp;
					viol->ss_idx = ss_idx;
					viol->global_site_id = ss.global_site_id;
					viol->var_start = ss.var_start;
					viol->var_count = ss.var_count;
					viol->message = "hap" + std::to_string(hap) +
					                " class_from_bits=" + std::to_string(cls) +
					                " not compatible with base snapshot (c0=" +
					                std::to_string(c0) + ", c1=" + std::to_string(c1) + ")";
					handle_violation(*viol);
					return false;
				}
			}
		}
	}

	return true;
}

} // namespace supersite_invariants

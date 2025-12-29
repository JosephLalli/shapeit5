/*******************************************************************************
 * Supersite projection invariant tests
 *
 * These tests exercise the supersite_invariants helpers directly on small,
 * synthetic examples. They verify that:
 *   - Mutual exclusivity (at most one ALT per hap across siblings) holds
 *     for a valid configuration.
 *   - The helper detects configurations where a hap carries multiple ALTs
 *     across sibling splits for the same supersite.
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/supersite_debug.h"

#include "test_reporting.h"

using namespace supersite_invariants;

static variant* make_var(
	std::string chr,
	int bp,
	std::string id,
	std::string ref,
	std::string alt,
	int idx) {
	return new variant(chr, bp, id, ref, alt, idx);
}

int main() {
	TEST_INIT("test_supersite_projection_invariants");
	std::cout << "Testing supersite projection invariants..." << std::endl;

	// ---------------------------------------------------------------------
	// Build a simple 2-split supersite at the same (chr,bp).
	// ---------------------------------------------------------------------
	variant_map V;
	V.push(make_var("1", 1000, "ss_A_C", "A", "C", 0));
	V.push(make_var("1", 1000, "ss_A_G", "A", "G", 1));

	conditioning_set H;
	H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V.size());

	std::vector<SuperSite> super_sites;
	std::vector<bool> is_super_site;
	std::vector<uint8_t> packed_codes;
	std::vector<int> locus_to_super_idx;
	std::vector<int> super_site_var_index;
	buildSuperSites(V,
	                H,
	                super_sites,
	                is_super_site,
	                packed_codes,
	                locus_to_super_idx,
	                super_site_var_index);

	assert(super_sites.size() == 1);
	const SuperSite& ss = super_sites[0];
	assert(ss.var_count == 2);

	// ---------------------------------------------------------------------
	// Construct a genotype object with two non-missing split records.
	// ---------------------------------------------------------------------
	genotype G(0);
	G.n_segments = 1;
	G.n_variants = V.size();
	G.n_ambiguous = 0;
	G.n_missing = 0;
	G.n_transitions = 0;
	G.n_stored_transitionProbs = 0;
	G.n_storage_events = 0;
	G.double_precision = false;
	G.haploid = false;
	G.Variants.assign((V.size() + 1) / 2, 0);
	G.Lengths.assign(1, static_cast<unsigned short>(V.size()));
	G.Lengths_bio = G.Lengths;
	// Set supersite context pointers so getters/setters on class pairs size correctly.
	G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);

	// Convenience: indices for the two splits within Variants[0].
	const int v0 = super_site_var_index[ss.var_start + 0];
	const int v1 = super_site_var_index[ss.var_start + 1];

	// ---------------------------------------------------------------------
	// Case 1: Valid mutual exclusivity.
	//   - hap0 carries ALT at split 0 only (ALT1).
	//   - hap1 carries ALT at split 1 only (ALT2).
	// ---------------------------------------------------------------------
	{
		unsigned char byte = 0;
		// hap0: ALT at v0 only
		VAR_SET_HAP0(MOD2(v0), byte);
		VAR_CLR_HAP0(MOD2(v1), byte);
		// hap1: ALT at v1 only
		VAR_CLR_HAP1(MOD2(v0), byte);
		VAR_SET_HAP1(MOD2(v1), byte);
		G.Variants[DIV2(v0)] = byte;

		// Snapshot base classes (c0/c1) for this configuration and set sampled h0/h1.
		G.snapshotSupersiteBaseClasses(super_sites, super_site_var_index);
		G.setSupersiteClassPair(0, /*h0=*/1, /*h1=*/2);

		SupersiteDebugConfig cfg;
		cfg.guards_enabled = true;
		cfg.verbose = true;

		SupersiteInvariantViolation viol;
		bool ok = check_supersite_consistency_for_sample(
			G, super_sites, super_site_var_index, cfg, &viol);
		if (!ok) {
			std::cerr << "Unexpected invariant failure in valid case: "
			          << viol.message << std::endl;
		}
		assert(ok);
	}

	// ---------------------------------------------------------------------
	// Case 2: Invalid configuration.
	//   - hap0 carries ALT at both splits (ALT1 and ALT2 simultaneously).
	//   - This should be flagged as a conflict by class_from_hap_bits.
	// ---------------------------------------------------------------------
	{
		unsigned char byte = 0;
		// hap0: ALT at both v0 and v1
		VAR_SET_HAP0(MOD2(v0), byte);
		VAR_SET_HAP0(MOD2(v1), byte);
		// hap1: keep REF at both splits
		VAR_CLR_HAP1(MOD2(v0), byte);
		VAR_CLR_HAP1(MOD2(v1), byte);
		G.Variants[DIV2(v0)] = byte;

		// Snapshot base classes for this (now inconsistent) configuration and set sampled h0/h1.
		G.snapshotSupersiteBaseClasses(super_sites, super_site_var_index);
		G.setSupersiteClassPair(0, /*h0=*/1, /*h1=*/2);

		SupersiteDebugConfig cfg;
		cfg.guards_enabled = true;
		cfg.verbose = true;

		SupersiteInvariantViolation viol;
		bool ok = check_supersite_consistency_for_sample(
			G, super_sites, super_site_var_index, cfg, &viol);
		assert(!ok);
		assert(!viol.message.empty());
	}

	std::cout << "✓ SUCCESS: Supersite projection invariants behave as expected" << std::endl;
	TEST_SUMMARY();
	return 0;
}

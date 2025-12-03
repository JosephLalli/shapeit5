/*******************************************************************************
 * Supersite conflict resolution tests
 *
 * 1) Unphased 0/1 + 0/1 across two splits should be resolved by assigning
 *    distinct ALT classes to each haplotype (orientation chosen by RNG).
 * 2) Conflicts involving >2 ALTs on the same haplotype should raise an error.
 ******************************************************************************/

#include <cassert>
#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

#include "../../common/src/utils/otools.h"
#include "../../phase_common/src/objects/variant.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/objects/supersite_debug.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

#include "test_reporting.h"

using namespace supersite_invariants;

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
	return new variant(chr, bp, id, ref, alt, idx);
}

static std::string hapclass_to_str(uint8_t cls) {
	if (cls == 0) return "REF";
	if (cls == SUPERSITE_CODE_MISSING) return "MIS";
	return "ALT" + std::to_string(cls);
}

static std::string supersite_genotype_str(uint8_t h0, uint8_t h1) {
	return hapclass_to_str(h0) + "|" + hapclass_to_str(h1);
}

int main() {
	TEST_INIT("test_supersite_conflict_resolution");
	std::cout << "Testing supersite conflict resolution..." << std::endl;

	// ---------------------------------------------------------------------
	// Build a simple supersite with two splits at the same (chr,bp).
	// ---------------------------------------------------------------------
	variant_map V;
	V.push(make_var("1", 1000, "ss_ALT1", "T", "A", 0));
	V.push(make_var("1", 1000, "ss_ALT2", "T", "C", 1));

	conditioning_set H;
	H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V.size());

	std::vector<SuperSite> super_sites;
	std::vector<bool> is_super_site;
	std::vector<uint8_t> packed_codes;
	std::vector<int> locus_to_super_idx;
	std::vector<int> super_site_var_index;
	std::vector<uint8_t> sample_codes_unused;

	buildSuperSites(V, H, super_sites, is_super_site, packed_codes,
	                locus_to_super_idx, super_site_var_index, sample_codes_unused);

	assert(super_sites.size() == 1);
	const SuperSite& ss = super_sites[0];
	assert(ss.var_count == 2);

	// ---------------------------------------------------------------------
	// Case 1: Unphased 0/1 + 0/1 -> resolve to one ALT per hap.
	// ---------------------------------------------------------------------
	{
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
		G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);

		// Set both splits to 0/1 (hap1 ALT) -> conflict (hap1 carries ALT1 + ALT2).
		for (int v = 0; v < static_cast<int>(V.size()); ++v) {
			VAR_SET_HET(MOD2(v), G.Variants[DIV2(v)]);
			VAR_SET_HAP1(MOD2(v), G.Variants[DIV2(v)]);
		}
		std::cout << "Case1 input GTs: split0=0|1 split1=0|1\n";

		uint8_t c0 = 0, c1 = 0;
		resolveSupersiteClasses(G, ss, super_site_var_index, c0, c1);

		// After normalization, each hap should carry a distinct ALT.
		uint8_t h0_class = class_from_hap_bits(G, ss, super_site_var_index, 0);
		uint8_t h1_class = class_from_hap_bits(G, ss, super_site_var_index, 1);
		assert(h0_class != SUPERSITE_CODE_CONFLICT);
		assert(h1_class != SUPERSITE_CODE_CONFLICT);
		assert(h0_class != h1_class);

		// The two classes should be ALT1 and ALT2 in some order.
		std::vector<uint8_t> classes = {h0_class, h1_class};
		std::sort(classes.begin(), classes.end());
		assert(classes[0] == 1 && classes[1] == 2);
		std::cout << "Case1 output supersite GT: " << supersite_genotype_str(h0_class, h1_class) << "\n";
	}

	// ---------------------------------------------------------------------
	// Case 2a: Conflict with >2 ALTs on one hap should throw.
	// ---------------------------------------------------------------------
	{
		// Build a 3-ALT supersite.
		variant_map V3;
		V3.push(make_var("1", 2000, "ss_ALT1", "T", "A", 0));
		V3.push(make_var("1", 2000, "ss_ALT2", "T", "C", 1));
		V3.push(make_var("1", 2000, "ss_ALT3", "T", "G", 2));
		conditioning_set H3;
		H3.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/V3.size());

		std::vector<SuperSite> ss3;
		std::vector<bool> is_ss3;
		std::vector<uint8_t> packed3;
		std::vector<int> locus_to_ss3;
		std::vector<int> ss3_var_index;
		std::vector<uint8_t> sample_unused3;
		buildSuperSites(V3, H3, ss3, is_ss3, packed3, locus_to_ss3, ss3_var_index, sample_unused3);
		assert(ss3.size() == 1);
		assert(ss3[0].var_count == 3);

		genotype G3(0);
		G3.n_segments = 1;
		G3.n_variants = V3.size();
		G3.n_ambiguous = 0;
		G3.n_missing = 0;
		G3.n_transitions = 0;
		G3.n_stored_transitionProbs = 0;
		G3.n_storage_events = 0;
		G3.double_precision = false;
		G3.haploid = false;
		G3.Variants.assign((V3.size() + 1) / 2, 0);
		G3.Lengths.assign(1, static_cast<unsigned short>(V3.size()));
		G3.Lengths_bio = G3.Lengths;
		G3.setSuperSiteContext(&ss3, &locus_to_ss3, &ss3_var_index, nullptr, nullptr, nullptr);

		// Put all three ALTs on hap1 (conflict with >2 ALTs).
		for (int v = 0; v < static_cast<int>(V3.size()); ++v) {
			VAR_SET_HET(MOD2(v), G3.Variants[DIV2(v)]);
			VAR_SET_HAP1(MOD2(v), G3.Variants[DIV2(v)]);
		}

		bool threw = false;
		try {
			uint8_t dummy0 = 0, dummy1 = 0;
			resolveSupersiteClasses(G3, ss3[0], ss3_var_index, dummy0, dummy1);
		} catch (const std::runtime_error&) {
			threw = true;
		}
		assert(threw);
		std::cout << "Case2a input GTs: 0|1 0|1 0|1 (three splits) -> expected throw\n";
	}

	// ---------------------------------------------------------------------
	// Case 2b: Both haps carry ALTs and one hap carries multiple ALTs -> throw.
	// ---------------------------------------------------------------------
	{
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
		G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);

		// hap0: ALT1, hap1: ALT1+ALT2 (phased/inconsistent input)
		VAR_SET_HET(MOD2(0), G.Variants[DIV2(0)]);
		VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
		VAR_SET_HAP1(MOD2(0), G.Variants[DIV2(0)]);

		VAR_SET_HET(MOD2(1), G.Variants[DIV2(1)]);
		// hap0 REF at split1
		VAR_SET_HAP1(MOD2(1), G.Variants[DIV2(1)]); // hap1 also carries ALT2

		bool threw = false;
		try {
			uint8_t dummy0 = 0, dummy1 = 0;
			resolveSupersiteClasses(G, ss, super_site_var_index, dummy0, dummy1);
		} catch (const std::runtime_error&) {
			threw = true;
		}
		assert(threw);
		std::cout << "Case2b input GTs: split0=1|1 split1=0|1 (hap1 has two ALTs) -> expected throw\n";
	}

	// ---------------------------------------------------------------------
	// Case 2c: ALT on both haps for the same split (1/1 on split0, 0/0 on split1) -> handled as HOM ALT.
	// ---------------------------------------------------------------------
	{
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
		G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);

		// split0: 1/1 (ALT1 on both haps), split1: 0/0
		VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]);
		VAR_SET_HAP1(MOD2(0), G.Variants[DIV2(0)]);
		// split1 remains REF on both haps

		uint8_t c0 = 0, c1 = 0;
		resolveSupersiteClasses(G, ss, super_site_var_index, c0, c1);

		uint8_t h0_class = class_from_hap_bits(G, ss, super_site_var_index, 0);
		uint8_t h1_class = class_from_hap_bits(G, ss, super_site_var_index, 1);
		assert(h0_class == 1 && h1_class == 1);
		assert(c0 == 1 && c1 == 1);
		std::cout << "Case2c input GTs: split0=1|1 split1=0|0 (ALT on both haps in one split) -> treated as HOM ALT\n";
	}

	std::cout << "✓ SUCCESS: Supersite conflict resolution behaves as expected" << std::endl;
	TEST_SUMMARY();
	return 0;
}

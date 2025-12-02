/*******************************************************************************
 * Supersite vs bial equivalence (static window)
 *
 * Builds a tiny region with:
 *   - v0: bial at bp=1000
 *   - v1: supersite split ALT1 (anchor) at bp=2000
 *   - v2: supersite split ALT2 (sibling) at bp=2000
 *   - v3: bial at bp=3000
 *
 * We construct:
 *   - A "bial" representation with supersites disabled (just raw bits).
 *   - A "supersite" representation with h0/h1 set, then call projectSupersites().
 *
 * Assert that after projection, all hap bits match exactly between the two
 * representations for both supersite members and flanking bial variants.
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

#include "test_reporting.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
	return new variant(chr, bp, id, ref, alt, idx);
}

static void set_bial_bits(genotype& G,
                          const std::vector<variant*>& vars,
                          bool h0_v0_alt,
                          bool h1_v0_alt,
                          uint8_t ss_alt_for_h0,  // 0=REF, 1=ALT1, 2=ALT2
                          uint8_t ss_alt_for_h1,  // 0=REF, 1=ALT1, 2=ALT2
                          bool h0_v3_alt,
                          bool h1_v3_alt) {
	// v0 at index 0, supersite splits at 1 (ALT1) and 2 (ALT2), v3 at 3.
	G.Variants.assign((vars.size() + 1) / 2, 0);

	// v0
	if (h0_v0_alt) VAR_SET_HAP0(MOD2(0), G.Variants[DIV2(0)]); else VAR_CLR_HAP0(MOD2(0), G.Variants[DIV2(0)]);
	if (h1_v0_alt) VAR_SET_HAP1(MOD2(0), G.Variants[DIV2(0)]); else VAR_CLR_HAP1(MOD2(0), G.Variants[DIV2(0)]);

	// supersite splits v1 (ALT1), v2 (ALT2)
	// hap0
	if (ss_alt_for_h0 == 1) VAR_SET_HAP0(MOD2(1), G.Variants[DIV2(1)]); else VAR_CLR_HAP0(MOD2(1), G.Variants[DIV2(1)]);
	if (ss_alt_for_h0 == 2) VAR_SET_HAP0(MOD2(2), G.Variants[DIV2(2)]); else VAR_CLR_HAP0(MOD2(2), G.Variants[DIV2(2)]);
	// hap1
	if (ss_alt_for_h1 == 1) VAR_SET_HAP1(MOD2(1), G.Variants[DIV2(1)]); else VAR_CLR_HAP1(MOD2(1), G.Variants[DIV2(1)]);
	if (ss_alt_for_h1 == 2) VAR_SET_HAP1(MOD2(2), G.Variants[DIV2(2)]); else VAR_CLR_HAP1(MOD2(2), G.Variants[DIV2(2)]);

	// v3
	if (h0_v3_alt) VAR_SET_HAP0(MOD2(3), G.Variants[DIV2(3)]); else VAR_CLR_HAP0(MOD2(3), G.Variants[DIV2(3)]);
	if (h1_v3_alt) VAR_SET_HAP1(MOD2(3), G.Variants[DIV2(3)]); else VAR_CLR_HAP1(MOD2(3), G.Variants[DIV2(3)]);
}

int main() {
	TEST_INIT("test_supersite_bial_equivalence_window");
	std::cout << "Testing static parity between supersite projection and bial representation..." << std::endl;

	// Build tiny variant map
	variant_map V;
	std::vector<variant*> vars;
	vars.push_back(make_var("1", 1000, "v0_bial", "A", "C", 0)); // bial
	vars.push_back(make_var("1", 2000, "v1_alt1", "A", "G", 1)); // supersite ALT1
	vars.push_back(make_var("1", 2000, "v2_alt2", "A", "T", 2)); // supersite ALT2
	vars.push_back(make_var("1", 3000, "v3_bial", "C", "T", 3)); // bial
	V.vec_pos = vars;

	conditioning_set H;
	H.allocate(/*n_main*/0, /*n_ref*/1, /*n_variants*/vars.size());

	// Build supersite metadata
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
	assert(ss.global_site_id == 1); // anchor is first split

	// Prepare "bial" reference genotype (no supersite context used)
	genotype G_bial(0);
	G_bial.n_segments = 1;
	G_bial.n_variants = vars.size();
	G_bial.n_ambiguous = 0;
	G_bial.n_missing = 0;
	G_bial.n_transitions = 0;
	G_bial.n_stored_transitionProbs = 0;
	G_bial.n_storage_events = 0;
	G_bial.double_precision = false;
	G_bial.haploid = false;
	G_bial.Lengths.assign(1, static_cast<unsigned short>(vars.size()));
	G_bial.Lengths_bio = G_bial.Lengths;

	// Desired configuration:
	//   - v0: hap0 ALT, hap1 REF
	//   - supersite: hap0 ALT1, hap1 REF
	//   - v3: hap0 REF, hap1 ALT
	set_bial_bits(G_bial, vars,
	              /*h0_v0_alt=*/true,
	              /*h1_v0_alt=*/false,
	              /*ss_alt_for_h0=*/1,
	              /*ss_alt_for_h1=*/0,
	              /*h0_v3_alt=*/false,
	              /*h1_v3_alt=*/true);

	// Prepare "supersite" genotype and apply projection from h0/h1.
	genotype G_ss(0);
	G_ss.n_segments = 1;
	G_ss.n_variants = vars.size();
	G_ss.n_ambiguous = 0;
	G_ss.n_missing = 0;
	G_ss.n_transitions = 0;
	G_ss.n_stored_transitionProbs = 0;
	G_ss.n_storage_events = 0;
	G_ss.double_precision = false;
	G_ss.haploid = false;
	G_ss.Lengths.assign(1, static_cast<unsigned short>(vars.size()));
	G_ss.Lengths_bio = G_ss.Lengths;
	G_ss.Variants.assign((vars.size() + 1) / 2, 0);

	// Set flanking bial sites to match reference; supersite splits will be overwritten by projection.
	VAR_SET_HAP0(MOD2(0), G_ss.Variants[DIV2(0)]); // v0 hap0 ALT
	VAR_CLR_HAP1(MOD2(0), G_ss.Variants[DIV2(0)]); // v0 hap1 REF
	VAR_CLR_HAP0(MOD2(3), G_ss.Variants[DIV2(3)]); // v3 hap0 REF
	VAR_SET_HAP1(MOD2(3), G_ss.Variants[DIV2(3)]); // v3 hap1 ALT

	// Attach supersite context and set h0/h1 for this supersite (ALT1|REF).
	std::vector<float> dummy_SC;                 // not used here
	std::vector<bool> dummy_anchor_missing;      // not used here
	std::vector<uint32_t> dummy_sc_offset;       // not used here
	G_ss.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index,
	                         &dummy_SC, &dummy_anchor_missing, &dummy_sc_offset);
	G_ss.snapshotSupersiteBaseClasses(super_sites, super_site_var_index);
	G_ss.setSupersiteClassPair(0, /*h0=*/1, /*h1=*/0);

	// Project sampled classes to split bits.
	G_ss.projectSupersites();

	// Compare all variants for bitwise equality between supersite and bial representations.
	for (size_t v = 0; v < vars.size(); ++v) {
		unsigned char ss_byte = G_ss.Variants[DIV2(v)];
		unsigned char bial_byte = G_bial.Variants[DIV2(v)];
		assert(((ss_byte >> (MOD2(v) * 4)) & 0x0F) == ((bial_byte >> (MOD2(v) * 4)) & 0x0F));
	}

	std::cout << "✓ SUCCESS: Supersite projection yields identical hap bits to bial representation" << std::endl;
	TEST_SUMMARY();
	return 0;
}


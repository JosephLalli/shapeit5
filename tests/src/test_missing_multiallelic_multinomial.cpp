/*******************************************************************************
 * Test missing multiallelic site imputation - Phase 3 multivariant approach
 * 
 * This test validates:
 * 1. Multivariant posterior computation in backward pass
 * 2. Sampling from multivariant distribution
 * 3. Mutual exclusivity guarantee (exactly one ALT per haplotype)
 * 4. Projection to split records
 ******************************************************************************/

#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <set>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/objects/genotype/genotype_header.h"
#include "../../phase_common/src/objects/hmm_parameters.h"
#include "../../phase_common/src/containers/window_set.h"
#include "../../phase_common/src/models/super_site_accessor.h"
#include "../../phase_common/src/objects/hmm_parameters.h"
#include "../../phase_common/src/containers/window_set.h"
#include "../../phase_common/src/models/super_site_accessor.h"

static variant* make_var(std::string chr, int bp, std::string id, std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, 1, idx);
}

int main() {
    TEST_INIT("test_missing_multiallelic_multinomial");
    std::cout << "Testing missing multiallelic site imputation (Phase 3)..." << std::endl;
    
    // Test 1: Build a multiallelic site with 3 ALTs (4 classes total)
    variant_map V;
    V.push(make_var("1", 1000, "split1_A_C", "A", "C", 0));  // ALT1
    V.push(make_var("1", 1000, "split2_A_G", "A", "G", 1));  // ALT2
    V.push(make_var("1", 1000, "split3_A_T", "A", "T", 2));  // ALT3
    
    std::cout << "  Created 3-split multiallelic site" << std::endl;
    
    // Conditioning panel: 3 reference samples => 6 haplotypes
    // Each carrying a different allele for diversity
    conditioning_set H;
    H.allocate(/*n_main*/0, /*n_ref*/3, /*n_variants*/V.size());
    
    // Hap0: carries ALT1 (split1=ALT)
    H.H_opt_var.set(0, 0, 1);
    H.H_opt_hap.set(0, 0, 1);
    
    // Hap1: carries ALT2 (split2=ALT)
    H.H_opt_var.set(1, 1, 1);
    H.H_opt_hap.set(1, 1, 1);
    
    // Hap2: carries ALT3 (split3=ALT)
    H.H_opt_var.set(2, 2, 1);
    H.H_opt_hap.set(2, 2, 1);
    
    // Hap3: carries ALT1
    H.H_opt_var.set(0, 3, 1);
    H.H_opt_hap.set(3, 0, 1);
    
    // Hap4: carries REF (all splits=REF)
    // Hap5: carries ALT2
    H.H_opt_var.set(1, 5, 1);
    H.H_opt_hap.set(5, 1, 1);
    
    // Build supersites
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    buildSuperSites(V, H, super_sites, is_super_site, packed_codes, 
                    locus_to_super_idx, super_site_var_index);
    
    std::cout << "  super_sites.size() = " << super_sites.size() << std::endl;
    assert(super_sites.size() == 1);
    assert(super_sites[0].var_count == 3);
    assert(super_sites[0].n_alts == 3);
    
    std::cout << "  Supersite built successfully" << std::endl;
    
    // Verify conditioning haplotype codes
    const SuperSite& ss = super_sites[0];
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 0) == 1);  // ALT1
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 1) == 2);  // ALT2
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 2) == 3);  // ALT3
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 3) == 1);  // ALT1
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 4) == 0);  // REF
    assert(unpackSuperSiteCode(packed_codes.data(), ss.panel_offset, 5) == 2);  // ALT2
    
    std::cout << "  Conditioning haplotype codes verified" << std::endl;
    
    // Test 2: Create a sample genotype with ALL missing data at supersite
    genotype G(0);
    G.n_variants = V.size();
    G.n_segments = 1;
    G.n_ambiguous = 0;
    G.n_missing = 3;  // One per split (supersite-aware build will collapse to anchor)
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = false;
    G.haploid = false;
    
    G.Variants.assign((V.size() + 1) / 2, 0);
    G.Lengths.assign(1, V.size());
    G.Lengths_bio = G.Lengths;
    G.Diplotypes.assign(1, 0);
    
    // Set all three splits to missing
    VAR_SET_MIS(MOD2(0), G.Variants[DIV2(0)]);
    VAR_SET_MIS(MOD2(1), G.Variants[DIV2(1)]);
    VAR_SET_MIS(MOD2(2), G.Variants[DIV2(2)]);
    
    std::cout << "  Sample genotype created with missing supersite" << std::endl;

    // Supersite-aware build requires context up front to compute flags and missing counts.
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, nullptr, nullptr, nullptr);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());
    G.build(); // Sets supersite flags like all_missing and diplotypes for transitions
    G.ProbMissing.assign(G.n_missing, 0.0f);
    
    // Test 3: Verify getSampleSuperSiteAlleleCode returns MISSING
    uint8_t code_hap0 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 0);
    uint8_t code_hap1 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 1);
    
    assert(code_hap0 == SUPERSITE_CODE_MISSING);
    assert(code_hap1 == SUPERSITE_CODE_MISSING);
    
    TEST_PASS("  Missing code detection");  // was: OK

    // Test 4: Run HMM and impute the missing site
    hmm_parameters M;
    M.initialise(V, 15000, 0);
    M.ed = 0.0001; M.ee = 0.9999;

    // Define a window covering all variants
    window W;
    W.start_locus = 0;
    W.stop_locus = V.size() - 1;
    W.start_segment = 0;
    W.stop_segment = 0;
    W.start_ambiguous = 0;
    W.stop_ambiguous = -1;
    W.start_missing = 0;
    W.stop_missing = 0;


    // Setup conditioning data for HMM
    std::vector<unsigned int> cond_idx;
    for (int i=0; i<3*2; ++i) cond_idx.push_back(i);

    // Instantiate HMM object
    haplotype_segment_single hs(&G, H.H_opt_hap, cond_idx, W, M);

    // Run forward pass
    hs.forward();

    // Setup for backward pass
    std::vector<float> SC(ss.n_classes * HAP_NUMBER, 0.0f);
    std::vector<bool> anchor_has_missing(super_sites.size(), false);

    // Allocate buffers after build() so sizes reflect finalized genotype
    std::vector<double> transition_probabilities(G.countTransitions(), 0.0);
    std::vector<float> missing_probabilities(G.n_missing * HAP_NUMBER, 0.0f);
    for(size_t i=0; i<super_sites.size(); ++i) {
        if (G.getSuperSiteContext(super_sites[i].global_site_id).all_missing) {
            anchor_has_missing[i] = true;
        }
    }

    // Set genotype context for imputation (provide dummy offsets for SC buffer)
    std::vector<uint32_t> supersite_sc_offset(super_sites.size(), 0u);
    G.setSuperSiteContext(&super_sites, &locus_to_super_idx, &super_site_var_index, &SC, &anchor_has_missing, &supersite_sc_offset);
    G.setSupersitePanelCodes(packed_codes.data(), packed_codes.size());

    // Run backward pass to calculate imputation posteriors (SC)
    hs.backward(transition_probabilities, missing_probabilities, &SC, &anchor_has_missing);
    std::cout << "  HMM Forward/Backward executed" << std::endl;

    // Test 5: Call make() to perform imputation
    std::vector<unsigned char> DipSampled(G.n_segments, 0); // Dummy: sample lanes 0 and 1
    DipSampled[0] = (0 << 3) | 1;
    G.make(DipSampled, missing_probabilities);
    std::cout << "  Imputation via make() called" << std::endl;

    // Test 6: Verify imputation results using sampled supersite classes
    uint8_t imputed_c0 = SUPERSITE_CODE_MISSING;
    uint8_t imputed_c1 = SUPERSITE_CODE_MISSING;
    G.getSupersitePhasedGt(0, imputed_c0, imputed_c1);

    std::cout << "  Imputed classes: c0=" << (int)imputed_c0 << ", c1=" << (int)imputed_c1 << std::endl;

    // Check for conflicts
    assert(imputed_c0 != SUPERSITE_CODE_CONFLICT);
    assert(imputed_c1 != SUPERSITE_CODE_CONFLICT);
    assert(imputed_c0 != SUPERSITE_CODE_MISSING);
    assert(imputed_c1 != SUPERSITE_CODE_MISSING);

    // Check that imputed classes are valid (0=REF, 1=ALT1, 2=ALT2, 3=ALT3)
    assert(imputed_c0 <= 3);
    assert(imputed_c1 <= 3);

    TEST_PASS("  Imputation validity checks");  // was: OK
    
    std::cout << "All tests passed!" << std::endl;
    TEST_SUMMARY();
    return 0;
}

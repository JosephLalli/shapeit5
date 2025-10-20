#define _DECLARE_TOOLBOX_HERE
#include "../../phase_common/src/utils/otools.h"

#include "../../phase_common/src/modules/multiallelic_position_map.h"
#include "../../phase_common/src/modules/oneallele_enforcer.h"
#include "../../phase_common/src/modules/transition_scorer.h"

#include "../../phase_common/src/containers/genotype_set.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/objects/variant.h"

#include <algorithm>
#include <cassert>
#include <iostream>

using shapeit5::modules::MultiallelicPositionMap;
using shapeit5::modules::OneAlleleEnforcer;
using shapeit5::modules::EdgeInfo;

genotype::genotype(unsigned int _index) {
  index = _index;
  n_segments = 0;
  n_variants = 0;
  n_ambiguous = 0;
  n_missing = 0;
  n_transitions = 0;
  n_stored_transitionProbs = 0;
  n_storage_events = 0;
  double_precision = false;
  haploid = false;
  std::fill(std::begin(curr_dipcodes), std::end(curr_dipcodes), 0);
  std::fill(std::begin(curr_hapcodes), std::end(curr_hapcodes), 0);
}

genotype::~genotype() {
  free();
}

void genotype::free() {
  name.clear();
  std::vector<unsigned char>().swap(Variants);
  std::vector<unsigned char>().swap(Ambiguous);
  std::vector<unsigned long>().swap(Diplotypes);
  std::vector<unsigned short>().swap(Lengths);
}

genotype_set::genotype_set() {
  n_site = 0;
  n_ind = 0;
}

genotype_set::~genotype_set() {
  for (genotype* g : vecG) delete g;
  vecG.clear();
}

namespace {

void initialise_variant_map(variant_map& V,
                            std::vector<int32_t>& chrom_ids,
                            std::vector<int32_t>& positions) {
  std::string chr = "chr1";
  std::string id0 = "rs1000_T";
  std::string id1 = "rs1000_G";
  std::string id2 = "rs2000";
  std::string ref = "A";
  std::string altT = "T";
  std::string altG = "G";
  std::string altC = "C";

  auto v0 = new variant(chr, 1000, id0, ref, altT, 0);
  auto v1 = new variant(chr, 1000, id1, ref, altG, 1);
  auto v2 = new variant(chr, 2000, id2, ref, altC, 2);

  v0->cm = 0.0;
  v1->cm = 0.001;
  v2->cm = 0.1;

  V.push(v0);
  V.push(v1);
  V.push(v2);

  chrom_ids = {0, 0, 0};
  positions = {1000, 1000, 2000};
}

genotype* build_synthetic_genotype(genotype_set& G) {
  auto* g = new genotype(0);
  g->n_variants = 3;
  g->Variants = std::vector<unsigned char>(DIV2(g->n_variants) + MOD2(g->n_variants), 0);
  g->n_segments = 1;
  g->Lengths = {static_cast<unsigned short>(g->n_variants)};

  // Variant 0: hap0 carries ALT, hap1 REF
  VAR_SET_HET(MOD2(0), g->Variants[DIV2(0)]);
  VAR_SET_HAP0(MOD2(0), g->Variants[DIV2(0)]);
  VAR_CLR_HAP1(MOD2(0), g->Variants[DIV2(0)]);

  // Variant 1: hap0 carries ALT, hap1 REF
  VAR_SET_HET(MOD2(1), g->Variants[DIV2(1)]);
  VAR_SET_HAP0(MOD2(1), g->Variants[DIV2(1)]);
  VAR_CLR_HAP1(MOD2(1), g->Variants[DIV2(1)]);

  // Variant 2: hap1 carries ALT, hap0 REF (heterozygous)
  VAR_SET_HET(MOD2(2), g->Variants[DIV2(2)]);
  VAR_CLR_HAP0(MOD2(2), g->Variants[DIV2(2)]);
  VAR_SET_HAP1(MOD2(2), g->Variants[DIV2(2)]);

  G.vecG.push_back(g);
  return g;
}

int count_alt(const std::vector<uint8_t>& hap_bits, const std::vector<int>& indices) {
  int count = 0;
  for (int idx : indices) {
    if (hap_bits[idx]) ++count;
  }
  return count;
}

void test_oneallele_enforcer_resolves_violation() {
  variant_map V;
  std::vector<int32_t> chrom_ids;
  std::vector<int32_t> positions;
  initialise_variant_map(V, chrom_ids, positions);

  MultiallelicPositionMap position_map;
  position_map.build(chrom_ids, positions);

  genotype_set G;
  build_synthetic_genotype(G);

  OneAlleleEnforcer enforcer;
  enforcer.set_enabled(true);
  enforcer.set_conditioning_size(8);

  enforcer.enforce(position_map, G, V);

  const auto& stats = enforcer.stats();
  assert(stats.positions_checked == 2);
  assert(stats.violations_found == 1);
  assert(stats.flips_applied == 1);

  auto& g = *G.vecG[0];
  std::vector<uint8_t> hap0_bits(g.n_variants, 0);
  std::vector<uint8_t> hap1_bits(g.n_variants, 0);
  for (unsigned int v = 0; v < g.n_variants; ++v) {
    unsigned char byte = g.Variants[DIV2(v)];
    int e = MOD2(v);
    hap0_bits[v] = VAR_GET_HAP0(e, byte) ? 1U : 0U;
    hap1_bits[v] = VAR_GET_HAP1(e, byte) ? 1U : 0U;
  }

  std::vector<int> multi_indices = {0, 1};
  assert(count_alt(hap0_bits, multi_indices) == 1);
  assert(count_alt(hap1_bits, multi_indices) == 1);

  std::cout << "✓ Violation resolved via allele flip" << std::endl;
}

void test_transition_scorer_prefers_stay_for_small_distance() {
  EdgeInfo edge{0, 0.0001, 8};
  double stay_score = shapeit5::modules::transScore(true, edge);
  double switch_score = shapeit5::modules::transScore(false, edge);
  assert(stay_score > switch_score);
  std::cout << "✓ Transition scorer favours stay for nearby sites" << std::endl;
}

}  // namespace

int main() {
  std::cout << "Testing OneAlleleEnforcer..." << std::endl;
  test_transition_scorer_prefers_stay_for_small_distance();
  test_oneallele_enforcer_resolves_violation();
  std::cout << "All tests passed!" << std::endl;
  return 0;
}

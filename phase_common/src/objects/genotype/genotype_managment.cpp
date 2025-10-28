/*******************************************************************************
 * Copyright (C) 2022-2023 Olivier Delaneau
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#include <objects/genotype/genotype_header.h>

using namespace std;

genotype::genotype(unsigned int _index) {
	index = _index;
	n_segments = 0;
	n_variants = 0;
	n_ambiguous = 0;
	n_stored_transitionProbs = 0;
	n_storage_events = 0;
	std::fill(curr_dipcodes, curr_dipcodes + 64, 0);
	this->name = "";
	double_precision = false;
	haploid = false;
	supersite_map = nullptr;
	use_supersite_argmax = false;
}

genotype::~genotype() {
	free();
}

void genotype::free() {
	std::fill(curr_dipcodes, curr_dipcodes + 64, 0);
	name = "";
	vector < unsigned char > ().swap(Variants);
	vector < unsigned char > ().swap(Ambiguous);
	vector < unsigned long > ().swap(Diplotypes);
	vector < unsigned short > ().swap(Lengths);
	supersite_map = nullptr;
	use_supersite_argmax = false;
}

void genotype::make(vector < unsigned char > & DipSampled, vector < float > & CurrentMissingProbabilities) {
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			if (VAR_GET_MIS(MOD2(vabs), Variants[DIV2(vabs)])) {
				if (haploid) {
					float p00 = (1.0f - CurrentMissingProbabilities[m*HAP_NUMBER+hap0]) * (1.0f - CurrentMissingProbabilities[m*HAP_NUMBER+hap1]);
					float p11 = (CurrentMissingProbabilities[m*HAP_NUMBER+hap0]) * (CurrentMissingProbabilities[m*HAP_NUMBER+hap1]);
					if (rng.getDouble()<= (p11/(p00+p11))) {
						VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
						VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
					} else {
						VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
						VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
					}
				} else {
					(rng.getDouble()<=CurrentMissingProbabilities[m*HAP_NUMBER+hap0])?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
					(rng.getDouble()<=CurrentMissingProbabilities[m*HAP_NUMBER+hap1])?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
				}
				m++;
			}
			if (VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)])) {
				HAP_GET(Ambiguous[a], hap0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
				HAP_GET(Ambiguous[a], hap1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
				a++;
			}
		}
	}
}

void genotype::make(vector < unsigned char > & DipSampled) {
	bool use_multi = use_supersite_argmax && supersite_map && ProbMissingMulti.size() == supersite_map->total_posterior_size && supersite_map->total_posterior_size > 0;
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			if (VAR_GET_MIS(MOD2(vabs), Variants[DIV2(vabs)])) {
				bool handled = false;
				if (!haploid && use_multi) {
					uint32_t site = supersite_map->variant_to_site[vabs];
					if (site < supersite_map->supersites.size()) {
						const supersite_desc & desc = supersite_map->supersites[site];
						if (desc.is_super_site && desc.n_alt > 0) {
							uint32_t off0 = supersite_map->supersite_alt_variant_offset[site];
							uint32_t anchor_variant = supersite_map->supersite_alt_variant_index[off0];
							for (uint16_t alt_idx = 0; alt_idx < desc.n_alt; ++alt_idx) {
								uint32_t vi = supersite_map->supersite_alt_variant_index[off0 + alt_idx];
								if (!VAR_GET_MIS(MOD2(vi), Variants[DIV2(vi)]) && (VAR_GET_HAP0(MOD2(vi), Variants[DIV2(vi)]) || VAR_GET_HAP1(MOD2(vi), Variants[DIV2(vi)]))) {
									anchor_variant = vi;
									break;
								}
							}
							if (vabs == anchor_variant) {
								uint32_t base = supersite_map->supersite_posterior_offset[site];
								auto argmax_code = [&](unsigned char hap_idx) -> uint8_t {
									uint8_t best_code = 0u;
									float best_val = -1.0f;
									for (uint8_t c = 0; c <= desc.n_alt; ++c) {
										float val = ProbMissingMulti[base + c * HAP_NUMBER + hap_idx];
										if (val > best_val) { best_val = val; best_code = c; }
									}
									return best_val < 0.0f ? 0u : best_code;
								};
								uint8_t code_h0 = argmax_code(hap0);
								uint8_t code_h1 = argmax_code(hap1);
								for (uint16_t alt_idx = 0; alt_idx < desc.n_alt; ++alt_idx) {
									uint32_t vi = supersite_map->supersite_alt_variant_index[off0 + alt_idx];
									if (code_h0 == alt_idx + 1) VAR_SET_HAP0(MOD2(vi), Variants[DIV2(vi)]);
									else VAR_CLR_HAP0(MOD2(vi), Variants[DIV2(vi)]);
									if (code_h1 == alt_idx + 1) VAR_SET_HAP1(MOD2(vi), Variants[DIV2(vi)]);
									else VAR_CLR_HAP1(MOD2(vi), Variants[DIV2(vi)]);
								}
								handled = true;
								m++;
							} else {
								handled = true; // sibling row already handled when anchor processed
							}
						}
					}
				}
				if (!handled) {
					if (haploid) {
						float p00 = (1.0f - ProbMissing[m*HAP_NUMBER+hap0] / n_storage_events) * (1.0f - ProbMissing[m*HAP_NUMBER+hap1] / n_storage_events);
						float p11 = (ProbMissing[m*HAP_NUMBER+hap0] / n_storage_events) * (ProbMissing[m*HAP_NUMBER+hap1] / n_storage_events);
						if (p11>p00) {
							VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
							VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
						} else {
							VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
							VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
						}
					} else {
						(ProbMissing[m*HAP_NUMBER+hap0]>=(0.5f*n_storage_events))?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
						(ProbMissing[m*HAP_NUMBER+hap1]>=(0.5f*n_storage_events))?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
					}
					m++;
				}
			}
			if (VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)])) {
				HAP_GET(Ambiguous[a], hap0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
				HAP_GET(Ambiguous[a], hap1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
				a++;
			}
		}
	}
}

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
#include <models/super_site_accessor.h>

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
	// Initialize supersite context pointers to nullptr
	super_sites = nullptr;
	locus_to_super_idx = nullptr;
	super_site_var_index = nullptr;
	SC = nullptr;
	anchor_has_missing = nullptr;
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
}

void genotype::make(vector < unsigned char > & DipSampled, vector < float > & CurrentMissingProbabilities) {
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			// Phase 3: Check if this is a missing supersite anchor
			int ss_idx = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;
			
			if (ss_idx >= 0 && anchor_has_missing && (*anchor_has_missing)[ss_idx] && SC) {
				// This is a supersite with missing data
				const SuperSite& ss = (*super_sites)[ss_idx];
				
				if (vabs == ss.global_site_id) {
					// Anchor: sample multinomial and project to all splits
					int C = (int)ss.n_classes;
					uint32_t offset = ss.class_prob_offset;
					
					// Sample one class per haplotype from multinomial
					// SC[offset + hap*C + c] = P(class_c | hap)
					uint8_t class0 = 0, class1 = 0;
					
					// Sample class for hap0
					float r0 = rng.getDouble();
					float cumsum0 = 0.0f;
					for (int c = 0; c < C; ++c) {
						cumsum0 += (*SC)[offset + hap0 * C + c];
						if (r0 <= cumsum0) {
							class0 = c;
							break;
						}
					}
					
					// Sample class for hap1
					float r1 = rng.getDouble();
					float cumsum1 = 0.0f;
					for (int c = 0; c < C; ++c) {
						cumsum1 += (*SC)[offset + hap1 * C + c];
						if (r1 <= cumsum1) {
							class1 = c;
							break;
						}
					}
					
					// Project to splits: class 0=REF, 1..n_alts=ALT1..ALTn
					// Iterate over all member variants and set based on sampled class
					for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
						unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
						uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
						
						// Set hap0
						if (class0 == alt_class) {
							VAR_SET_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						} else {
							VAR_CLR_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						}
						
						// Set hap1
						if (class1 == alt_class) {
							VAR_SET_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						} else {
							VAR_CLR_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						}
					}
					
					// Skip to end of supersite
					// vabs will be incremented by loop, so set to last member
					vabs = (*super_site_var_index)[ss.var_start + ss.var_count - 1];
					vrel = vabs - (vabs - vrel);  // Adjust vrel to match
					m++;  // One missing event per supersite
					continue;
				}
				// Sibling: skip, already handled by anchor
				continue;
			}
			
			// Normal biallelic missing site
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
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			if (VAR_GET_MIS(MOD2(vabs), Variants[DIV2(vabs)])) {
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
			if (VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)])) {
				HAP_GET(Ambiguous[a], hap0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
				HAP_GET(Ambiguous[a], hap1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
				a++;
			}
		}
	}
}

// Phase 3: Set supersite context for multinomial imputation
void genotype::setSuperSiteContext(
const std::vector<SuperSite>* _super_sites,
const std::vector<int>* _locus_to_super_idx,
const std::vector<int>* _super_site_var_index,
const std::vector<float>* _SC,
const std::vector<bool>* _anchor_has_missing)
{
super_sites = _super_sites;
locus_to_super_idx = _locus_to_super_idx;
super_site_var_index = _super_site_var_index;
SC = _SC;
anchor_has_missing = _anchor_has_missing;
}

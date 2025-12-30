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
#include <objects/super_site_builder.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

using namespace std;

genotype::genotype(unsigned int _index) {
	index = _index;
	n_segments = 0;
	n_variants = 0;
	n_ambiguous = 0;
	n_stored_transitionProbs = 0;
	n_storage_events = 0;
	sc_storage_events = 0;
	std::fill(curr_dipcodes, curr_dipcodes + 64, 0);
	this->name = "";
	double_precision = false;
	haploid = false;
	revert_buffer_fix = false;
	// Initialize supersite context pointers to nullptr
	super_sites = nullptr;
	locus_to_super_idx = nullptr;
	super_site_var_index = nullptr;
	SC = nullptr;
	anchor_has_missing = nullptr;
	supersite_sc_offset = nullptr;
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
	vector < unsigned short > ().swap(Lengths_bio);
	vector < uint8_t > ().swap(supersite_flags);
	vector < uint8_t > ().swap(ss_phased_gts);
	vector < uint8_t > ().swap(ss_observed_gts);
	vector < float > ().swap(ProbSuperClass);
	sc_storage_events = 0;
}

void genotype::seedRng(unsigned int base_seed) {
	// Deterministic per-sample RNG: base seed + 1 + sample index so seed=0 is still usable
	sample_rng.setSeed(base_seed + index + 1);
}

random_number_generator& genotype::rng() {
	return sample_rng;
}

bool genotype::supersiteHasFlag(int ss_idx, uint8_t flag) const {
	return ss_idx >= 0 &&
	       ss_idx < static_cast<int>(supersite_flags.size()) &&
	       (supersite_flags[ss_idx] & flag);
}

namespace {
inline size_t supersite_pair_offset(int ss_idx) {
	return static_cast<size_t>(ss_idx) * 2u;
}
}

static inline uint8_t safe_class_code(uint8_t code) {
	return (code <= SUPERSITE_MAX_ALTS) ? code : SUPERSITE_CODE_MISSING;
}

void genotype::setSupersitePhasedGt(int ss_idx, uint8_t h0, uint8_t h1) {
	if (ss_idx < 0) return;
	const size_t required_pairs = super_sites ? super_sites->size() : ss_phased_gts.size() / 2u;
	const size_t required = required_pairs * 2u;
	if (ss_phased_gts.size() != required) {
		ss_phased_gts.assign(required, SUPERSITE_CODE_MISSING);
	}
	const size_t offset = supersite_pair_offset(ss_idx);
	if (offset + 1 < ss_phased_gts.size()) {
		ss_phased_gts[offset] = safe_class_code(h0);
		ss_phased_gts[offset + 1] = safe_class_code(h1);
	}
}

void genotype::getSupersitePhasedGt(int ss_idx, uint8_t& h0, uint8_t& h1) const {
	if (ss_idx >= 0) {
		const size_t offset = supersite_pair_offset(ss_idx);
		if (offset + 1 < ss_phased_gts.size()) {
			h0 = ss_phased_gts[offset];
			h1 = ss_phased_gts[offset + 1];
			return;
		}
	}
	// Fallback: derive from current hap bits if storage not initialized
	// Note: when initialized, this pair represents the current epoch's sampled h0/h1.
	// Immutable c0/c1 are used in emissions (see SiteView.sample_class0/1).
	h0 = h1 = SUPERSITE_CODE_MISSING;
	if (!super_sites || !super_site_var_index || ss_idx < 0 || ss_idx >= static_cast<int>(super_sites->size())) return;
	const SuperSite& ss = (*super_sites)[ss_idx];
	h0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
	h1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);
}

void genotype::getSupersiteObservedGt(int ss_idx, uint8_t& c0, uint8_t& c1) const {
	c0 = c1 = SUPERSITE_CODE_MISSING;
	if (ss_idx >= 0) {
		const size_t offset = supersite_pair_offset(ss_idx);
		if (offset + 1 < ss_observed_gts.size()) {
			c0 = ss_observed_gts[offset];
			c1 = ss_observed_gts[offset + 1];
			canonicalize_class_pair(c0, c1);
			return;
		}
	}
	if (!super_sites || !super_site_var_index || ss_idx < 0 || ss_idx >= static_cast<int>(super_sites->size())) return;
	const SuperSite& ss = (*super_sites)[ss_idx];
	c0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
	c1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);
	canonicalize_class_pair(c0, c1);
}

bool genotype::supersiteIsAmbiguous(int ss_idx) const {
	uint8_t h0 = SUPERSITE_CODE_MISSING;
	uint8_t h1 = SUPERSITE_CODE_MISSING;
	getSupersitePhasedGt(ss_idx, h0, h1);
	return (h0 != SUPERSITE_CODE_MISSING &&
	        h1 != SUPERSITE_CODE_MISSING &&
	        h0 != h1);
}

void genotype::snapshotSupersitePhasedGts(const std::vector<SuperSite>& super_sites_ref,
                                          const std::vector<int>& super_site_var_index_ref) {
	const size_t required = super_sites_ref.size() * 2u;
	if (ss_phased_gts.size() != required) {
		ss_phased_gts.assign(required, SUPERSITE_CODE_MISSING);
	}
	for (size_t ss_idx = 0; ss_idx < super_sites_ref.size(); ++ss_idx) {
		const SuperSite& ss = super_sites_ref[ss_idx];
		uint8_t h0 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 0);
		uint8_t h1 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 1);
		const size_t offset = supersite_pair_offset(static_cast<int>(ss_idx));
		ss_phased_gts[offset] = h0;
		ss_phased_gts[offset + 1] = h1;
	}
}

void genotype::snapshotSupersiteObservedGts(const std::vector<SuperSite>& super_sites_ref,
                                            const std::vector<int>& super_site_var_index_ref) {
	const size_t required = super_sites_ref.size() * 2u;
	if (ss_observed_gts.size() != required) {
		ss_observed_gts.assign(required, SUPERSITE_CODE_MISSING);
	}
	for (size_t ss_idx = 0; ss_idx < super_sites_ref.size(); ++ss_idx) {
		const SuperSite& ss = super_sites_ref[ss_idx];
		uint8_t c0 = SUPERSITE_CODE_REF;
		uint8_t c1 = SUPERSITE_CODE_REF;
		resolveSupersiteClasses(*this, ss, super_site_var_index_ref, c0, c1);
		canonicalize_class_pair(c0, c1);
		const size_t offset = supersite_pair_offset(static_cast<int>(ss_idx));
		ss_observed_gts[offset] = c0;
		ss_observed_gts[offset + 1] = c1;
	}
}

genotype::SuperSiteContext genotype::getSuperSiteContext(unsigned int locus) const {
	SuperSiteContext ctx;
	if (!super_sites || !locus_to_super_idx) return ctx;
	if (locus >= locus_to_super_idx->size()) return ctx;
	ctx.ss_idx = (*locus_to_super_idx)[locus];
	if (ctx.ss_idx < 0) return ctx;
	ctx.is_member = true;
	const SuperSite& ss = (*super_sites)[ctx.ss_idx];
	ctx.is_anchor = (static_cast<uint32_t>(locus) == ss.global_site_id);
	ctx.has_het = supersiteHasFlag(ctx.ss_idx, SS_FLAG_HET);
	ctx.has_sca = supersiteHasFlag(ctx.ss_idx, SS_FLAG_SCA);
	ctx.all_missing = supersiteHasFlag(ctx.ss_idx, SS_FLAG_ALL_MIS);
	return ctx;
}

void genotype::make(vector < unsigned char > & DipSampled, vector < float > & CurrentMissingProbabilities) {
	if (DipSampled.size() != n_segments) {
		std::fprintf(stderr,
		             "[MAKE_OOB] sample=%s DipSampled.size()=%zu n_segments=%u\n",
		             name.c_str(), DipSampled.size(), n_segments);
		std::fflush(stderr);
		std::abort();
	}
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		if (!DIP_GET(Diplotypes[s], DipSampled[s])) {
			std::fprintf(stderr,
			             "[MAKE_OOB] sample=%s seg=%u dipcode=%u not in mask dipcount=%u mask=0x%016llx\n",
			             name.c_str(), s, static_cast<unsigned>(DipSampled[s]), countDiplotypes(Diplotypes[s]),
			             static_cast<unsigned long long>(Diplotypes[s]));
		std::fflush(stderr);
		std::abort();
	}
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			// Phase 3: Check if this is a missing supersite anchor
			int ss_idx = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;
			
			if (ss_idx >= 0 && anchor_has_missing && (*anchor_has_missing)[ss_idx] && SC) {
				// This is a supersite with missing data
				const SuperSite& ss = (*super_sites)[ss_idx];

				if (vabs == ss.global_site_id) {
					// Anchor: sample multivariant and project to all splits
					const int C = static_cast<int>(ss.n_classes); // 1 (REF) + n_alts
					uint32_t offset = supersite_sc_offset ? (*supersite_sc_offset)[ss_idx] : 0;

					// Sample one class per haplotype from multivariant, but couple
					// the draw to the ALT mass only (REF implicit), to mirror the
					// biallelic "r <= p_alt" semantics.
					auto sample_supersite_class = [&](int hap_lane) -> uint8_t {
						if (C <= 1) return 0u; // Only REF available

						// Row pointer for this haplotype: [REF, ALT1..ALTn]
						const float* row = &(*SC)[offset + hap_lane * C];

						float alt_sum = 0.0f;
						for (int c = 1; c < C; ++c) {
							alt_sum += row[c];
						}

						const float r = static_cast<float>(sample_rng.getDouble());

						// If no ALT mass (or r falls outside ALT CDF), choose REF (class 0)
						if (alt_sum <= 0.0f || r > alt_sum) {
							return 0u;
						}

						// Walk ALT-only CDF (classes 1..C-1). This keeps REF implicit and
						// aligns with the bial path, where the draw is compared against p_alt.
						float cumsum = 0.0f;
						for (int c = 1; c < C; ++c) {
							cumsum += row[c];
							if (r <= cumsum) {
								return static_cast<uint8_t>(c);
							}
						}
						// Fallback for rounding error: last ALT
						return static_cast<uint8_t>(C - 1);
					};

					uint8_t h0 = sample_supersite_class(hap0);
					uint8_t h1 = sample_supersite_class(hap1);

                    // Project to splits from sampled h0/h1: class 0=REF, 1..n_alts=ALT1..ALTn
					// Iterate over all member variants and set based on sampled class
                    int alt_count_h0 = 0;
                    int alt_count_h1 = 0;
                    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
                        unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
                        uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
                        unsigned char& split_byte = Variants[DIV2(split_vabs)];
                        
                        // Set hap0
                        if (h0 == alt_class) {
                            VAR_SET_HAP0(MOD2(split_vabs), split_byte);
                            alt_count_h0++;
                        } else {
                            VAR_CLR_HAP0(MOD2(split_vabs), split_byte);
                        }
                        
                        // Set hap1
                        if (h1 == alt_class) {
                            VAR_SET_HAP1(MOD2(split_vabs), split_byte);
                            alt_count_h1++;
                        } else {
                            VAR_CLR_HAP1(MOD2(split_vabs), split_byte);
                        }
                    }
                    
                    // Persist the sampled classes for this epoch as the current h0/h1 pair
                    setSupersitePhasedGt(ss_idx, h0, h1);

					// Do not attempt range skip: supersite members are not guaranteed contiguous
					// Count a single missing event for the supersite and proceed; siblings are skipped below
					m++;  // One missing event per supersite
					continue;
				}
				// Sibling: skip, already handled by anchor
				{
					bool is_member = false;
					for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
						if ((*super_site_var_index)[ss.var_start + ai] == vabs) { is_member = true; break; }
					}
					if (is_member && (int)vabs != (int)ss.global_site_id) {
						continue;
					}
				}
			}
			
			// Normal biallelic missing site
			if (VAR_GET_MIS(MOD2(vabs), Variants[DIV2(vabs)])) {
				if (haploid) {
					float p00 = (1.0f - CurrentMissingProbabilities[m*HAP_NUMBER+hap0]) * (1.0f - CurrentMissingProbabilities[m*HAP_NUMBER+hap1]);
					float p11 = (CurrentMissingProbabilities[m*HAP_NUMBER+hap0]) * (CurrentMissingProbabilities[m*HAP_NUMBER+hap1]);
					if (sample_rng.getDouble()<= (p11/(p00+p11))) {
						VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
						VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
					} else {
						VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
						VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
					}
				} else {
					float r0 = sample_rng.getDouble();
					float p0 = CurrentMissingProbabilities[m*HAP_NUMBER+hap0];
					(r0<=p0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
					
					float r1 = sample_rng.getDouble();
					float p1 = CurrentMissingProbabilities[m*HAP_NUMBER+hap1];
					(r1<=p1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);

				}
				m++;
			}
				if (VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)])) {
					// Recompute ss_idx for this specific variant
				int ss_idx_amb = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;
				
				// Check if this is a heterozygous supersite
				if (ss_idx_amb >= 0 && super_sites && super_site_var_index) {
					const SuperSite& ss = (*super_sites)[ss_idx_amb];
					if (vabs == ss.global_site_id) {
						// This is a heterozygous supersite anchor
                    // Determine h0/h1 (sampled classes for this epoch) from lanes and immutable c0/c1
                    
							// Get the current allele codes from the immutable supersite snapshot
                    uint8_t current_c0 = SUPERSITE_CODE_MISSING;
                    uint8_t current_c1 = SUPERSITE_CODE_MISSING;
                    getSupersiteObservedGt(ss_idx_amb, current_c0, current_c1);
						
						// The Ambiguous array tells us the phase orientation for this site
						// Use it to determine which allele goes to which lane
						unsigned char amb_code = Ambiguous[a];
						
						// Determine which allele class each sampled lane should get
						// If HAP_GET(amb_code, hap) == 0, this lane gets c0's allele
						// If HAP_GET(amb_code, hap) == 1, this lane gets c1's allele
                    uint8_t h0 = HAP_GET(amb_code, hap0) ? current_c1 : current_c0;
                    uint8_t h1 = HAP_GET(amb_code, hap1) ? current_c1 : current_c0;

	                    // Project to all splits based on sampled classes
                    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
                        unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
                        uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
                        
                        // Set hap0
                        if (h0 == alt_class) {
                            VAR_SET_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
                        } else {
                            VAR_CLR_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
                        }
                        
                        // Set hap1
                        if (h1 == alt_class) {
                            VAR_SET_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
                        } else {
                            VAR_CLR_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
                        }
                    }

                    // Persist sampled h0/h1 for this supersite anchor
                    setSupersitePhasedGt(ss_idx_amb, h0, h1);
						
							// Do not perform range skip: supersite members are not guaranteed contiguous
							a++;  // One ambiguous event per supersite
							continue;
					} else {
						// This is a sibling of a supersite (not the anchor). Members may be non-consecutive
						bool is_member = false;
						for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
							if ((*super_site_var_index)[ss.var_start + ai] == vabs) { is_member = true; break; }
						}
						if (is_member) {
							// Skip this sibling - already handled via anchor
							continue;
						}
					}
				}
				
				// Normal biallelic heterozygous site
				bool hap0_bit = HAP_GET(Ambiguous[a], hap0);
				bool hap1_bit = HAP_GET(Ambiguous[a], hap1);
				hap0_bit ? VAR_SET_HAP0(MOD2(vabs), Variants[DIV2(vabs)]) : VAR_CLR_HAP0(MOD2(vabs), Variants[DIV2(vabs)]);
				hap1_bit ? VAR_SET_HAP1(MOD2(vabs), Variants[DIV2(vabs)]) : VAR_CLR_HAP1(MOD2(vabs), Variants[DIV2(vabs)]);
				a++;
			}
		}

	}
}

void genotype::make(vector < unsigned char > & DipSampled) {
	if (DipSampled.size() != n_segments) {
		std::fprintf(stderr,
		             "[MAKE_OOB] sample=%s DipSampled.size()=%zu n_segments=%u (no missing probs)\n",
		             name.c_str(), DipSampled.size(), n_segments);
		std::fflush(stderr);
		std::abort();
	}
	bool needs_supersite_projection = false;
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		if (!DIP_GET(Diplotypes[s], DipSampled[s])) {
			std::fprintf(stderr,
			             "[MAKE_OOB] sample=%s seg=%u dipcode=%u not in mask dipcount=%u mask=0x%016llx (no missing probs)\n",
			             name.c_str(), s, static_cast<unsigned>(DipSampled[s]), countDiplotypes(Diplotypes[s]),
			             static_cast<unsigned long long>(Diplotypes[s]));
			std::fflush(stderr);
			std::abort();
		}
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			int ss_idx_mis = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;

			if (VAR_GET_MIS(MOD2(vabs), Variants[DIV2(vabs)])) {
				if (ss_idx_mis >= 0 && super_sites && super_site_var_index) {
					const SuperSite& ss = (*super_sites)[ss_idx_mis];

					// Sibling: anchor will project
					if (static_cast<uint32_t>(vabs) != ss.global_site_id) {
						continue;
					}

					// Anchor: pick class per hap from aggregated SC posteriors
					const int C = static_cast<int>(ss.n_classes);
					auto avg_supersite_class = [&](int lane) -> uint8_t {
						if (ProbSuperClass.empty() || sc_storage_events == 0) {
							return SUPERSITE_CODE_REF;
						}
						float best_p = -1.0f;
						int best_c = 0;
						for (int c = 0; c < C; ++c) {
							const size_t idx = supersite_class_index(ss_idx_mis, lane, c);
							const float p = ProbSuperClass[idx] / static_cast<float>(sc_storage_events);
							if (p > best_p) {
								best_p = p;
								best_c = c;
							}
						}
						return static_cast<uint8_t>(best_c);
					};

					uint8_t final_c0 = avg_supersite_class(hap0);
					uint8_t final_c1 = avg_supersite_class(hap1);

					for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
						unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
						uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
						unsigned char& vbyte = Variants[DIV2(split_vabs)];

						if (final_c0 == alt_class) VAR_SET_HAP0(MOD2(split_vabs), vbyte);
						else                       VAR_CLR_HAP0(MOD2(split_vabs), vbyte);

						if (final_c1 == alt_class) VAR_SET_HAP1(MOD2(split_vabs), vbyte);
						else                       VAR_CLR_HAP1(MOD2(split_vabs), vbyte);
					}

					setSupersitePhasedGt(ss_idx_mis, final_c0, final_c1);
					needs_supersite_projection = true;
					continue;
				}
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
			bool is_amb = VAR_GET_AMB(MOD2(vabs), Variants[DIV2(vabs)]);

				if (is_amb) {
					// Recompute ss_idx for this specific variant
				int ss_idx_amb = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;

				// Check if this is a heterozygous supersite
				if (ss_idx_amb >= 0 && super_sites && super_site_var_index) {
					const SuperSite& ss = (*super_sites)[ss_idx_amb];

					if (vabs == ss.global_site_id) {
						// This is a heterozygous supersite anchor
						// Determine h0/h1 (sampled classes for this epoch) from lanes and immutable c0/c1

						// Get the current allele codes from the immutable supersite snapshot
						uint8_t current_c0 = SUPERSITE_CODE_MISSING;
						uint8_t current_c1 = SUPERSITE_CODE_MISSING;
						getSupersiteObservedGt(ss_idx_amb, current_c0, current_c1);

						// The Ambiguous array tells us the phase orientation for this site
						// Use it to determine which allele goes to which lane
						unsigned char amb_code = Ambiguous[a];

						// Determine which allele class each sampled lane should get
						// If HAP_GET(amb_code, hap) == 0, this lane gets c0's allele
						// If HAP_GET(amb_code, hap) == 1, this lane gets c1's allele
						uint8_t h0 = HAP_GET(amb_code, hap0) ? current_c1 : current_c0;
						uint8_t h1 = HAP_GET(amb_code, hap1) ? current_c1 : current_c0;

						// Project to all splits based on sampled classes
						for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
							unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
							uint8_t alt_class = ai + 1;

							// Set hap0
							if (h0 == alt_class) {
								VAR_SET_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							} else {
								VAR_CLR_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							}

							// Set hap1
							if (h1 == alt_class) {
								VAR_SET_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							} else {
								VAR_CLR_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							}
						}
						// Persist sampled h0/h1 for this supersite anchor
						setSupersitePhasedGt(ss_idx_amb, h0, h1);

						// Do not perform range skip: supersite members are not guaranteed contiguous
						a++;  // One ambiguous event per supersite
						continue;
					} else {
						// This is a sibling of a supersite (not the anchor). Members may be non-consecutive
						bool is_member = false;
						for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
							if ((*super_site_var_index)[ss.var_start + ai] == vabs) { is_member = true; break; }
						}
						if (is_member) {
							// Skip this sibling - already handled via anchor
							continue;
						}
					}
				}

					bool hap0_bit = HAP_GET(Ambiguous[a], hap0);
				bool hap1_bit = HAP_GET(Ambiguous[a], hap1);
				hap0_bit ? VAR_SET_HAP0(MOD2(vabs), Variants[DIV2(vabs)]) : VAR_CLR_HAP0(MOD2(vabs), Variants[DIV2(vabs)]);
				hap1_bit ? VAR_SET_HAP1(MOD2(vabs), Variants[DIV2(vabs)]) : VAR_CLR_HAP1(MOD2(vabs), Variants[DIV2(vabs)]);
				a++;
			}
		}
	}
	if (needs_supersite_projection && super_sites) {
		projectSupersites();
	}
}

// Phase 3: Set supersite context for multivariant imputation
void genotype::setSuperSiteContext(
const std::vector<SuperSite>* _super_sites,
const std::vector<int>* _locus_to_super_idx,
const std::vector<int>* _super_site_var_index,
const std::vector<float>* _SC,
const std::vector<bool>* _anchor_has_missing,
const std::vector<uint32_t>* _supersite_sc_offset)
{
super_sites = _super_sites;
locus_to_super_idx = _locus_to_super_idx;
super_site_var_index = _super_site_var_index;
SC = _SC;
anchor_has_missing = _anchor_has_missing;
supersite_sc_offset = _supersite_sc_offset;
	if (super_sites) {
		const size_t required = super_sites->size() * 2u;
		if (ss_phased_gts.size() != required) {
			ss_phased_gts.assign(required, SUPERSITE_CODE_MISSING);
		}
	}
}

void genotype::projectSupersites() {
	if (!super_sites || !locus_to_super_idx || !super_site_var_index) {
		return; // No supersite data available
	}

	// Iterate through all supersites to perform post-HMM projection
	for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
		const SuperSite& ss = (*super_sites)[ss_idx];
		int anchor_locus = static_cast<int>(ss.global_site_id);

		// Skip if anchor is missing
		if (VAR_GET_MIS(MOD2(anchor_locus), Variants[DIV2(anchor_locus)])) {
			continue;
		}

        // Determine current sampled classes (h0/h1) for this anchor from stored pair
        // Note: getSupersitePhasedGt returns the current epoch's h0/h1; c0/c1 remain immutable in emissions.
        uint8_t h0 = SUPERSITE_CODE_MISSING;
        uint8_t h1 = SUPERSITE_CODE_MISSING;
        getSupersitePhasedGt(static_cast<int>(ss_idx), h0, h1);

        // Project to splits from sampled h0/h1: class 0=REF, 1..n_alts=ALT1..ALTn
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            unsigned int split_locus = (*super_site_var_index)[ss.var_start + ai];
            uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.

            // Set hap0 for this split
            if (h0 == alt_class) {
                VAR_SET_HAP0(MOD2(split_locus), Variants[DIV2(split_locus)]);
            } else {
                VAR_CLR_HAP0(MOD2(split_locus), Variants[DIV2(split_locus)]);
            }

            // Set hap1 for this split
            if (h1 == alt_class) {
                VAR_SET_HAP1(MOD2(split_locus), Variants[DIV2(split_locus)]);
            } else {
                VAR_CLR_HAP1(MOD2(split_locus), Variants[DIV2(split_locus)]);
            }
        }

	}
}

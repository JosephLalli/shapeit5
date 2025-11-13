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
#include <models/supersite_trace_utils.h>
#include <models/super_site_accessor.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>

// K-inflation debug instrumentation
namespace debug {
std::string SUPERDEBUG_SAMPLENAME = "";
int SUPERDEBUG_BP = 0;

void load_debug_settings() {
    const char* sample_env = std::getenv("SHAPEIT5_SUPERDEBUG_SAMPLENAME");
    if (sample_env) {
        SUPERDEBUG_SAMPLENAME = std::string(sample_env);
    }
    const char* bp_env = std::getenv("SHAPEIT5_SUPERDEBUG_BP");
    if (bp_env) {
        SUPERDEBUG_BP = std::atoi(bp_env);
    }
}

void print_supersite_state(const genotype* G, const SuperSite& ss, const std::vector<int>& super_site_var_index, const std::string& context) {
    if (G->name != SUPERDEBUG_SAMPLENAME || (int)ss.global_site_id != SUPERDEBUG_BP) return;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "[SUPERDEBUG] Sample=" << G->name << " Pos=" << ss.global_site_id << " Context='" << context << "'" << std::endl;

    std::string hap0_str = "  HAP0: ";
    std::string hap1_str = "  HAP1: ";

    for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
        unsigned int split_locus = super_site_var_index[ss.var_start + ai];
        unsigned char v_byte = G->Variants[DIV2(split_locus)];
        hap0_str += std::to_string(VAR_GET_HAP0(MOD2(split_locus), v_byte)) + " ";
        hap1_str += std::to_string(VAR_GET_HAP1(MOD2(split_locus), v_byte)) + " ";
    }
    std::cout << hap0_str << std::endl;
    std::cout << hap1_str << std::endl;
    std::cout.unsetf(std::ios_base::floatfield);
}
} // namespace debug
// End K-inflation debug

#if !defined(__OPTIMIZE__)
#include <cstdlib>
#include <iostream>

namespace {
template <typename Vec>
inline void supersite_debug_require(const Vec* vec, size_t required, const char* vec_name, const char* func) {
	if (!vec) {
		std::cerr << "[supersite debug] " << func << ": " << vec_name << " pointer is null" << std::endl;
		std::abort();
	}
	const size_t size = vec->size();
	if (size < required) {
		std::cerr << "[supersite debug] " << func << ": " << vec_name
				  << " size=" << size << " required>=" << required << std::endl;
		std::abort();
	}
}

inline void supersite_debug_check_var_count(uint16_t var_count, const char* func) {
	if (var_count == 0) {
		std::cerr << "[supersite debug] " << func << ": supersite var_count is zero" << std::endl;
		std::abort();
	}
}
} // namespace
#endif

using namespace std;

genotype::genotype(unsigned int _index) {
	if (debug::SUPERDEBUG_BP == 0) debug::load_debug_settings();
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
	vector < uint8_t > ().swap(supersite_flags);
}

bool genotype::supersiteHasFlag(int ss_idx, uint8_t flag) const {
	return ss_idx >= 0 &&
	       ss_idx < static_cast<int>(supersite_flags.size()) &&
	       (supersite_flags[ss_idx] & flag);
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
	for (unsigned int s = 0, vabs = 0, a = 0, m = 0 ; s < n_segments ; s ++) {
		unsigned char hap0 = DIP_HAP0(DipSampled[s]);
		unsigned char hap1 = DIP_HAP1(DipSampled[s]);
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel++, vabs++) {
			// Phase 3: Check if this is a missing supersite anchor
			int ss_idx = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;
			
#if !defined(__OPTIMIZE__)
			if (ss_idx >= 0) {
				const char* dbg_func = __func__;
				supersite_debug_require(super_sites, static_cast<size_t>(ss_idx) + 1, "super_sites", dbg_func);
				supersite_debug_require(locus_to_super_idx, static_cast<size_t>(vabs) + 1, "locus_to_super_idx", dbg_func);
				supersite_debug_require(super_site_var_index, static_cast<size_t>((*super_sites)[ss_idx].var_start) + (*super_sites)[ss_idx].var_count, "super_site_var_index", dbg_func);
				supersite_debug_check_var_count((*super_sites)[ss_idx].var_count, dbg_func);
				if (anchor_has_missing) {
					supersite_debug_require(anchor_has_missing, static_cast<size_t>(ss_idx) + 1, "anchor_has_missing", dbg_func);
				}
				if (SC && anchor_has_missing && (*anchor_has_missing)[ss_idx]) {
					std::cout << "[supersite debug] " << dbg_func << ": checking SC for ss_idx=" << ss_idx 
					          << " n_classes=" << (*super_sites)[ss_idx].n_classes 
					          << " anchor_has_missing=" << (*anchor_has_missing)[ss_idx] 
					          << " SC.size()=" << SC->size() << std::endl;
					std::cout.flush();
					supersite_debug_require(supersite_sc_offset, static_cast<size_t>(ss_idx) + 1, "supersite_sc_offset", dbg_func);
					const size_t offset_debug = (*supersite_sc_offset)[ss_idx];
					const size_t required_debug = offset_debug + static_cast<size_t>(HAP_NUMBER) * static_cast<size_t>((*super_sites)[ss_idx].n_classes);
					supersite_debug_require(SC, required_debug, "SC", dbg_func);
				}
			}
#endif

			if (ss_idx >= 0 && anchor_has_missing && (*anchor_has_missing)[ss_idx] && SC) {
				// This is a supersite with missing data
				const SuperSite& ss = (*super_sites)[ss_idx];
				
				if (vabs == ss.global_site_id) {
					// Anchor: sample multivariant and project to all splits
					int C = (int)ss.n_classes;
					uint32_t offset = supersite_sc_offset ? (*supersite_sc_offset)[ss_idx] : 0;
					
					// Sample one class per haplotype from multivariant
					// SC[offset + hap*C + c] = P(class_c | hap)
					uint8_t class0 = 0, class1 = 0;

                    // HYPOTHESIS 3 DEBUGGING
                    if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                        debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo3: Enter missing anchor");
                        std::cout << "[SUPERDEBUG] H3: Imputing missing data at Pos=" << ss.global_site_id << " for sample " << this->name << std::endl;
                        std::cout << "[SUPERDEBUG] H3: SC offset=" << offset << " hap0_idx=" << hap0 << " hap1_idx=" << hap1 << " n_classes=" << C << std::endl;
                        std::string sc_probs0 = "  SC[hap0]: ";
                        std::string sc_probs1 = "  SC[hap1]: ";
                        for (int c = 0; c < C; ++c) {
                            sc_probs0 += std::to_string((*SC)[offset + hap0 * C + c]) + " ";
                            sc_probs1 += std::to_string((*SC)[offset + hap1 * C + c]) + " ";
                        }
                        std::cout << sc_probs0 << std::endl;
                        std::cout << sc_probs1 << std::endl;
                    }
					
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
                    if (cumsum0 < r0) {
                        class0 = C > 0 ? static_cast<uint8_t>(C - 1) : 0;
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
                    if (cumsum1 < r1) {
                        class1 = C > 0 ? static_cast<uint8_t>(C - 1) : 0;
                    }

                    // HYPOTHESIS 3 DEBUGGING
                    if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                        std::cout << "[SUPERDEBUG] H3: r0=" << r0 << " -> sampled_class0=" << (int)class0 << std::endl;
                        std::cout << "[SUPERDEBUG] H3: r1=" << r1 << " -> sampled_class1=" << (int)class1 << std::endl;
                    }
                    if (supersite_trace_enabled()) {
                        supersite_trace_log("[SupersiteSample] sample=%s ss_idx=%d anchor=%u C=%d class0=%u class1=%u offset=%u\n",
                                            this->name.c_str(),
                                            ss_idx,
                                            ss.global_site_id,
                                            C,
                                            static_cast<unsigned>(class0),
                                            static_cast<unsigned>(class1),
                                            offset);
                    }
					
					// Project to splits: class 0=REF, 1..n_alts=ALT1..ALTn
					// Iterate over all member variants and set based on sampled class
					int alt_count_h0 = 0;
					int alt_count_h1 = 0;
					for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
						unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
						uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
						
						// Set hap0
						if (class0 == alt_class) {
							VAR_SET_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							alt_count_h0++;
						} else {
							VAR_CLR_HAP0(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						}
						
						// Set hap1
						if (class1 == alt_class) {
							VAR_SET_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
							alt_count_h1++;
						} else {
							VAR_CLR_HAP1(MOD2(split_vabs), Variants[DIV2(split_vabs)]);
						}
					}
					if (supersite_trace_enabled()) {
						supersite_trace_log("[SupersiteSample] projection sample=%s ss_idx=%d alt_count_h0=%d alt_count_h1=%d\n",
						                    this->name.c_str(),
						                    ss_idx,
						                    alt_count_h0,
						                    alt_count_h1);
					}

                    // HYPOTHESIS 3 DEBUGGING
                    if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                        debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo3: After projection in missing");
                    }
					if (supersite_debug::guard_checks_enabled()) {
						int alt_count_h0 = 0;
						int alt_count_h1 = 0;
						for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
							unsigned int split_vabs = (*super_site_var_index)[ss.var_start + ai];
							unsigned char vbyte = Variants[DIV2(split_vabs)];
							alt_count_h0 += VAR_GET_HAP0(MOD2(split_vabs), vbyte) ? 1 : 0;
							alt_count_h1 += VAR_GET_HAP1(MOD2(split_vabs), vbyte) ? 1 : 0;
						}
						if (alt_count_h0 > 1 || alt_count_h1 > 1 || (class0 > 0 && alt_count_h0 == 0) || (class1 > 0 && alt_count_h1 == 0)) {
							std::string msg = "projection mismatch sample=" + name +
								" ss_idx=" + std::to_string(ss_idx) +
								" class0=" + std::to_string(class0) +
								" class1=" + std::to_string(class1) +
								" alt_h0=" + std::to_string(alt_count_h0) +
								" alt_h1=" + std::to_string(alt_count_h1);
							supersite_debug::report_guard_violation("genotype::make", msg.c_str());
							assert(alt_count_h0 <= 1);
							assert(alt_count_h1 <= 1);
							assert(class0 == 0 || alt_count_h0 > 0);
							assert(class1 == 0 || alt_count_h1 > 0);
						}
					}
					
					// Do not attempt range skip: members may be non-consecutive after MAC pruning
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
						if (supersite_trace_enabled()) {
							supersite_trace_log("[SupersiteSample] skip sibling sample=%s ss_idx=%d locus=%u\n",
							                    this->name.c_str(), ss_idx, vabs);
						}
						continue;
					}
				}
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
				// Recompute ss_idx for this specific variant
				int ss_idx_amb = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1;
				
				// Check if this is a heterozygous supersite
				if (ss_idx_amb >= 0 && super_sites && super_site_var_index) {
					const SuperSite& ss = (*super_sites)[ss_idx_amb];
					
					if (vabs == ss.global_site_id) {
						// This is a heterozygous supersite anchor
						// Determine which allele class each sampled lane represents
						
                        // HYPOTHESIS 2 DEBUGGING
                        if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                            debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: Enter AMB block");
                        }

						// Get the current allele codes from the Variants data
						uint8_t current_c0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
						uint8_t current_c1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);
						
						// The Ambiguous array tells us the phase orientation for this site
						// Use it to determine which allele goes to which lane
						unsigned char amb_code = Ambiguous[a];
						
						// Determine which allele class each sampled lane should get
						// If HAP_GET(amb_code, hap) == 0, this lane gets c0's allele
						// If HAP_GET(amb_code, hap) == 1, this lane gets c1's allele
						uint8_t class0 = HAP_GET(amb_code, hap0) ? current_c1 : current_c0;
						uint8_t class1 = HAP_GET(amb_code, hap1) ? current_c1 : current_c0;

                        // HYPOTHESIS 2 DEBUGGING
                        if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                            std::cout << "[SUPERDEBUG] H2: Resolving het at Pos=" << ss.global_site_id << " for sample " << this->name << std::endl;
                            std::cout << "[SUPERDEBUG] H2: current_c0=" << (int)current_c0 << " current_c1=" << (int)current_c1 << " amb_code=" << (int)amb_code << std::endl;
                            std::cout << "[SUPERDEBUG] H2: hap0_lane=" << (int)hap0 << " hap1_lane=" << (int)hap1 << std::endl;
                            std::cout << "[SUPERDEBUG] H2: sampled_class0=" << (int)class0 << " sampled_class1=" << (int)class1 << std::endl;
                        }
						
						// Project to all splits based on sampled classes
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

                        // HYPOTHESIS 2 DEBUGGING
                        if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                            debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After projection in AMB");
                        }
						
							// Do not perform range skip: members may be non-consecutive after MAC pruning
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

		// Determine which allele class is on each haplotype from the anchor's perspective
		uint8_t class0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
		uint8_t class1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);

		// Project to splits: class 0=REF, 1..n_alts=ALT1..ALTn
		for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
			unsigned int split_locus = (*super_site_var_index)[ss.var_start + ai];
			uint8_t alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.

			// Set hap0 for this split
			if (class0 == alt_class) {
				VAR_SET_HAP0(MOD2(split_locus), Variants[DIV2(split_locus)]);
			} else {
				VAR_CLR_HAP0(MOD2(split_locus), Variants[DIV2(split_locus)]);
			}

			// Set hap1 for this split
			if (class1 == alt_class) {
				VAR_SET_HAP1(MOD2(split_locus), Variants[DIV2(split_locus)]);
			} else {
				VAR_CLR_HAP1(MOD2(split_locus), Variants[DIV2(split_locus)]);
			}
		}
	}
}

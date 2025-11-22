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
	vector < uint8_t > ().swap(supersite_class_pairs);
	vector < uint8_t > ().swap(supersite_class_pairs_base);
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

void genotype::setSupersiteClassPair(int ss_idx, uint8_t class0, uint8_t class1) {
	if (ss_idx < 0) return;
	const size_t required_pairs = super_sites ? super_sites->size() : supersite_class_pairs.size() / 2u;
	const size_t required = required_pairs * 2u;
	if (supersite_class_pairs.size() != required) {
		supersite_class_pairs.assign(required, SUPERSITE_CODE_MISSING);
	}
	const size_t offset = supersite_pair_offset(ss_idx);
	if (offset + 1 < supersite_class_pairs.size()) {
		supersite_class_pairs[offset] = safe_class_code(class0);
		supersite_class_pairs[offset + 1] = safe_class_code(class1);
	}
}

void genotype::getSupersiteClassPair(int ss_idx, uint8_t& class0, uint8_t& class1) const {
	if (ss_idx >= 0) {
		const size_t offset = supersite_pair_offset(ss_idx);
		if (offset + 1 < supersite_class_pairs.size()) {
			class0 = supersite_class_pairs[offset];
			class1 = supersite_class_pairs[offset + 1];
			return;
		}
	}
    // Fallback: derive from current hap bits if storage not initialized
    // Note: when initialized, this pair represents the current epoch's sampled h0/h1.
    // Immutable c0/c1 are used in emissions (see SiteView.sample_class0/1).
	class0 = class1 = SUPERSITE_CODE_MISSING;
	if (!super_sites || !super_site_var_index || ss_idx < 0 || ss_idx >= static_cast<int>(super_sites->size())) return;
	const SuperSite& ss = (*super_sites)[ss_idx];
	class0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
	class1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);
}

void genotype::getSupersiteBaseClassPair(int ss_idx, uint8_t& class0, uint8_t& class1) const {
	class0 = class1 = SUPERSITE_CODE_MISSING;
	if (ss_idx >= 0) {
		const size_t offset = supersite_pair_offset(ss_idx);
		if (offset + 1 < supersite_class_pairs_base.size()) {
			class0 = supersite_class_pairs_base[offset];
			class1 = supersite_class_pairs_base[offset + 1];
            canonicalize_class_pair(class0, class1);
			return;
		}
	}
	if (!super_sites || !super_site_var_index || ss_idx < 0 || ss_idx >= static_cast<int>(super_sites->size())) return;
	const SuperSite& ss = (*super_sites)[ss_idx];
	class0 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 0);
	class1 = getSampleSuperSiteAlleleCode(this, ss, *super_site_var_index, 1);
    canonicalize_class_pair(class0, class1);
}

bool genotype::supersiteIsAmbiguous(int ss_idx) const {
	uint8_t class0 = SUPERSITE_CODE_MISSING;
	uint8_t class1 = SUPERSITE_CODE_MISSING;
	getSupersiteClassPair(ss_idx, class0, class1);
	return (class0 != SUPERSITE_CODE_MISSING &&
	        class1 != SUPERSITE_CODE_MISSING &&
	        class0 != class1);
}

void genotype::snapshotSupersiteClasses(const std::vector<SuperSite>& super_sites_ref,
                                        const std::vector<int>& super_site_var_index_ref) {
	const size_t required = super_sites_ref.size() * 2u;
	if (supersite_class_pairs.size() != required) {
		supersite_class_pairs.assign(required, SUPERSITE_CODE_MISSING);
	}
	for (size_t ss_idx = 0; ss_idx < super_sites_ref.size(); ++ss_idx) {
		const SuperSite& ss = super_sites_ref[ss_idx];
		uint8_t class0 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 0);
		uint8_t class1 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 1);
		const size_t offset = supersite_pair_offset(static_cast<int>(ss_idx));
		supersite_class_pairs[offset] = class0;
		supersite_class_pairs[offset + 1] = class1;
	}
}

void genotype::snapshotSupersiteBaseClasses(const std::vector<SuperSite>& super_sites_ref,
                                            const std::vector<int>& super_site_var_index_ref) {
	const size_t required = super_sites_ref.size() * 2u;
	if (supersite_class_pairs_base.size() != required) {
		supersite_class_pairs_base.assign(required, SUPERSITE_CODE_MISSING);
	}
	for (size_t ss_idx = 0; ss_idx < super_sites_ref.size(); ++ss_idx) {
		const SuperSite& ss = super_sites_ref[ss_idx];
		uint8_t class0 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 0);
		uint8_t class1 = getSampleSuperSiteAlleleCode(this, ss, super_site_var_index_ref, 1);
        canonicalize_class_pair(class0, class1);
		const size_t offset = supersite_pair_offset(static_cast<int>(ss_idx));
		supersite_class_pairs_base[offset] = class0;
		supersite_class_pairs_base[offset + 1] = class1;
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
                    // Here, h0/h1 are the sampled supersite classes for this epoch (mutable);
                    // c0/c1 remain the immutable site classes used for emissions.
                    uint8_t h0 = 0, h1 = 0;

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
					
                    // Sample h0 for hap0
                    float r0 = rng.getDouble();
                    float cumsum0 = 0.0f;
                    for (int c = 0; c < C; ++c) {
                        cumsum0 += (*SC)[offset + hap0 * C + c];
                        if (r0 <= cumsum0) {
                            h0 = c;
                            break;
                        }
                    }
                    if (cumsum0 < r0) {
                        h0 = C > 0 ? static_cast<uint8_t>(C - 1) : 0;
                    }
                    
                    // Sample h1 for hap1
                    float r1 = rng.getDouble();
                    float cumsum1 = 0.0f;
                    for (int c = 0; c < C; ++c) {
                        cumsum1 += (*SC)[offset + hap1 * C + c];
                        if (r1 <= cumsum1) {
                            h1 = c;
                            break;
                        }
                    }
					if (cumsum1 < r1) {
                        h1 = C > 0 ? static_cast<uint8_t>(C - 1) : 0;
                    }

                    // HYPOTHESIS 3 DEBUGGING
                    if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && (int)ss.global_site_id == debug::SUPERDEBUG_BP) {
                        std::cout << "[SUPERDEBUG] H3: r0=" << r0 << " -> sampled_h0=" << (int)h0 << std::endl;
                        std::cout << "[SUPERDEBUG] H3: r1=" << r1 << " -> sampled_h1=" << (int)h1 << std::endl;
                    }
                    if (supersite_trace_enabled()) {
                        supersite_trace_log("[SupersiteSample] sample=%s ss_idx=%d anchor=%u C=%d h0=%u h1=%u offset=%u\n",
                                            this->name.c_str(),
                                            ss_idx,
                                            ss.global_site_id,
                                            C,
                                            static_cast<unsigned>(h0),
                                            static_cast<unsigned>(h1),
                                            offset);
                    }
					
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

                        // Keep anchor flagged missing if desired for emission semantics; clear MIS on siblings so
                        // downstream class readers can see the imputed allele.
                        if (split_vabs != ss.global_site_id) {
                            if (h0 == alt_class && h1 == alt_class) {
                                VAR_SET_HOM(MOD2(split_vabs), split_byte);
                            } else if (h0 == alt_class || h1 == alt_class) {
                                VAR_SET_HET(MOD2(split_vabs), split_byte);
                            } else {
                                VAR_SET_HOM(MOD2(split_vabs), split_byte);
                            }
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
						if (alt_count_h0 > 1 || alt_count_h1 > 1 || (h0 > 0 && alt_count_h0 == 0) || (h1 > 0 && alt_count_h1 == 0)) {
                        std::string msg = "projection mismatch sample=" + name +
                            " ss_idx=" + std::to_string(ss_idx) +
                            " h0=" + std::to_string(h0) +
                            " h1=" + std::to_string(h1) +
                            " alt_h0=" + std::to_string(alt_count_h0) +
                            " alt_h1=" + std::to_string(alt_count_h1);
                        supersite_debug::report_guard_violation("genotype::make", msg.c_str());
                        assert(alt_count_h0 <= 1);
                        assert(alt_count_h1 <= 1);
                        assert(h0 == 0 || alt_count_h0 > 0);
                        assert(h1 == 0 || alt_count_h1 > 0);
                    }
                    }
                    
                    // Persist the sampled classes for this epoch as the current h0/h1 pair
                    setSupersiteClassPair(ss_idx, h0, h1);

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
                    // Determine h0/h1 (sampled classes for this epoch) from lanes and immutable c0/c1
                    
                    const bool is_superdebug_target =
							(!debug::SUPERDEBUG_SAMPLENAME.empty() &&
							 this->name == debug::SUPERDEBUG_SAMPLENAME &&
							 static_cast<int>(ss.global_site_id) == debug::SUPERDEBUG_BP);
                        if (is_superdebug_target) {
							unsigned long dip_mask = Diplotypes[s];
							std::cout << "[SUPERDEBUG] H2: Enter AMB block segment=" << s
									  << " dip_mask=0x" << std::hex << dip_mask << std::dec
									  << " DipSampled=" << static_cast<unsigned>(DipSampled[s])
									  << " hap0_lane=" << static_cast<int>(hap0)
									  << " hap1_lane=" << static_cast<int>(hap1)
									  << std::endl;
                            debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: Enter AMB block");
                        }

						// Get the current allele codes from the immutable supersite snapshot
                    uint8_t current_c0 = SUPERSITE_CODE_MISSING;
                    uint8_t current_c1 = SUPERSITE_CODE_MISSING;
                    getSupersiteBaseClassPair(ss_idx_amb, current_c0, current_c1);
						
						// The Ambiguous array tells us the phase orientation for this site
						// Use it to determine which allele goes to which lane
						unsigned char amb_code = Ambiguous[a];
						
						// Determine which allele class each sampled lane should get
						// If HAP_GET(amb_code, hap) == 0, this lane gets c0's allele
						// If HAP_GET(amb_code, hap) == 1, this lane gets c1's allele
                    uint8_t h0 = HAP_GET(amb_code, hap0) ? current_c1 : current_c0;
                    uint8_t h1 = HAP_GET(amb_code, hap1) ? current_c1 : current_c0;

					// TRACE: Log amb_code interpretation for burn3 debugging
					if (supersite_trace_enabled() && vabs == ss.global_site_id && ss.global_site_id == 0) {
						bool hap0_bit = HAP_GET(amb_code, hap0);
						bool hap1_bit = HAP_GET(amb_code, hap1);
						std::fprintf(stderr, "[MAKE_SS_AMB] sample=%s locus=%u ss_idx=%d\n",
						             name.c_str(), vabs, ss_idx_amb);
						std::fprintf(stderr, "  sampled_lanes: hap0=%u hap1=%u\n",
						             (unsigned)hap0, (unsigned)hap1);
						std::fprintf(stderr, "  amb_code=0x%02x hap0_bit=%u hap1_bit=%u\n",
						             (unsigned)amb_code, (unsigned)hap0_bit, (unsigned)hap1_bit);
						std::fprintf(stderr, "  allele_classes: c0=%u c1=%u\n",
						             (unsigned)current_c0, (unsigned)current_c1);
						std::fprintf(stderr, "  result: h0=%u h1=%u\n",
						             (unsigned)h0, (unsigned)h1);
					}

                        // HYPOTHESIS 2 DEBUGGING
                        if (is_superdebug_target) {
                            std::cout << "[SUPERDEBUG] H2: Resolving het at Pos=" << ss.global_site_id << " for sample " << this->name << std::endl;
                            std::cout << "[SUPERDEBUG] H2: current_c0=" << (int)current_c0 << " current_c1=" << (int)current_c1 << " amb_code=" << (int)amb_code << std::endl;
                            std::cout << "[SUPERDEBUG] H2: hap0_lane=" << (int)hap0 << " hap1_lane=" << (int)hap1 << std::endl;
                            std::cout << "[SUPERDEBUG] H2: sampled_h0=" << (int)h0 << " sampled_h1=" << (int)h1 << std::endl;
                        }
                    
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

                        // HYPOTHESIS 2 DEBUGGING
						if (is_superdebug_target) {
                            debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After projection in AMB");
                        }
                    // Persist sampled h0/h1 for this supersite anchor
                    setSupersiteClassPair(ss_idx_amb, h0, h1);
						
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
				if (!debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME && static_cast<int>(vabs) == debug::SUPERDEBUG_BP) {
					unsigned long dip_mask = Diplotypes[s];
					std::cout << "[SUPERDEBUG] BIAL: segment=" << s
					          << " dip_mask=0x" << std::hex << dip_mask << std::dec
					          << " DipSampled=" << static_cast<unsigned>(DipSampled[s])
					          << " HapLanes(" << static_cast<int>(hap0) << "," << static_cast<int>(hap1) << ")"
					          << " amb_code=" << static_cast<int>(Ambiguous[a])
					          << std::endl;
				}
			// TRACE: Log biallelic ambiguous interpretation for burn3 debugging
			if (supersite_trace_enabled() && vabs == 0) {
				unsigned char amb_code = Ambiguous[a];
				bool hap0_bit = HAP_GET(amb_code, hap0);
				bool hap1_bit = HAP_GET(amb_code, hap1);
				std::fprintf(stderr, "[MAKE_BIAL_AMB] sample=%s locus=%u\n",
				             name.c_str(), vabs);
				std::fprintf(stderr, "  sampled_lanes: hap0=%u hap1=%u\n",
				             (unsigned)hap0, (unsigned)hap1);
				std::fprintf(stderr, "  amb_code=0x%02x hap0_bit=%u hap1_bit=%u\n",
				             (unsigned)amb_code, (unsigned)hap0_bit, (unsigned)hap1_bit);
				std::fprintf(stderr, "  result: HAP0=%u HAP1=%u\n",
				             (unsigned)hap0_bit, (unsigned)hap1_bit);
			}
				HAP_GET(Ambiguous[a], hap0)?VAR_SET_HAP0(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP0(MOD2(vabs),Variants[DIV2(vabs)]);
				HAP_GET(Ambiguous[a], hap1)?VAR_SET_HAP1(MOD2(vabs),Variants[DIV2(vabs)]):VAR_CLR_HAP1(MOD2(vabs),Variants[DIV2(vabs)]);
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
	if (super_sites) {
		const size_t required = super_sites->size() * 2u;
		if (supersite_class_pairs.size() != required) {
			supersite_class_pairs.assign(required, SUPERSITE_CODE_MISSING);
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

		const bool trace_enabled = supersite_trace_enabled();
		const bool is_superdebug_target =
			(!debug::SUPERDEBUG_SAMPLENAME.empty() &&
			 this->name == debug::SUPERDEBUG_SAMPLENAME &&
			 static_cast<int>(ss.global_site_id) == debug::SUPERDEBUG_BP);

		if (is_superdebug_target) {
			debug::print_supersite_state(this, ss, *super_site_var_index, "projectSupersites:entry");
		}

        // Determine current sampled classes (h0/h1) for this anchor from stored pair
        // Note: getSupersiteClassPair returns the current epoch's h0/h1; c0/c1 remain immutable in emissions.
        uint8_t h0 = SUPERSITE_CODE_MISSING;
        uint8_t h1 = SUPERSITE_CODE_MISSING;
        getSupersiteClassPair(static_cast<int>(ss_idx), h0, h1);

        if (trace_enabled) {
            supersite_trace_log("[SupersiteProject] sample=%s ss_idx=%zu locus=%u h0=%u h1=%u\n",
                                name.c_str(),
                                ss_idx,
                                static_cast<unsigned>(ss.global_site_id),
                                static_cast<unsigned>(h0),
                                static_cast<unsigned>(h1));
        }
        if (is_superdebug_target) {
            std::cout << "[SUPERDEBUG] Sample=" << name
                      << " Pos=" << ss.global_site_id
                      << " Context='projectSupersites:anchor-classes'"
                      << " h0=" << static_cast<unsigned>(h0)
                      << " h1=" << static_cast<unsigned>(h1)
                      << std::endl;
        }

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

		if (is_superdebug_target) {
			debug::print_supersite_state(this, ss, *super_site_var_index, "projectSupersites:exit");
		}
	}
}

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

#ifndef _GENOTYPE_H
#define _GENOTYPE_H

#include <cstdint>

#include <utils/otools.h>

// Forward declaration to avoid circular dependency
struct SuperSite;

//Macros for packing/unpacking diplotypes
#define DIP_GET(dip,idx)	(((dip)>>(idx))&1UL)
#define DIP_SET(dip,idx)	((dip)|=(1UL<<(idx)))
#define DIP_HAP0(idx)		((idx)>>3)
#define DIP_HAP1(idx)		((idx)&7)

//Macros for packing/unpacking haplotypes
#define HAP_GET(hap, idx)	(((hap)>>(idx))&1U)
#define HAP_SET(hap, idx)	((hap)|=(1U<<(idx)))

//Macros for packing/unpacking variants
#define VAR_GET_HOM(e,v)	((((v)>>((e)<<2))&3)==0)
#define VAR_GET_MIS(e,v)	((((v)>>((e)<<2))&3)==1)
#define VAR_GET_HET(e,v)	((((v)>>((e)<<2))&3)==2)
#define VAR_GET_SCA(e,v)	((((v)>>((e)<<2))&3)==3)
//#define VAR_GET_AMB(e,v)	((((v)>>((e)<<2))&3)!=0)
#define VAR_GET_AMB(e,v)	((((v)>>((e)<<2))&3)>1)
#define VAR_SET_HOM(e,v)	((e)?((v)&=207):((v)&=252))
#define VAR_SET_MIS(e,v)	((v)|=(1<<((e)<<2)))
#define VAR_SET_HET(e,v)	((v)|=(2<<((e)<<2)))
#define VAR_SET_SCA(e,v)	((v)|=(3<<((e)<<2)))

#define VAR_GET_HAP0(e,v)	(((v)&(4<<((e)<<2)))!=0)
#define VAR_SET_HAP0(e,v)	((v)|=(4<<((e)<<2)))
#define VAR_CLR_HAP0(e,v)	((e)?((v)&=191):((v)&=251))
#define VAR_GET_HAP1(e,v)	(((v)&(8<<((e)<<2)))!=0)
#define VAR_SET_HAP1(e,v)	((v)|=(8<<((e)<<2)))
#define VAR_CLR_HAP1(e,v)	((e)?((v)&=127):((v)&=247))

// Macro to check if a variant is a supersite member
// Returns supersite index if variant is part of a supersite, -1 otherwise
// Usage: int ss_idx = VAR_GET_SS_IDX(locus_to_super_idx, vabs);
#define VAR_GET_SS_IDX(locus_to_super_idx, vabs) \
    ((locus_to_super_idx) ? (*locus_to_super_idx)[vabs] : -1)

// Macro to check if a variant is a supersite member (boolean result)
// Usage: if (VAR_IS_SUPERSITE(locus_to_super_idx, vabs)) { ... }
#define VAR_IS_SUPERSITE(locus_to_super_idx, vabs) \
    ((locus_to_super_idx) && ((*locus_to_super_idx)[vabs] >= 0))

// Macro to check if a variant is a supersite anchor (not a sibling split)
// Returns true if variant is the anchor position where HMM DP runs
// Usage: if (VAR_IS_SS_ANCHOR(super_sites, locus_to_super_idx, vabs)) { ... }
#define VAR_IS_SS_ANCHOR(super_sites, locus_to_super_idx, vabs) \
    ((super_sites) && (locus_to_super_idx) && \
     ((*locus_to_super_idx)[vabs] >= 0) && \
     ((int)(vabs) == (int)(*super_sites)[(*locus_to_super_idx)[vabs]].global_site_id))

// Macros for genotype::build() to handle supersite anchors vs siblings
// Supersite-related macros (now unused - kept for reference/documentation)
// Refactored build() now uses inline var_len pattern instead of these macros

#define MASK_INIT	0xFFFFFFFFFFFFFFFFUL
#define MASK_SCAF	0x00AA00AA00AA00AAUL
#define MASK_UNF0	0x55AA55AA55AA55AAUL
#define MASK_UNF1	0x3333CCCC3333CCCCUL
#define MASK_UNF2	0x0F0F0F0FF0F0F0F0UL

// Supersite flag bits stored per anchor (genotype::supersite_flags entries)
enum SupersiteFlagBits : uint8_t {
	SS_FLAG_HET      = 1u << 0,
	SS_FLAG_SCA      = 1u << 1,
	SS_FLAG_ALL_MIS  = 1u << 2
};


class genotype {
public:
	struct SuperSiteContext {
		int ss_idx = -1;
		bool is_member = false;
		bool is_anchor = false;
		bool has_het = false;
		bool has_sca = false;
		bool all_missing = false;
	};

	// INTERNAL DATA
	std::string name;
	unsigned int index;						// Index in containers
	unsigned int n_segments;				// Number of segments
	unsigned int n_variants;				// Number of variants	(to iterate over Variants)
	unsigned int n_ambiguous;				// Number of ambiguous variants
	unsigned int n_missing;					// Number of missing
	unsigned int n_transitions;				// Number of transitions
	unsigned int n_stored_transitionProbs;	// Number of transition probabilities stored in memory
	unsigned int n_storage_events;			// Number of storage having been done
	unsigned char curr_dipcodes [64];		// List of diplotypes in a given segment (buffer style variable)
	unsigned char curr_hapcodes [16];		// List of diplotypes in a given segment (buffer style variable)
	bool double_precision;					//If I get underflows using float, move to double
	bool haploid;							//Is this sample haploid?

	// VARIANT / HAPLOTYPE / DIPLOTYPE DATA
	std::vector < unsigned char > Variants;		// 0.5 byte per variant
	std::vector < unsigned char > Ambiguous;	// 1 byte per ambiguous variant
	std::vector < unsigned long > Diplotypes;	// 8 bytes per segment
	std::vector < unsigned short > Lengths;		// Raw variant span per segment (anchors count all splits)
	std::vector < unsigned short > Lengths_bio;	// Biological span per segment (anchors=1)

	//PHASE PROBS
	std::vector < bool > ProbMask;
	std::vector < float > ProbStored;
	std::vector < float > ProbMissing;

	// Supersite aggregate flags (per supersite index; cached in build())
	std::vector<uint8_t> supersite_flags;

	// SUPERSITE CONTEXT (Phase 3: multivariant imputation)
	const std::vector<SuperSite>* super_sites;
	const std::vector<int>* locus_to_super_idx;
	const std::vector<int>* super_site_var_index;
	const std::vector<float>* SC;  // CurrentSuperClassPosteriors from compute_job
	const std::vector<bool>* anchor_has_missing;
	const std::vector<uint32_t>* supersite_sc_offset;  // Thread-local SC offsets
	// Mutable per-epoch sampled classes (h0,h1) per supersite
	std::vector<uint8_t> supersite_class_pairs;
	// Immutable snapshot of site classes (c0,c1) per supersite, used for emissions
		std::vector<uint8_t> supersite_class_pairs_base;
		// Per-sample RNG for deterministic multithreaded runs
		random_number_generator sample_rng;

		//METHODS
		genotype(unsigned int);
		~genotype();
		void free();
		void seedRng(unsigned int base_seed);
		random_number_generator& rng();
		void make(std::vector < unsigned char > &, std::vector < float > &);
		void make(std::vector < unsigned char > &);
		void setSuperSiteContext(const std::vector<SuperSite>* _super_sites,
		                          const std::vector<int>* _locus_to_super_idx,
		                          const std::vector<int>* _super_site_var_index,
	                          const std::vector<float>* _SC,
	                          const std::vector<bool>* _anchor_has_missing,
	                          const std::vector<uint32_t>* _supersite_sc_offset = nullptr);
	void snapshotSupersiteClasses(const std::vector<SuperSite>& super_sites,
	                              const std::vector<int>& super_site_var_index);
	// Capture immutable c0/c1 snapshot for emissions
	void snapshotSupersiteBaseClasses(const std::vector<SuperSite>& super_sites,
	                                  const std::vector<int>& super_site_var_index);
	void setSupersiteClassPair(int ss_idx, uint8_t class0, uint8_t class1);
	void getSupersiteClassPair(int ss_idx, uint8_t& class0, uint8_t& class1) const;
	void getSupersiteBaseClassPair(int ss_idx, uint8_t& class0, uint8_t& class1) const;
	bool supersiteIsAmbiguous(int ss_idx) const;
	void build();
	void sample(std::vector < double > &, std::vector < float > &);
	void sampleForward(std::vector < double > &, std::vector < float > &);
	void sampleBackward(std::vector < double > &, std::vector < float > &);
	void mapMerges(std::vector < double > &, double , std::vector < bool > &);
	void performMerges(std::vector < double > &, std::vector < bool > &);
	void store(std::vector < double > &, std::vector < float > &);
	void solve();
	void projectSupersites();
	void scaffoldTrio(genotype *, genotype *, std::vector < unsigned int > &);
	void scaffoldDuoFather(genotype *, std::vector < unsigned int > &);
	void scaffoldDuoMother(genotype *, std::vector < unsigned int > &);
	uint32_t setHetsAsMissing();

	//INLINES
	unsigned int countDiplotypes(unsigned long);
	void makeDiplotypes(unsigned long);
	unsigned int countTransitions();
	bool isOrdered(unsigned long _dip);
	bool supersiteHasFlag(int ss_idx, uint8_t flag) const;
	SuperSiteContext getSuperSiteContext(unsigned int locus) const;
};

inline
bool genotype::isOrdered(unsigned long _dip) {
    std::fill(std::begin(curr_hapcodes), std::begin(curr_hapcodes)+16, 0);
	for (unsigned int d = 0, i = 0 ; d < 64 ; ++d) {
		if (DIP_GET(_dip, d)) {
			unsigned char hap0 = DIP_HAP0(d);
			unsigned char hap1 = DIP_HAP1(d);
			curr_hapcodes[hap0] = 1;
			curr_hapcodes[HAP_NUMBER + hap1] = 1;
		}
	}
	for (int h = 0 ; h < HAP_NUMBER ; h++) {
		if (curr_hapcodes[h] && curr_hapcodes[HAP_NUMBER+h]) {
			return false;
		}
	}
	return true;
}

inline
unsigned int genotype::countDiplotypes(unsigned long _dip) {
	unsigned int c = 0;
	for (unsigned long dip = _dip; dip; c++) dip &= dip - 1;
	return c;
}

inline
void genotype::makeDiplotypes(unsigned long _dip) {
	for (unsigned int d = 0, i = 0 ; d < 64 ; ++d) if (DIP_GET(_dip, d)) curr_dipcodes[i++] = d;
}

inline
unsigned int genotype::countTransitions() {
	unsigned int prev_dipcount = 1, c = 0;
	for (unsigned int s = 0 ; s < n_segments ; s++) {
		unsigned int curr_dipcount = countDiplotypes(Diplotypes[s]);
		c+= prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}
	return c;
}


#include <string>
#include <vector>

// Forward declarations
class genotype;
struct SuperSite;

namespace debug {
// Global debug settings
extern std::string SUPERDEBUG_SAMPLENAME;
extern int SUPERDEBUG_BP;

// Function to load settings from environment variables
void load_debug_settings();

// Function to print the state of a supersite for a genotype
void print_supersite_state(const genotype* G, const SuperSite& ss, const std::vector<int>& super_site_var_index, const std::string& context);
} // namespace debug

#endif

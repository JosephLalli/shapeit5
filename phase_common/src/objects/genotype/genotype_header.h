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

#include <cstddef>
#include <cstdint>

#include <utils/otools.h>

// Supersite constants (fallback if super_site_accessor.h not included yet)
#ifndef SUPERSITE_MAX_ALTS
#define SUPERSITE_MAX_ALTS 255
#endif

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
#define VAR_SET_MIS(e,v)	((v)=((v)&~(3<<((e)<<2)))|(1<<((e)<<2)))
#define VAR_SET_HET(e,v)	((v)=((v)&~(3<<((e)<<2)))|(2<<((e)<<2)))
#define VAR_SET_SCA(e,v)	((v)=((v)&~(3<<((e)<<2)))|(3<<((e)<<2)))

#define VAR_GET_HAP0(e,v)	(((v)&(4<<((e)<<2)))!=0)
#define VAR_SET_HAP0(e,v)	((v)|=(4<<((e)<<2)))
#define VAR_CLR_HAP0(e,v)	((e)?((v)&=191):((v)&=251))
#define VAR_GET_HAP1(e,v)	(((v)&(8<<((e)<<2)))!=0)
#define VAR_SET_HAP1(e,v)	((v)|=(8<<((e)<<2)))
#define VAR_CLR_HAP1(e,v)	((e)?((v)&=127):((v)&=247))

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

	// === SUPERSITE INDEXING SEMANTICS ===
	// When supersites are enabled, arrays are either:
	//   PER-VARIANT: Includes all biallelic splits (anchor + siblings)
	//   PER-BASEPAIR: Only biological positions (anchors count as 1, siblings excluded)

	// PER-VARIANT: Total variant count including all supersite siblings
	// Used to iterate over Variants[], which stores raw genotype data for every variant
	unsigned int n_variants;

	// PER-BASEPAIR: Count of ambiguous biological positions (het or sca)
	// Siblings are excluded because:
	// 1. Loop 3 in build() skips siblings: if (ctx.is_member && !ctx.is_anchor) continue;
	// 2. HMM excludes siblings from cursor advancement: data_amb = hmm_amb && !is_sibling;
	// Used to size Ambiguous[], which has 1 entry per ambiguous biological position
	unsigned int n_ambiguous;

	unsigned int n_missing;					// Number of missing
	unsigned int n_transitions;				// Number of transitions
	unsigned int n_stored_transitionProbs;	// Number of transition probabilities stored in memory
	unsigned int n_storage_events;			// Number of storage having been done
	unsigned char curr_dipcodes [64];		// List of diplotypes in a given segment (buffer style variable)
	unsigned char curr_hapcodes [16];		// List of diplotypes in a given segment (buffer style variable)
	bool double_precision;					//If I get underflows using float, move to double
	bool haploid;							//Is this sample haploid?

	// VARIANT / HAPLOTYPE / DIPLOTYPE DATA
	// === SUPERSITE INDEXING (see n_variants/n_ambiguous comments above) ===

	// PER-VARIANT: Raw genotype data, 0.5 byte per variant (including siblings)
	// Indexed by variant index [0, n_variants), packed 2 variants per byte
	std::vector < unsigned char > Variants;

	// PER-BASEPAIR: Ambiguous site masks, 1 byte per ambiguous biological position
	// Indexed by ambiguous index [0, n_ambiguous), siblings excluded
	// HMM cursor curr_abs_ambiguous indexes into this array
	std::vector < unsigned char > Ambiguous;

	std::vector < unsigned long > Diplotypes;	// 8 bytes per segment

	// PER-SEGMENT variant spans (both stored per segment, but measure different things):
	// Lengths[]: PER-VARIANT span - raw variant count including all sibling splits
	// Lengths_bio[]: PER-BASEPAIR span - biological position count (anchors=1, siblings=0)
	std::vector < unsigned short > Lengths;
	std::vector < unsigned short > Lengths_bio;

	//PHASE PROBS
	std::vector < bool > ProbMask;
	std::vector < float > ProbStored;
	std::vector < float > ProbMissing;
	std::vector < float > ProbSuperClass;   // Aggregated supersite class posteriors

	unsigned int sc_storage_events;         // Number of supersite SC aggregation events

	// Supersite aggregate flags (per supersite index; cached in build())
	std::vector<uint8_t> supersite_flags;

	// SUPERSITE CONTEXT (Phase 3: multivariant imputation)
	const std::vector<SuperSite>* super_sites;
	const std::vector<int>* locus_to_super_idx;
	const std::vector<int>* super_site_var_index;
	const uint8_t* supersite_panel_codes;
	size_t supersite_panel_codes_size;
	const std::vector<float>* SC;  // CurrentSuperClassPosteriors from compute_job
	const std::vector<bool>* anchor_has_missing;
	const std::vector<uint32_t>* supersite_sc_offset;  // Thread-local SC offsets
	// Mutable per-epoch sampled classes (h0,h1) per supersite
	std::vector<uint8_t> ss_phased_gts;
	// Immutable snapshot of site classes (c0,c1) per supersite, used for emissions
	std::vector<uint8_t> ss_observed_gts;
	// Missing mask per supersite: bit0=hap0 missing, bit1=hap1 missing
	std::vector<uint8_t> ss_missing_mask;
		// Per-sample RNG for deterministic multithreaded runs
		random_number_generator sample_rng;
		bool revert_buffer_fix; // Opt-in legacy sampler behavior

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
		void setSupersitePanelCodes(const uint8_t* _panel_codes, size_t _panel_codes_size);
	void snapshotSupersitePhasedGts(const std::vector<SuperSite>& super_sites,
	                                const std::vector<int>& super_site_var_index);
	// Capture immutable c0/c1 snapshot for emissions
	void snapshotSupersiteObservedGts(const std::vector<SuperSite>& super_sites,
	                                  const std::vector<int>& super_site_var_index);
	void setSupersiteObservedGt(int ss_idx, uint8_t c0, uint8_t c1, uint8_t missing_mask);
	void setSupersiteMissingMask(int ss_idx, uint8_t mask);
	uint8_t getSupersiteMissingMask(int ss_idx) const;
	bool supersiteIsMissing(int ss_idx) const;
	void setRevertBufferFix(bool value) { revert_buffer_fix = value; }
	void setSupersitePhasedGt(int ss_idx, uint8_t h0, uint8_t h1);
	void getSupersitePhasedGt(int ss_idx, uint8_t& h0, uint8_t& h1) const;
	void getSupersiteObservedGt(int ss_idx, uint8_t& c0, uint8_t& c1) const;
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
	size_t supersite_class_index(int ss_idx, int lane, int cls) const;

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

inline size_t genotype::supersite_class_index(int ss_idx, int lane, int cls) const {
	constexpr int MAX_CLASSES = SUPERSITE_MAX_ALTS + 1;
	return static_cast<size_t>(ss_idx) * HAP_NUMBER * MAX_CLASSES +
	       static_cast<size_t>(lane) * MAX_CLASSES +
	       static_cast<size_t>(cls);
}


#endif

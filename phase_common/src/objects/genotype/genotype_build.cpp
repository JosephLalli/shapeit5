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

#define MASK_INIT	0xFFFFFFFFFFFFFFFFUL
#define MASK_SCAF	0x00AA00AA00AA00AAUL
#define MASK_UNF0	0x55AA55AA55AA55AAUL
#define MASK_UNF1	0x3333CCCC3333CCCCUL
#define MASK_UNF2	0x0F0F0F0FF0F0F0F0UL

#include <objects/genotype/genotype_header.h>
#include <models/super_site_accessor.h>

using namespace std;

uint32_t genotype::setHetsAsMissing() {
	uint32_t nreset = 0;
	for (uint32_t v = 0 ; v < n_variants ; v++) {
		if (VAR_GET_AMB(MOD2(v), Variants[DIV2(v)])) {
			VAR_SET_MIS(MOD2(v), Variants[DIV2(v)]);
			nreset ++;
		}
	}
	return nreset;
}

void genotype::build() {
	//1. Count number of segments
	unsigned n_rel_unf = 0, n_rel_var = 0, n_rel_sca = 0, n_abs_seg = 0, n_abs_amb = 0, n_rel_amb = 0, n_abs_mis = 0;
	unsigned var_len = 0;

	const bool has_supersites = super_sites && locus_to_super_idx && super_site_var_index;
	if (has_supersites) {
		supersite_flags.assign(super_sites->size(), 0u);
		for (size_t ss_idx = 0; ss_idx < super_sites->size(); ++ss_idx) {
			const SuperSite &ss = (*super_sites)[ss_idx];
			bool any_het = false;
			bool any_sca = false;
			bool all_missing = (ss.var_count > 0);
			for (uint32_t offset = 0; offset < ss.var_count; ++offset) {
				int v_idx = (*super_site_var_index)[ss.var_start + offset];
				unsigned char v_code = Variants[DIV2(v_idx)];
				bool is_het = VAR_GET_HET(MOD2(v_idx), v_code);
				any_het |= is_het;
				any_sca |= VAR_GET_SCA(MOD2(v_idx), v_code);
				if (!VAR_GET_MIS(MOD2(v_idx), v_code)) {
					all_missing = false;
				}
				if (std::getenv("SHAPEIT5_TEST_TRACE")) {
					std::fprintf(stderr, "[build] ss_idx=%zu offset=%u v_idx=%d byte=0x%02x is_het=%d\n",
					            ss_idx, offset, v_idx, v_code, (int)is_het);
				}
			}
			uint8_t flags = 0u;
			if (any_het) flags |= SS_FLAG_HET;
			if (any_sca) flags |= SS_FLAG_SCA;
			if (all_missing) flags |= SS_FLAG_ALL_MIS;
			supersite_flags[ss_idx] = flags;
			if (std::getenv("SHAPEIT5_TEST_TRACE")) {
				std::fprintf(stderr, "[build] ss_idx=%zu any_het=%d flags=0x%02x\n",
				            ss_idx, (int)any_het, flags);
			}
		}
	} else {
		supersite_flags.clear();
	}
	// First loop: Count segments while preventing empty segments from supersite siblings
	for (unsigned int v = 0 ; v < n_variants ;) {
		unsigned char var_code = Variants[DIV2(v)];
		bool f_sca = VAR_GET_SCA(MOD2(v), var_code);
		bool f_het = VAR_GET_HET(MOD2(v), var_code);
		bool f_mis = VAR_GET_MIS(MOD2(v), var_code);
		SuperSiteContext ctx = getSuperSiteContext(v);
		int ss_idx = ctx.ss_idx;
		bool is_anchor = ctx.is_member && ctx.is_anchor;

		// Determine how many variants to process (1 for normal, var_count for supersite anchor)
		if (is_anchor) {
			var_len = (*super_sites)[ss_idx].var_count;
			f_sca = ctx.has_sca;
			f_het = ctx.has_het;
			f_mis = ctx.all_missing;
		} else {
			var_len = 1;
		}
		const unsigned eff_len = is_anchor ? 1u : var_len; // treat supersite anchor as a single biological site

		bool is_sibling = ctx.is_member && !ctx.is_anchor; // Identify sibling variants
		bool f_het_for_boundary = f_het;
		bool f_sca_for_boundary = f_sca;
		if (is_sibling) {
			// CRITICAL FIX: Exclude siblings from boundary calculation to prevent empty segments
			// Siblings represent alternate alleles at the same genomic position as their anchor.
			// They do not add to n_rel_var independently (anchor processes all var_count variants).
			// If siblings contributed to predicted_unfold, a boundary could be triggered with n_rel_var=0,
			// creating an empty segment that causes AlphaSumSum underflow in backward pass TRANS_HAP().
			f_het_for_boundary = false; // Treat sibling as non-het for boundary decision
			f_sca_for_boundary = false; // Treat sibling as non-sca for boundary decision
		}

		unsigned int predicted_unfold = n_rel_unf + f_het_for_boundary + (n_rel_sca||f_sca_for_boundary);
		if (predicted_unfold == 4 || (n_rel_var >= (std::numeric_limits<unsigned short>::max() - eff_len + 1)) || (n_rel_amb == MAX_AMB)) {
			// Segment boundary
			n_rel_unf = 0;
			n_rel_sca = 0;
			n_rel_var = 0;
			n_rel_amb = 0;
			n_abs_seg++;
		} else {
			n_rel_unf += f_het;
			n_rel_sca += f_sca;
			n_abs_amb += (f_het||f_sca);
			n_rel_amb += (f_het||f_sca);
			n_abs_mis += f_mis;
			n_rel_var += eff_len;
			v += var_len;
		}
	}
	n_segments = n_abs_seg + 1;
	n_ambiguous = n_abs_amb;
	n_missing = n_abs_mis;

	//2. Build Segments (same logic as loop 1, with same sibling boundary fix)
	n_rel_unf = 0; n_rel_var = 0; n_rel_sca = 0; n_abs_seg = 0; n_abs_amb = 0; n_rel_amb = 0; n_abs_mis = 0;
	unsigned int n_rel_var_eff = 0; // effective biological length (anchors count as 1)
	Lengths = vector < unsigned short > (n_segments, 0U);      // raw variant span
	Lengths_bio = vector < unsigned short > (n_segments, 0U);  // biological span

	for (unsigned int v = 0 ; v < n_variants ;) {
		unsigned char var_code = Variants[DIV2(v)];
		bool f_sca = VAR_GET_SCA(MOD2(v), var_code);
		bool f_het = VAR_GET_HET(MOD2(v), var_code);
		bool f_mis = VAR_GET_MIS(MOD2(v), var_code);
		SuperSiteContext ctx = getSuperSiteContext(v);
		int ss_idx = ctx.ss_idx;
		bool is_anchor = ctx.is_member && ctx.is_anchor;
		
		// Determine how many variants to process (1 for normal, var_count for supersite anchor)
		if (is_anchor) {
			var_len = (*super_sites)[ss_idx].var_count;
			f_sca = ctx.has_sca;
			f_het = ctx.has_het;
			f_mis = ctx.all_missing;
		} else {
			var_len = 1;
		}
		const unsigned eff_len = is_anchor ? 1u : var_len;

		bool is_sibling = ctx.is_member && !ctx.is_anchor; // Identify sibling variants
		bool f_het_for_boundary = f_het;
		bool f_sca_for_boundary = f_sca;
		if (is_sibling) {
			// CRITICAL FIX: Exclude siblings from boundary calculation (same as loop 1)
			// See detailed comment in loop 1 above for full explanation.
			f_het_for_boundary = false; // Treat sibling as non-het for boundary decision
			f_sca_for_boundary = false; // Treat sibling as non-sca for boundary decision
		}

		unsigned int predicted_unfold = n_rel_unf + f_het_for_boundary + (n_rel_sca||f_sca_for_boundary);
		if (predicted_unfold == 4 || (n_rel_var_eff >= (std::numeric_limits<unsigned short>::max() - eff_len + 1)) || (n_rel_amb == MAX_AMB)) {
			// Segment boundary
#if !defined(__OPTIMIZE__)
			if (n_rel_var > std::numeric_limits<unsigned short>::max()) {
				std::fprintf(stderr, "[SEGMENT_OVERFLOW] Sample=%s: raw segment length=%u exceeds uint16_t at v=%u\n",
				             name.c_str(), n_rel_var, v);
				std::fflush(stderr);
				std::abort();
			}
#endif
			Lengths[n_abs_seg] = n_rel_var;            // raw span (includes siblings)
			Lengths_bio[n_abs_seg] = n_rel_var_eff;    // biological span (anchors=1)
			// ASSERTION: Detect empty segment creation that would cause underflow
			if (n_rel_var == 0 && n_abs_seg > 0) {
				std::fprintf(stderr, "[EMPTY_SEGMENT_ERROR] Sample=%s: Created empty segment %u at variant v=%u\n",
				             name.c_str(), n_abs_seg, v);
				std::fprintf(stderr, "  Boundary triggered by: predicted_unfold=%u n_rel_var=%u n_rel_amb=%u\n",
				             predicted_unfold, n_rel_var, n_rel_amb);
				std::fprintf(stderr, "  This causes AlphaSumSum underflow in backward pass\n");
				std::fprintf(stderr, "  Previous segment: Lengths[%u]=%u\n",
				             n_abs_seg > 0 ? n_abs_seg - 1 : 0,
				             n_abs_seg > 0 ? Lengths[n_abs_seg - 1] : 0);
				assert(false && "Empty segment created - will cause underflow in TRANS_HAP()");
			}
			n_rel_unf = 0;
			n_rel_sca = 0;
			n_rel_var = 0;
			n_rel_var_eff = 0;
			n_rel_amb = 0;
			n_abs_seg++;
		} else {
			n_rel_unf += f_het;
			n_rel_sca += f_sca;
			n_abs_amb += (f_het||f_sca);
			n_rel_amb += (f_het||f_sca);
			n_abs_mis += f_mis;
			n_rel_var += var_len;       // real span for storage and indexing
			n_rel_var_eff += eff_len;   // effective span for boundary decisions
			v += var_len;
		}
	}
	Lengths[n_abs_seg] = n_rel_var;
	Lengths_bio[n_abs_seg] = n_rel_var_eff;
#if !defined(__OPTIMIZE__)
	if (n_rel_var > std::numeric_limits<unsigned short>::max()) {
		std::fprintf(stderr, "[SEGMENT_OVERFLOW] Sample=%s: final raw segment length=%u exceeds uint16_t (n_variants=%u)\n",
		             name.c_str(), n_rel_var, n_variants);
		std::fflush(stderr);
		std::abort();
	}
#endif
	// ASSERTION: Detect final empty segment creation
	if (n_rel_var == 0 && n_abs_seg > 0) {
		std::fprintf(stderr, "[EMPTY_SEGMENT_ERROR] Sample=%s: Created final empty segment %u (total n_variants=%u)\n",
		             name.c_str(), n_abs_seg, n_variants);
		std::fprintf(stderr, "  Loop terminated with n_rel_var=0 after last boundary\n");
		std::fprintf(stderr, "  This causes AlphaSumSum underflow in backward pass\n");
		std::fprintf(stderr, "  Previous segment: Lengths[%u]=%u\n",
		             n_abs_seg > 0 ? n_abs_seg - 1 : 0,
		             n_abs_seg > 0 ? Lengths[n_abs_seg - 1] : 0);
		assert(false && "Final empty segment created - will cause underflow in TRANS_HAP()");
	}

	//3. Build Ambiguous
	Ambiguous = vector < unsigned char >(n_ambiguous, 0U);
	vector < unsigned char > orderedSegments = vector < unsigned char >(n_segments, 0);
	for (unsigned int s = 0, a0 = 0, a1 = 0, a2 = 0, vabs = 0 ; s < n_segments ; s ++) {
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) {
			unsigned int v_idx = vabs + vrel;
			unsigned char var_code = Variants[DIV2(v_idx)];
			bool f_sca = VAR_GET_SCA(MOD2(v_idx), var_code);
			bool f_het = VAR_GET_HET(MOD2(v_idx), var_code);
			SuperSiteContext ctx = getSuperSiteContext(v_idx);
			if (ctx.is_member && !ctx.is_anchor) continue;
			if (ctx.is_anchor) {
				f_sca = f_sca || ctx.has_sca;
				f_het = ctx.has_het;
			}
			
			if (f_sca) {
				for (unsigned int h = 0 ; h < HAP_NUMBER ; h ++) {
					bool allele = (h%2)?VAR_GET_HAP1(MOD2(v_idx), var_code):VAR_GET_HAP0(MOD2(v_idx), var_code);
					if (allele) HAP_SET(Ambiguous[a0], h);
				}
				orderedSegments[s] = 1;
			}
			a0 += (f_sca||f_het);
		}
		unsigned int n_unf = orderedSegments[s];
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) {
			unsigned int v_idx = vabs + vrel;
			unsigned char var_code = Variants[DIV2(v_idx)];
			bool f_sca = VAR_GET_SCA(MOD2(v_idx), var_code);
			bool f_het = VAR_GET_HET(MOD2(v_idx), var_code);
			SuperSiteContext ctx = getSuperSiteContext(v_idx);
			if (ctx.is_member && !ctx.is_anchor) continue;
			if (ctx.is_anchor) {
				f_sca = f_sca || ctx.has_sca;
				f_het = ctx.has_het;
			}
			
			if (f_het) {
				for (unsigned int h = 0 ; h < HAP_NUMBER ; h ++) {
					bool allele = ((h>>n_unf)%2);
					if (allele) HAP_SET(Ambiguous[a1], h);
				}
				n_unf++;
			}
			a1 += (f_sca||f_het);
		}
		vabs += Lengths[s];
	}

	//4. Build Diplotypes
	Diplotypes = vector < unsigned long > (n_segments);
	for (unsigned int s = 0, vabs = 0, a = 0 ; s < n_segments ; s ++) {
		unsigned int n_unf = orderedSegments[s];
		Diplotypes[s] = n_unf ? MASK_SCAF : MASK_INIT;
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) {
			unsigned int v_idx = vabs + vrel;
			unsigned char var_code = Variants[DIV2(v_idx)];
			bool f_het = VAR_GET_HET(MOD2(v_idx), var_code);
			SuperSiteContext ctx = getSuperSiteContext(v_idx);
			if (ctx.is_member && !ctx.is_anchor) continue;
			if (ctx.is_anchor) {
				f_het = ctx.has_het;
			}
			
			if (f_het) {
				switch (n_unf) {
				case 0: Diplotypes[s] &= MASK_UNF0; break;
				case 1: Diplotypes[s] &= MASK_UNF1; break;
				case 2: Diplotypes[s] &= MASK_UNF2; break;
				}
			}
			n_unf += f_het;
		}
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) {
			unsigned int v_idx = vabs + vrel;
			unsigned char var_code = Variants[DIV2(v_idx)];
			SuperSiteContext ctx = getSuperSiteContext(v_idx);
			if (ctx.is_member && !ctx.is_anchor) continue;
			bool is_amb = VAR_GET_AMB(MOD2(v_idx), var_code);
			if (ctx.is_anchor) {
				is_amb = ctx.has_sca || ctx.has_het;
			}
			a += is_amb;
		}

		vabs += Lengths[s];
	}

	//5. Count transitions
	n_transitions = countTransitions();

	//6. ASSERTION: Final validation that no empty segments exist
	for (unsigned int s = 0; s < n_segments; s++) {
		if (Lengths[s] == 0) {
			std::fprintf(stderr, "[EMPTY_SEGMENT_ERROR] Sample=%s: CRITICAL - Empty segment detected at index %u/%u\n",
			             name.c_str(), s, n_segments);
			std::fprintf(stderr, "  n_variants=%u n_segments=%u n_ambiguous=%u\n",
			             n_variants, n_segments, n_ambiguous);
			std::fprintf(stderr, "  Segment lengths:");
			for (unsigned int i = 0; i < n_segments && i < 20; i++) {
				std::fprintf(stderr, " [%u]=%u", i, Lengths[i]);
			}
			std::fprintf(stderr, "\n");
			std::fprintf(stderr, "  This WILL cause underflow in backward pass TRANS_HAP()\n");
			assert(false && "Empty segment in built genotype - BUG in genotype::build()");
		}
	}
}

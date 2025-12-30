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
	const char* trace_build_loops = std::getenv("SHAPEIT5_TRACE_BUILD_LOOPS");
	const bool trace_build_loops_on = trace_build_loops && trace_build_loops[0] != '\0' && trace_build_loops[0] != '0';

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
			}
			uint8_t flags = 0u;
			if (any_het) flags |= SS_FLAG_HET;
			if (any_sca) flags |= SS_FLAG_SCA;
			if (all_missing) flags |= SS_FLAG_ALL_MIS;
			supersite_flags[ss_idx] = flags;
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
		const bool is_sibling = ctx.is_member && !ctx.is_anchor; // Identify sibling variants
		// Biological length counts anchors as 1 and siblings as 0.
		const unsigned eff_len = is_sibling ? 0u : (is_anchor ? 1u : var_len);

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
		if (trace_build_loops_on && name == "HG00107" && v >= 0 && v <= 83810) {
			std::fprintf(stderr, "[LOOP1_BOUNDARY] v=%u n_abs_amb=%u predicted_unfold=%u f_het=%d is_sibling=%d is_member=%d is_anchor=%d ss_idx=%d\n",
							v, n_abs_amb, predicted_unfold, f_het, is_sibling, ctx.is_member, ctx.is_anchor, ctx.ss_idx);
		}
		if (predicted_unfold == 4 || (n_rel_var >= (std::numeric_limits<unsigned short>::max() - eff_len + 1)) || (n_rel_amb == MAX_AMB)) {
			// Segment boundary
			// DEBUG: Log boundary trigger for failing sample - check v=83803
			n_rel_unf = 0;
			n_rel_sca = 0;
			n_rel_var = 0;
			n_rel_amb = 0;
			n_abs_seg++;
		} else {
			// n_abs_amb counts PER-BASEPAIR ambiguous sites (used to size Ambiguous[])
			// Siblings are excluded because:
			// 1. Loop 3 skips siblings: if (ctx.is_member && !ctx.is_anchor) continue;
			// 2. HMM excludes siblings from cursor advancement: data_amb = hmm_amb && !is_sibling;
			// Therefore we use _for_boundary versions which are zeroed for siblings
			const bool is_sibling = (ctx.is_member && !ctx.is_anchor);
			// CRITICAL FIX: For anchors with scaffold status, f_sca_eff is false.
			// This matches Loop 3's advancement logic at lines 291 and 331: a0/a1 += (f_sca_eff||f_het)
			// Without this, Loop 1 counts more ambiguous sites than Loop 3 fills, causing heap overflow.
			const bool f_sca_eff = f_sca && !(ctx.is_anchor && ctx.has_sca);
			const bool amb_for_count = is_sibling ? false : (f_het || f_sca_eff);
			// DEBUG: Check v=83803 specifically
			if (trace_build_loops_on && name == "HG00107" && v >= 83800 && v <= 83810) {
				std::fprintf(stderr, "[LOOP1_ELSE] v=%u n_abs_amb=%u f_het=%d amb_for_count=%d is_sibling=%d is_member=%d is_anchor=%d ss_idx=%d\n",
				             v, n_abs_amb, f_het, amb_for_count, is_sibling, ctx.is_member, ctx.is_anchor, ctx.ss_idx);
			}
			n_rel_unf += f_het_for_boundary;
			n_rel_sca += f_sca_for_boundary;
			n_abs_amb += amb_for_count;
			n_rel_amb += amb_for_count;
			n_abs_mis += f_mis;
			n_rel_var += eff_len;
			v += var_len;
		}
	}
	n_segments = n_abs_seg + 1;
	n_ambiguous = n_abs_amb;
	n_missing = n_abs_mis;

	// DEBUG: Log final state for the failing sample
	if (trace_build_loops_on && name == "HG00107") {
		std::fprintf(stderr, "[LOOP1_END] sample=%s n_variants=%u n_ambiguous=%u n_segments=%u\n",
		             name.c_str(), n_variants, n_ambiguous, n_segments);
	}

	//2. Build Segments (same logic as loop 1, with same sibling boundary fix)
	n_rel_unf = 0; n_rel_var = 0; n_rel_sca = 0; n_abs_seg = 0; n_abs_amb = 0; n_rel_amb = 0; n_abs_mis = 0;
	unsigned int n_rel_var_eff = 0; // effective biological length (anchors count as 1, siblings=0)
	unsigned int seg_start_v = 0;
	Lengths = vector < unsigned short > (n_segments, 0U);      // raw variant span
	Lengths_bio = vector < unsigned short > (n_segments, 0U);  // biological span
	const char* trace_lenbio_detail = std::getenv("SHAPEIT5_TRACE_LENBIO_DETAIL");
	const char* trace_lenbio_sample = std::getenv("SHAPEIT5_TRACE_LENBIO_SAMPLE");
	const bool trace_lenbio_detail_on = trace_lenbio_detail && trace_lenbio_detail[0] != '\0' && trace_lenbio_detail[0] != '0';
	const bool trace_lenbio_sample_ok = (!trace_lenbio_sample || trace_lenbio_sample[0] == '\0' || name == trace_lenbio_sample);
	std::vector<uint8_t> lenbio_detail_traced(trace_lenbio_detail_on ? n_segments : 0, 0);
	const char* trace_lenbio_dump = std::getenv("SHAPEIT5_TRACE_LENBIO_DUMP");
	const char* trace_lenbio_dump_all = std::getenv("SHAPEIT5_TRACE_LENBIO_DUMP_ALL");
	const char* trace_lenbio_seg = std::getenv("SHAPEIT5_TRACE_LENBIO_SEG");
	const char* trace_lenbio_trigger = std::getenv("SHAPEIT5_TRACE_LENBIO_TRIGGER");
	const char* trace_lenbio_ring = std::getenv("SHAPEIT5_TRACE_LENBIO_RING");
	const char* trace_lenbio_maxlines = std::getenv("SHAPEIT5_TRACE_LENBIO_MAXLINES");
	const bool trace_lenbio_dump_on = trace_lenbio_dump && trace_lenbio_dump[0] != '\0' && trace_lenbio_dump[0] != '0';
	const bool trace_lenbio_dump_all_on = trace_lenbio_dump_all && trace_lenbio_dump_all[0] != '\0' && trace_lenbio_dump_all[0] != '0';
	const unsigned int trace_lenbio_seg_idx = (trace_lenbio_seg && trace_lenbio_seg[0] != '\0') ? static_cast<unsigned int>(std::strtoul(trace_lenbio_seg, nullptr, 10)) : 0U;
	const bool trace_lenbio_seg_filter_on = (trace_lenbio_seg && trace_lenbio_seg[0] != '\0');
	const unsigned int trace_lenbio_trigger_val = (trace_lenbio_trigger && trace_lenbio_trigger[0] != '\0') ? static_cast<unsigned int>(std::strtoul(trace_lenbio_trigger, nullptr, 10)) : 0U;
	const unsigned int trace_lenbio_ring_size = (trace_lenbio_ring && trace_lenbio_ring[0] != '\0') ? static_cast<unsigned int>(std::strtoul(trace_lenbio_ring, nullptr, 10)) : 256U;
	const size_t trace_lenbio_maxlines_val = (trace_lenbio_maxlines && trace_lenbio_maxlines[0] != '\0') ? static_cast<size_t>(std::strtoull(trace_lenbio_maxlines, nullptr, 10)) : 200000U;
	const bool trace_lenbio_active = trace_lenbio_dump_on && trace_lenbio_sample_ok;
	size_t trace_lenbio_lines = 0;
	bool trace_lenbio_streaming = false;
	struct LenbioTraceEvent {
		unsigned int seg_idx;
		unsigned int v;
		unsigned int seg_start_v;
		unsigned int n_rel_var;
		unsigned int n_rel_var_eff;
		unsigned int n_rel_amb;
		unsigned int n_rel_unf;
		unsigned int n_rel_sca;
		unsigned int var_len;
		unsigned int eff_len;
		unsigned int predicted_unfold;
		int ss_idx;
		bool is_member;
		bool is_anchor;
		bool f_het;
		bool f_sca;
		bool f_mis;
	};
	std::vector<LenbioTraceEvent> lenbio_ring_buf;
	size_t lenbio_ring_pos = 0;
	bool lenbio_ring_full = false;
	if (trace_lenbio_active && trace_lenbio_ring_size > 0) {
		lenbio_ring_buf.resize(trace_lenbio_ring_size);
	}
	auto trace_lenbio_seg_ok = [&](unsigned int seg_idx) -> bool {
		return !trace_lenbio_seg_filter_on || seg_idx == trace_lenbio_seg_idx;
	};
	auto trace_lenbio_emit = [&](const LenbioTraceEvent& ev, const char* tag) {
		if (!trace_lenbio_active) return;
		if (!trace_lenbio_seg_ok(ev.seg_idx)) return;
		if (trace_lenbio_lines >= trace_lenbio_maxlines_val) return;
		std::fprintf(stderr,
		             "[LENBIO_DUMP][build-%s] sample=%s seg=%u v=%u seg_start=%u len=%u len_bio=%u n_rel_amb=%u n_rel_unf=%u n_rel_sca=%u "
		             "pred_unf=%u var_len=%u eff_len=%u ss_idx=%d member=%d anchor=%d f_het=%d f_sca=%d f_mis=%d\n",
		             tag, name.c_str(), ev.seg_idx, ev.v, ev.seg_start_v, ev.n_rel_var, ev.n_rel_var_eff,
		             ev.n_rel_amb, ev.n_rel_unf, ev.n_rel_sca, ev.predicted_unfold,
		             ev.var_len, ev.eff_len, ev.ss_idx, ev.is_member ? 1 : 0, ev.is_anchor ? 1 : 0,
		             ev.f_het ? 1 : 0, ev.f_sca ? 1 : 0, ev.f_mis ? 1 : 0);
		trace_lenbio_lines++;
	};
	auto trace_lenbio_ring_add = [&](const LenbioTraceEvent& ev) {
		if (!trace_lenbio_active || lenbio_ring_buf.empty()) return;
		if (!trace_lenbio_seg_ok(ev.seg_idx)) return;
		lenbio_ring_buf[lenbio_ring_pos] = ev;
		lenbio_ring_pos = (lenbio_ring_pos + 1) % lenbio_ring_buf.size();
		if (lenbio_ring_pos == 0) lenbio_ring_full = true;
	};
	auto trace_lenbio_ring_dump = [&](const char* tag) {
		if (!trace_lenbio_active || lenbio_ring_buf.empty()) return;
		const size_t n = lenbio_ring_full ? lenbio_ring_buf.size() : lenbio_ring_pos;
		const size_t start = lenbio_ring_full ? lenbio_ring_pos : 0;
		for (size_t i = 0; i < n; ++i) {
			const size_t idx = (start + i) % lenbio_ring_buf.size();
			trace_lenbio_emit(lenbio_ring_buf[idx], tag);
		}
	};

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
		const bool is_sibling = ctx.is_member && !ctx.is_anchor; // Identify sibling variants
		// Biological length counts anchors as 1 and siblings as 0.
		const unsigned eff_len = is_sibling ? 0u : (is_anchor ? 1u : var_len);

		bool f_het_for_boundary = f_het;
		bool f_sca_for_boundary = f_sca;
		if (is_sibling) {
			// CRITICAL FIX: Exclude siblings from boundary calculation (same as loop 1)
			// See detailed comment in loop 1 above for full explanation.
			f_het_for_boundary = false; // Treat sibling as non-het for boundary decision
			f_sca_for_boundary = false; // Treat sibling as non-sca for boundary decision
		}

		unsigned int predicted_unfold = n_rel_unf + f_het_for_boundary + (n_rel_sca||f_sca_for_boundary);
		if (trace_lenbio_active) {
			LenbioTraceEvent ev;
			ev.seg_idx = n_abs_seg;
			ev.v = v;
			ev.seg_start_v = seg_start_v;
			ev.n_rel_var = n_rel_var;
			ev.n_rel_var_eff = n_rel_var_eff;
			ev.n_rel_amb = n_rel_amb;
			ev.n_rel_unf = n_rel_unf;
			ev.n_rel_sca = n_rel_sca;
			ev.var_len = var_len;
			ev.eff_len = eff_len;
			ev.predicted_unfold = predicted_unfold;
			ev.ss_idx = ss_idx;
			ev.is_member = ctx.is_member;
			ev.is_anchor = ctx.is_anchor;
			ev.f_het = f_het;
			ev.f_sca = f_sca;
			ev.f_mis = f_mis;
			trace_lenbio_ring_add(ev);
			if (trace_lenbio_dump_all_on || trace_lenbio_streaming) {
				trace_lenbio_emit(ev, "step");
			} else if (trace_lenbio_trigger_val > 0 && n_rel_var_eff >= trace_lenbio_trigger_val) {
				trace_lenbio_streaming = true;
				trace_lenbio_ring_dump("ring");
				trace_lenbio_emit(ev, "trigger");
			}
		}
		if (predicted_unfold == 4 ||
		    (n_rel_var >= (std::numeric_limits<unsigned short>::max() - var_len + 1)) ||
		    (n_rel_var_eff >= (std::numeric_limits<unsigned short>::max() - eff_len + 1)) ||
		    (n_rel_amb == MAX_AMB)) {
			if (trace_lenbio_active) {
				LenbioTraceEvent ev;
				ev.seg_idx = n_abs_seg;
				ev.v = v;
				ev.seg_start_v = seg_start_v;
				ev.n_rel_var = n_rel_var;
				ev.n_rel_var_eff = n_rel_var_eff;
				ev.n_rel_amb = n_rel_amb;
				ev.n_rel_unf = n_rel_unf;
				ev.n_rel_sca = n_rel_sca;
				ev.var_len = var_len;
				ev.eff_len = eff_len;
				ev.predicted_unfold = predicted_unfold;
				ev.ss_idx = ss_idx;
				ev.is_member = ctx.is_member;
				ev.is_anchor = ctx.is_anchor;
				ev.f_het = f_het;
				ev.f_sca = f_sca;
				ev.f_mis = f_mis;
				trace_lenbio_ring_dump("ring");
				trace_lenbio_emit(ev, "boundary");
			}
			// One-off debug guard to explain len_bio overflow on segment closure; remove after root cause is found.
			if (trace_lenbio_detail_on && trace_lenbio_sample_ok &&
			    !lenbio_detail_traced.empty() && !lenbio_detail_traced[n_abs_seg]) {
				const bool hit_unf = (predicted_unfold == 4);
				const bool hit_len_raw = (n_rel_var >= (std::numeric_limits<unsigned short>::max() - var_len + 1));
				const bool hit_len_bio = (n_rel_var_eff >= (std::numeric_limits<unsigned short>::max() - eff_len + 1));
				const bool hit_amb = (n_rel_amb == MAX_AMB);
				unsigned int non_sib = 0, sib = 0, anchor = 0, bial = 0;
				for (unsigned int vrel = 0; vrel < n_rel_var; ++vrel) {
					unsigned int locus = seg_start_v + vrel;
					SuperSiteContext ctx2 = getSuperSiteContext(locus);
					if (ctx2.is_member) {
						if (ctx2.is_anchor) {
							anchor++;
							non_sib++;
						} else {
							sib++;
						}
					} else {
						bial++;
						non_sib++;
					}
				}
				const unsigned int seg_end_v = (v > 0) ? (v - 1) : 0;
				std::fprintf(stderr,
				             "[LENBIO_DETAIL][build-boundary] sample=%s seg=%u start_v=%u end_v=%u next_v=%u len=%u len_bio=%u "
				             "non_sib=%u sib=%u anchor=%u bial=%u "
				             "hit_unf=%u hit_len_raw=%u hit_len_bio=%u hit_amb=%u "
				             "pred_unf=%u n_rel_unf=%u n_rel_sca=%u n_rel_amb=%u "
				             "var_len=%u eff_len=%u ss_idx=%d member=%d anchor=%d f_het=%d f_sca=%d f_mis=%d\n",
				             name.c_str(), n_abs_seg, seg_start_v, seg_end_v, v, n_rel_var, n_rel_var_eff,
				             non_sib, sib, anchor, bial,
				             hit_unf ? 1u : 0u, hit_len_raw ? 1u : 0u, hit_len_bio ? 1u : 0u, hit_amb ? 1u : 0u,
				             predicted_unfold, n_rel_unf, n_rel_sca, n_rel_amb,
				             var_len, eff_len, ss_idx, ctx.is_member ? 1 : 0, ctx.is_anchor ? 1 : 0,
				             f_het ? 1 : 0, f_sca ? 1 : 0, f_mis ? 1 : 0);
				lenbio_detail_traced[n_abs_seg] = 1;
			}
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
			Lengths_bio[n_abs_seg] = n_rel_var_eff;    // biological span (anchors=1, siblings=0)
			// TEMP DEBUG: trace segments that hit the uint16 cap or look inconsistent.
			const char* trace_lenbio = std::getenv("SHAPEIT5_TRACE_LENBIO");
			if (trace_lenbio && trace_lenbio[0] != '\0' && trace_lenbio[0] != '0') {
				if (Lengths_bio[n_abs_seg] == std::numeric_limits<unsigned short>::max() ||
				    Lengths_bio[n_abs_seg] > Lengths[n_abs_seg]) {
					const unsigned int seg_end_v = (v > 0) ? (v - 1) : 0;
					std::fprintf(stderr,
					             "[LENBIO_TRACE][build] sample=%s seg=%u start_v=%u end_v=%u len=%u len_bio=%u eff_len=%u\n",
					             name.c_str(), n_abs_seg, seg_start_v, seg_end_v,
					             Lengths[n_abs_seg], Lengths_bio[n_abs_seg], n_rel_var_eff);
				}
			}
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
			seg_start_v = v;
			n_abs_seg++;
		} else {
			// Same PER-BASEPAIR counting logic as Loop 1 - exclude siblings
			const bool is_sibling_l2 = (ctx.is_member && !ctx.is_anchor);
			// CRITICAL FIX: Match Loop 3's f_sca_eff logic (see Loop 1 comment for details)
			const bool f_sca_eff_l2 = f_sca && !(ctx.is_anchor && ctx.has_sca);
			const bool amb_for_count_l2 = is_sibling_l2 ? false : (f_het || f_sca_eff_l2);
			n_rel_unf += f_het_for_boundary;
			n_rel_sca += f_sca_for_boundary;
			n_abs_amb += amb_for_count_l2;
			n_rel_amb += amb_for_count_l2;
			n_abs_mis += f_mis;
			n_rel_var += var_len;       // real span for storage and indexing
			n_rel_var_eff += eff_len;   // effective span for boundary decisions
			v += var_len;
			// One-off debug guard to pinpoint len_bio divergence; remove after root cause is found.
			if (trace_lenbio_detail_on && trace_lenbio_sample_ok && !lenbio_detail_traced.empty()) {
				const unsigned int seg_span = v - seg_start_v;
				const bool bad_span = (n_rel_var != seg_span);
				const bool bad_lenbio = (n_rel_var_eff > n_rel_var);
				if ((bad_span || bad_lenbio) && !lenbio_detail_traced[n_abs_seg]) {
					std::fprintf(stderr,
					             "[LENBIO_DETAIL][build] sample=%s seg=%u v=%u seg_start=%u seg_span=%u len=%u len_bio=%u "
					             "var_len=%u eff_len=%u ss_idx=%d member=%d anchor=%d\n",
					             name.c_str(), n_abs_seg, v, seg_start_v, seg_span, n_rel_var, n_rel_var_eff,
					             var_len, eff_len, ss_idx, ctx.is_member ? 1 : 0, ctx.is_anchor ? 1 : 0);
					lenbio_detail_traced[n_abs_seg] = 1;
				}
			}
		}
	}
	Lengths[n_abs_seg] = n_rel_var;
	Lengths_bio[n_abs_seg] = n_rel_var_eff;
	{
		const char* trace_lenbio = std::getenv("SHAPEIT5_TRACE_LENBIO");
		if (trace_lenbio && trace_lenbio[0] != '\0' && trace_lenbio[0] != '0') {
			if (Lengths_bio[n_abs_seg] == std::numeric_limits<unsigned short>::max() ||
			    Lengths_bio[n_abs_seg] > Lengths[n_abs_seg]) {
				const unsigned int seg_end_v = (n_variants > 0) ? (n_variants - 1) : 0;
				std::fprintf(stderr,
				             "[LENBIO_TRACE][build] sample=%s seg=%u start_v=%u end_v=%u len=%u len_bio=%u eff_len=%u\n",
				             name.c_str(), n_abs_seg, seg_start_v, seg_end_v,
				             Lengths[n_abs_seg], Lengths_bio[n_abs_seg], n_rel_var_eff);
			}
		}
	}

	// DEBUG: Verify sum of Lengths equals n_variants
	if (name == "HG00107") {
		unsigned int sum_lengths = 0;
		for (unsigned int s = 0; s < n_segments; s++) sum_lengths += Lengths[s];
		std::fprintf(stderr, "[LOOP2_END] sample=%s sum_lengths=%u n_variants=%u n_segments=%u match=%d\n",
		             name.c_str(), sum_lengths, n_variants, n_segments, sum_lengths == n_variants);
	}

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
	unsigned int loop3_het_count = 0;  // DEBUG: count hets seen in loop 3
	for (unsigned int s = 0, a0 = 0, a1 = 0, a2 = 0, vabs = 0 ; s < n_segments ; s ++) {
		for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++) {
			unsigned int v_idx = vabs + vrel;
			unsigned char var_code = Variants[DIV2(v_idx)];
			bool f_sca = VAR_GET_SCA(MOD2(v_idx), var_code);
			bool f_het = VAR_GET_HET(MOD2(v_idx), var_code);
			SuperSiteContext ctx = getSuperSiteContext(v_idx);
			if (ctx.is_member && !ctx.is_anchor) continue;
			if (ctx.is_anchor) {
				f_het = ctx.has_het;
				f_sca = ctx.has_sca;  // Match Loop 1: use aggregated sca for anchors
			} else {
				f_sca = f_sca || ctx.has_sca;
			}
			const bool f_sca_eff = f_sca && !(ctx.is_anchor && ctx.has_sca);
			
			if (f_sca_eff) {
				for (unsigned int h = 0 ; h < HAP_NUMBER ; h ++) {
					bool allele = (h%2)?VAR_GET_HAP1(MOD2(v_idx), var_code):VAR_GET_HAP0(MOD2(v_idx), var_code);
					if (allele) HAP_SET(Ambiguous[a0], h);
				}
				orderedSegments[s] = 1;
				if (!debug::SUPERDEBUG_SAMPLENAME.empty() &&
				    name == debug::SUPERDEBUG_SAMPLENAME &&
				    static_cast<int>(v_idx) == debug::SUPERDEBUG_BP) {
					std::cout << "[SUPERDEBUG] BUILD_SCA locus=" << v_idx
					          << " seg=" << s
					          << " amb_idx=" << a0
					          << " amb_code=0x" << std::hex
					          << static_cast<int>(Ambiguous[a0]) << std::dec
					          << " ctx_anchor=" << ctx.is_anchor
					          << " ctx_member=" << ctx.is_member
					          << " ctx_has_het=" << ctx.has_het
					          << " ctx_has_sca=" << ctx.has_sca
					          << std::endl;
				}
			}
			a0 += (f_sca_eff||f_het);
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
				f_het = ctx.has_het;
				f_sca = ctx.has_sca;  // Match Loop 1: use aggregated sca for anchors
			} else {
				f_sca = f_sca || ctx.has_sca;
			}
			const bool f_sca_eff = f_sca && !(ctx.is_anchor && ctx.has_sca);
			
			if (f_het) {
				loop3_het_count++;
				// DEBUG: Log het counting crossing threshold - search for divergence
				if (trace_build_loops_on && name == "HG00107" && loop3_het_count >= 0 && loop3_het_count <= 2265) {
					std::fprintf(stderr, "[LOOP3_COUNT] v_idx=%u het#%u a1=%u f_het=%d is_member=%d is_anchor=%d ss_idx=%d\n",
					             v_idx, loop3_het_count, a1, f_het, ctx.is_member, ctx.is_anchor, ctx.ss_idx);
				}
				if (a1 >= n_ambiguous) {
					std::fprintf(stderr, "[LOOP3_OVERFLOW] sample=%s a1=%u n_ambiguous=%u\n", name.c_str(), a1, n_ambiguous);
					std::fprintf(stderr, "  v_idx=%u seg=%u vrel=%u vabs=%u\n", v_idx, s, vrel, vabs);
					std::fprintf(stderr, "  ctx: is_member=%d is_anchor=%d has_het=%d has_sca=%d ss_idx=%d\n",
					             ctx.is_member, ctx.is_anchor, ctx.has_het, ctx.has_sca, ctx.ss_idx);
					std::fprintf(stderr, "  flags: f_het=%d f_sca=%d f_sca_eff=%d\n", f_het, f_sca, f_sca_eff);
					std::fprintf(stderr, "  loop3_het_count=%u (hets processed so far)\n", loop3_het_count);
					std::fflush(stderr);
					std::abort();
				}
				for (unsigned int h = 0 ; h < HAP_NUMBER ; h ++) {
					bool allele = ((h>>n_unf)%2);
					if (allele) HAP_SET(Ambiguous[a1], h);
				}
				if (!debug::SUPERDEBUG_SAMPLENAME.empty() &&
				    name == debug::SUPERDEBUG_SAMPLENAME &&
				    static_cast<int>(v_idx) == debug::SUPERDEBUG_BP) {
					std::cout << "[SUPERDEBUG] BUILD_AMB locus=" << v_idx
					          << " seg=" << s
					          << " n_unf=" << n_unf
					          << " amb_idx=" << a1
					          << " amb_code=0x" << std::hex
					          << static_cast<int>(Ambiguous[a1]) << std::dec
					          << " ctx_anchor=" << ctx.is_anchor
					          << " ctx_member=" << ctx.is_member
					          << " ctx_has_het=" << ctx.has_het
					          << " ctx_has_sca=" << ctx.has_sca
					          << std::endl;
				}
				n_unf++;
			}
			a1 += (f_sca_eff||f_het);
		}
		vabs += Lengths[s];
	}

	if (debug::SUPERDEBUG_BP > 0 &&
	    !debug::SUPERDEBUG_SAMPLENAME.empty() &&
	    name == debug::SUPERDEBUG_SAMPLENAME &&
	    static_cast<int>(debug::SUPERDEBUG_BP) < static_cast<int>(n_variants)) {
		// On rebuild, report the amb_code at the target locus
		int amb_idx = 0;
		int locus = debug::SUPERDEBUG_BP;
		for (unsigned int s = 0, vabs2 = 0, a_idx = 0; s < n_segments; ++s) {
			for (unsigned int vrel = 0; vrel < Lengths[s]; ++vrel) {
				unsigned int v_idx2 = vabs2 + vrel;
				unsigned char var_code = Variants[DIV2(v_idx2)];
				bool f_sca = VAR_GET_SCA(MOD2(v_idx2), var_code);
				bool f_het = VAR_GET_HET(MOD2(v_idx2), var_code);
				SuperSiteContext ctx = getSuperSiteContext(v_idx2);
				if (ctx.is_member && !ctx.is_anchor) continue;
				if (ctx.is_anchor) {
					f_sca = f_sca || ctx.has_sca;
					f_het = ctx.has_het;
				}
				if (f_sca || f_het) {
					if (static_cast<int>(v_idx2) == locus) amb_idx = a_idx;
					++a_idx;
				}
			}
			vabs2 += Lengths[s];
		}
		unsigned char amb_code = (amb_idx < Ambiguous.size()) ? Ambiguous[amb_idx] : 0;
		std::cout << "[BUILD_SUMMARY] sample=" << name
		          << " locus=" << locus
		          << " amb_idx=" << amb_idx
		          << " amb_code=0x" << std::hex << static_cast<int>(amb_code) << std::dec
		          << " n_segments=" << n_segments
		          << " n_ambiguous=" << n_ambiguous
		          << std::endl;
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

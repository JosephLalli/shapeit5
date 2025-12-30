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

#include <containers/window_set.h>

using namespace std;

window_set::window_set() {
	W.clear();
}

window_set::~window_set() {
	W.clear();
}

void window_set::clear() {
	W.clear();
}

int window_set::size() {
	return W.size();
}

bool window_set::split(double min_length_cm, int left_index, int right_index,
                       vector < int > & idx_sta, vector < int > & idx_sto,
                       vector < int > & bio_idx_sta, vector < int > & bio_idx_sto,
                       vector < double > & ccm_sta, vector < double > & ccm_sto,
                       vector < int > & output, random_number_generator & rng) {
	int number_of_segments = right_index-left_index+1;
	int number_of_variants = bio_idx_sto[right_index] - bio_idx_sta[left_index] + 1;
	double length_of_region = ccm_sto[right_index] - ccm_sta[left_index];

	//A phasing window must (i) span >=4 segments, (ii) contain >= 100 variants and (iii) span more than "min_length_cm" cM
	if (number_of_segments < 4 || number_of_variants < 100 || length_of_region < min_length_cm) return false;
	else {
		int split_point = rng.getInt(number_of_segments/2) + number_of_segments/4 + 1;
		vector <  int > left_output, right_output;
		bool ret1 = split(min_length_cm, left_index, left_index + split_point, idx_sta, idx_sto, bio_idx_sta, bio_idx_sto, ccm_sta, ccm_sto, left_output, rng);
		bool ret2 = split(min_length_cm, left_index + split_point, right_index, idx_sta, idx_sto, bio_idx_sta, bio_idx_sto, ccm_sta, ccm_sto, right_output, rng);

		if (ret1 && ret2) {
			//succesful split, so operate it
			output = vector < int >(left_output.size() + right_output.size());
			std::copy(left_output.begin(), left_output.end(), output.begin());
			std::copy(right_output.begin(), right_output.end(), output.begin() + left_output.size());
		} else {
			//unsuccesful split, so return current coordinates
			output.clear();
			output.push_back(left_index);
			output.push_back(right_index);
		}
		return true;
	}
}


int window_set::build (variant_map & V, genotype * g, float min_window_size, random_number_generator & rng) {

	// PER-BASEPAIR: Returns true only for ambiguous biological positions
	// Siblings are excluded (return false) to match genotype::n_ambiguous semantics
	auto is_locus_ambiguous = [&](unsigned int locus) -> bool {
		bool is_amb = VAR_GET_AMB(MOD2(locus), g->Variants[DIV2(locus)]);
		genotype::SuperSiteContext ctx = g->getSuperSiteContext(locus);
		if (ctx.is_member) {
			if (!ctx.is_anchor) return false;  // Siblings excluded
			return ctx.has_het || ctx.has_sca;
		}
		return is_amb;
	};

	// PER-BASEPAIR: Returns true only for missing biological positions
	// Siblings are excluded (return false) to match PER-BASEPAIR semantics
	auto is_locus_missing = [&](unsigned int locus) -> bool {
		bool is_mis = VAR_GET_MIS(MOD2(locus), g->Variants[DIV2(locus)]);
		genotype::SuperSiteContext ctx = g->getSuperSiteContext(locus);
		if (ctx.is_member) {
			if (!ctx.is_anchor) return false;  // Siblings excluded
			return ctx.all_missing;
		}
		return is_mis;
	};

	//1. Mapping coordinates of each segment
	vector < unsigned int > loc_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > loc_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > bio_loc_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > bio_loc_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > amb_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > amb_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > mis_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > mis_siz = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > tra_idx = vector < unsigned int >(g->n_segments, 0);
	vector < unsigned int > tra_siz = vector < unsigned int >(g->n_segments, 0);
	vector < double > ccm_sta = vector < double >(g->n_segments, 0);
	vector < double > ccm_sto = vector < double >(g->n_segments, 0);
	vector < int > idx_sta = vector < int >(g->n_segments, 0);
	vector < int > idx_sto = vector < int >(g->n_segments, 0);
	vector < int > bio_idx_sta = vector < int >(g->n_segments, 0);
	vector < int > bio_idx_sto = vector < int >(g->n_segments, 0);

	// Prefix counts for missing/ambiguous across the full variant span.
	// The prefix only advances on biological loci (anchors), so siblings are naturally skipped.
	std::vector<unsigned int> amb_prefix(g->n_variants + 1, 0);
	std::vector<unsigned int> mis_prefix(g->n_variants + 1, 0);
	for (unsigned int locus = 0; locus < g->n_variants; ++locus) {
		amb_prefix[locus + 1] = amb_prefix[locus] + (is_locus_ambiguous(locus) ? 1u : 0u);
		mis_prefix[locus + 1] = mis_prefix[locus] + (is_locus_missing(locus) ? 1u : 0u);
	}

	unsigned int prev_dipcounts = 1, curr_dipcounts = 0;
	for (unsigned int s = 0, a = 0, t = 0, v = 0, vb = 0, m = 0 ; s < g->n_segments ; s ++) {
		//update a
		amb_idx[s] = a;
		for (unsigned int vrel = 0 ; vrel < g->Lengths[s] ; vrel ++) {
			unsigned int locus = v + vrel;
			if (is_locus_ambiguous(locus)) amb_siz[s]++;
		}
		a += amb_siz[s];
		//update m
		mis_idx[s] = m;
		for (unsigned int vrel = 0 ; vrel < g->Lengths[s] ; vrel ++) {
			unsigned int locus = v + vrel;
			if (is_locus_missing(locus)) mis_siz[s]++;
		}
		m += mis_siz[s];
		//update v
		loc_idx[s] = v;
		loc_siz[s] = g->Lengths[s];
		v += loc_siz[s];
		//update bio v
		bio_loc_idx[s] = vb;
		bio_loc_siz[s] = g->Lengths_bio.empty() ? g->Lengths[s] : g->Lengths_bio[s];
		vb += bio_loc_siz[s];
		//update idx
		idx_sta[s] = loc_idx[s];
		idx_sto[s] = loc_idx[s]+loc_siz[s]-1;
		bio_idx_sta[s] = bio_loc_idx[s];
		bio_idx_sto[s] = bio_loc_idx[s] + bio_loc_siz[s] - 1;
		//update ccm
		ccm_sta[s] = V.vec_pos[idx_sta[s]]->cm;
		ccm_sto[s] = V.vec_pos[idx_sto[s]]->cm;
		//update t
		tra_idx[s] = t;
		curr_dipcounts = g->countDiplotypes(g->Diplotypes[s]);
		tra_siz[s] = prev_dipcounts * curr_dipcounts;
		t += tra_siz[s];
		prev_dipcounts = curr_dipcounts;
	}

	//2. Reccursive split
	vector < int > output;
	output.push_back(0);
	output.push_back(g->n_segments-1);
	split(min_window_size, 0, g->n_segments-1, idx_sta, idx_sto, bio_idx_sta, bio_idx_sto, ccm_sta, ccm_sto, output, rng);
	int n_windows = output.size()/2;

	//3. Update coordinates
	W = vector < window > (n_windows);
	for (unsigned int w = 0 ; w < n_windows ; w ++) {
		W[w].start_segment = output[2*w+0];
		W[w].stop_segment = output[2*w+1];
		W[w].start_locus = loc_idx[W[w].start_segment];
		W[w].stop_locus = loc_idx[W[w].stop_segment] + loc_siz[W[w].stop_segment] - 1;

		// Adjust ambiguous indices within boundary segments so we don't count loci outside the window
		const unsigned int start_seg_locus = loc_idx[W[w].start_segment];
		const unsigned int stop_seg_locus = loc_idx[W[w].stop_segment];

		// Ambiguous range: use prefix counts so sibling splits (non-anchors) are ignored
		const unsigned int amb_before_start = amb_prefix[W[w].start_locus];
		const unsigned int amb_before_stop_plus1 = amb_prefix[W[w].stop_locus + 1];
		const unsigned int amb_count = amb_before_stop_plus1 > amb_before_start ? (amb_before_stop_plus1 - amb_before_start) : 0;
		if (amb_count == 0) {
			W[w].start_ambiguous = 0;
			W[w].stop_ambiguous = -1;
		} else {
			W[w].start_ambiguous = static_cast<int>(amb_before_start);
			W[w].stop_ambiguous = static_cast<int>(amb_before_start + amb_count - 1);
		}

		// Missing range: same prefix logic (anchors only, siblings skipped)
		const unsigned int mis_before_start = mis_prefix[W[w].start_locus];
		const unsigned int mis_before_stop_plus1 = mis_prefix[W[w].stop_locus + 1];
		const unsigned int mis_count = mis_before_stop_plus1 > mis_before_start ? (mis_before_stop_plus1 - mis_before_start) : 0;
		if (mis_count == 0) {
			W[w].start_missing = 0;
			W[w].stop_missing = -1;
		} else {
			W[w].start_missing = static_cast<int>(mis_before_start);
			W[w].stop_missing = static_cast<int>(mis_before_start + mis_count - 1);
		}

		// Transition bounds for this window: include all transitions for segments within [start_segment, stop_segment]
		// A segment's transitions correspond to the boundary (prev_segment -> segment), with segment 0 using a dummy prev count of 1.
		W[w].start_transition = tra_idx[W[w].start_segment];
		W[w].stop_transition = tra_idx[W[w].stop_segment] + tra_siz[W[w].stop_segment] - 1;
		if (std::getenv("SHAPEIT5_DEBUG_TRANS_PARITY")) {
			// Compute expected transitions within this window to sanity-check bounds.
			// The backward pass skips the first segment's boundary for windows starting at segment > 0
			// (SET_OTHER_TRANS requires curr_segment_index > segment_first).
			// So for non-zero start segments, we count from start_segment+1 instead.
			unsigned int expected_trans = 0;
			unsigned int loop_start = (W[w].start_segment == 0) ? 0 : W[w].start_segment + 1;
			unsigned int prev_dip = (loop_start == 0)
				? 1u
				: g->countDiplotypes(g->Diplotypes[loop_start - 1]);
			for (unsigned int s = loop_start; s <= W[w].stop_segment; ++s) {
				const unsigned int curr_dip = g->countDiplotypes(g->Diplotypes[s]);
				const unsigned int trans_here = prev_dip * curr_dip;
				expected_trans += trans_here;
				prev_dip = curr_dip;
			}
			if (W[w].stop_transition < W[w].start_transition ||
			    static_cast<unsigned int>(W[w].stop_transition - W[w].start_transition + 1) < expected_trans) {
				std::fprintf(stderr,
					"[WINDOW_TRANS_BOUNDS] sample=%s w=%u seg=[%u,%u] locus=[%u,%u] "
					"start_tr=%d stop_tr=%d span=%d expected_trans=%u\n",
					g ? g->name.c_str() : "?", w,
					W[w].start_segment, W[w].stop_segment,
					W[w].start_locus, W[w].stop_locus,
					W[w].start_transition, W[w].stop_transition,
					(W[w].stop_transition - W[w].start_transition + 1),
					expected_trans);
			} else if (w < 3) {
				std::fprintf(stderr,
					"[WINDOW_TRANS_BOUNDS] sample=%s w=%u seg=[%u,%u] locus=[%u,%u] "
					"start_tr=%d stop_tr=%d span=%d expected_trans=%u\n",
					g ? g->name.c_str() : "?", w,
					W[w].start_segment, W[w].stop_segment,
					W[w].start_locus, W[w].stop_locus,
					W[w].start_transition, W[w].stop_transition,
					(W[w].stop_transition - W[w].start_transition + 1),
					expected_trans);
			}
		}
	}
	return n_windows;
}

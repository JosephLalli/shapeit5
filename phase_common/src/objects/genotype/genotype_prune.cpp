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
#include <cstdlib>
#include <cmath>

using namespace std;

class Transition {
public:
	double prob;
	unsigned int idx;
	Transition() { prob = 0.0; idx = 0;}
	~Transition() {}
	bool operator < (const Transition & t) const {
		// Fuzzy comparison to handle FP rounding differences between biallelic and supersite modes.
		// This is needed because supersite processes anchor+siblings (multiple FP ops per locus)
		// while biallelic processes single variants (one FP op per locus), causing ULP-level
		// differences in transition probabilities that can affect sort order.
		// Deterministic sort is enabled by default; disable with SHAPEIT5_NONDETERMINISTIC_SORT=1.
		static bool use_epsilon = (std::getenv("SHAPEIT5_NONDETERMINISTIC_SORT") == nullptr);
		// NOTE: epsilon is treated as an *absolute* tolerance here so that tiny
		// differences at any scale are considered ties and resolved by the
		// transition index. This is intentionally conservative: it prefers
		// stability of sort order (and thus MCMC trajectories) over preserving
		// microscopic FP distinctions between biallelic and supersite paths.
		static const double epsilon = 1e-9; // absolute FP drift tolerance

		if (use_epsilon) {
			double diff = std::abs(prob - t.prob);
			// If probabilities differ by less than epsilon (absolute), use index as tie-breaker
			if (diff <= epsilon) {
				return idx < t.idx; // Deterministic tie-breaker
			}
		}
		// Standard comparison: sort by decreasing probability
		return prob > t.prob;
	}
};

class TransStatistics {
public:
	double entropy;
	unsigned int idx;
	bool merged;
	TransStatistics() {entropy = 1000; idx = -1; merged = false; }
	~TransStatistics() {};
	bool operator < (const TransStatistics & s) const {
		// Deterministic ordering: break ties on idx so sort results do not
		// depend on platform-specific quicksort behavior when entropies are
		// effectively equal (observed to accumulate subtle divergence across
		// epochs). Keep epsilon to avoid churn from FP noise.
		static const double epsilon = 1e-12;
		double diff = std::abs(entropy - s.entropy);
		double max_val = std::max(std::abs(entropy), std::abs(s.entropy));
		if (diff <= epsilon * (max_val > 0.0 ? max_val : 1.0)) {
			return idx < s.idx;
		}
		return entropy < s.entropy;
	}
};

void genotype::mapMerges(vector < double > & currProbs, double thresholdProbMass, vector < bool > & flagMerges) {
	vector < TransStatistics > vecTransStatistics = vector < TransStatistics > (n_segments - 1);
	vector < Transition > vecTransitions = vector < Transition > (4096);

	//Step0: initialize cursors
	unsigned int prev_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned int curr_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned char prev_dipcodes [64];
	makeDiplotypes(Diplotypes[0]);
	std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
	unsigned int toffset = prev_dipcount;
	unsigned int n_curr_transitions = 0;
	unsigned int aoffset = 0, voffset = 0;

	for (int s = 1 ; s < n_segments ; s++) {
		//Step1: update cursors (1)
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		n_curr_transitions = prev_dipcount * curr_dipcount;
		makeDiplotypes(Diplotypes[s]);

		//Step2: intialize transition statistics
		vecTransStatistics[s-1].idx = s;
		vecTransStatistics[s-1].merged = false;
		vecTransStatistics[s-1].entropy = 4096;

		auto isAmbiguous = [&](unsigned int locus) -> bool {
			bool is_amb = VAR_GET_AMB(MOD2(locus), Variants[DIV2(locus)]);
			SuperSiteContext ctx = getSuperSiteContext(locus);
			if (ctx.is_member) {
				if (!ctx.is_anchor) return false;
				is_amb = ctx.has_het || ctx.has_sca;
			}
			return is_amb;
		};
		//Step3: check number of variants in merged segment (guard against uint16_t overflow).
		unsigned int segment_length_raw = static_cast<unsigned int>(Lengths[s-1]) +
		                                  static_cast<unsigned int>(Lengths[s]);
		const bool segment_length_ok = (segment_length_raw < std::numeric_limits< unsigned short >::max());
		const bool overflow_suspect = (Lengths_bio[s-1] > Lengths[s-1]) || (Lengths_bio[s] > Lengths[s]);
		if (segment_length_ok && !overflow_suspect) {
			unsigned int n_ambiguous_merged = 0;
			for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++)
				if (isAmbiguous(voffset + vrel))
					n_ambiguous_merged++;
			//Step4: check number of ambiguous variants in merged segment
			if (n_ambiguous_merged < MAX_AMB) {
				//Step5: sort transitions by decreasing order
				for (int t = 0 ; t < n_curr_transitions ; t ++) {
					vecTransitions[t].prob = currProbs[toffset + t];
					vecTransitions[t].idx = t;
				}

				sort(vecTransitions.begin(), vecTransitions.begin() + n_curr_transitions);
				//Step6: compute transition entropy
				vecTransStatistics[s-1].entropy = 0.0;
				for (int t = 0 ; t < n_curr_transitions ; t ++) {
					double cProb = vecTransitions[t].prob;
					double lProb = -1.0 * ((cProb==0.0f)?0:log10(cProb));
					vecTransStatistics[s-1].entropy += cProb * lProb;
				}
				//Step7: check that 8 haplotypes capture lots of the cumulative probability mass
				double cumSumProbs = 0.0;
				vector < int > Mhaps = vector < int >(HAP_NUMBER * HAP_NUMBER, -1);
				int final_n_haps = 0;
				for (int t = 0, n_haps = 0 ; t < n_curr_transitions ; t ++) {
					cumSumProbs += vecTransitions[t].prob;
					unsigned int prev_dip = prev_dipcodes[vecTransitions[t].idx/curr_dipcount];
					unsigned int next_dip = curr_dipcodes[vecTransitions[t].idx%curr_dipcount];
					unsigned int prev_h0 = DIP_HAP0(prev_dip);
					unsigned int prev_h1 = DIP_HAP1(prev_dip);
					unsigned int next_h0 = DIP_HAP0(next_dip);
					unsigned int next_h1 = DIP_HAP1(next_dip);
					unsigned int merged_h0 = prev_h0 * HAP_NUMBER + next_h0;
					unsigned int merged_h1 = prev_h1 * HAP_NUMBER + next_h1;
					bool new_h0 = (Mhaps[merged_h0] < 0);
					bool new_h1 = ((Mhaps[merged_h1] < 0) && (merged_h0 != merged_h1));
					//if ((n_haps + new_h0 + new_h1) <= HAP_NUMBER) {
						if (new_h0) Mhaps[merged_h0] = n_haps++;
						if (new_h1) Mhaps[merged_h1] = n_haps++;
					//}
					if (n_haps == HAP_NUMBER && cumSumProbs > thresholdProbMass) {
						vecTransStatistics[s-1].merged = true;
					}
					final_n_haps = n_haps;
				}
			}
		}

		//Step8: update cursors (2)
		for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++) if (isAmbiguous(voffset + vrel)) aoffset++;
		voffset += Lengths[s-1];
		std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
		prev_dipcount = curr_dipcount;
		toffset += n_curr_transitions;
	}
	//Step9: map acceptable merges
	sort(vecTransStatistics.begin(), vecTransStatistics.end());
	flagMerges = vector < bool > (n_segments+1, false);
	for (unsigned int s = 0 ; s < vecTransStatistics.size() ; s ++) {
		bool no_adjacent_merges = !flagMerges[vecTransStatistics[s].idx-1] && !flagMerges[vecTransStatistics[s].idx+1];
		bool can_be_merged = vecTransStatistics[s].merged;
		bool will_merge = no_adjacent_merges && can_be_merged;
		flagMerges[vecTransStatistics[s].idx] = will_merge;
	}
}

void genotype::performMerges(vector < double > & currProbs, vector < bool > & flagMerges) {
	// NOTE: Mhaps lane remapping is intentionally free-form. Given a pair of
	// segments s-1 and s, we build a merged hap basis by scanning the sorted
	// transition list and assigning up to 8 distinct (prev_h, next_h) pairs
	// to lanes 0..7. This means that, even when two different representations
	// (biallelic vs supersite) describe the same biological site, Mhaps is
	// allowed to pick different lane permutations if the underlying transition
	// matrices or diplotype masks differ. Some tests (e.g.
	// test_supersite_expansion_epochs*) rely on specially constructed dummy-alt
	// data where those inputs *should* match, so that any observed divergence
	// in lane assignments is interpreted as a supersite bug rather than a rule
	// that Mhaps must enforce.
	vector < Transition > vecTransitions = vector < Transition > (4096);

	//Step0: initialize duplicates
	vector < unsigned char > Ambiguous2 = vector < unsigned char > (Ambiguous.size(), 0);
	vector < unsigned long > Diplotypes2;
	vector < unsigned short > Lengths2;
	vector < unsigned short > Lengths_bio2;
	unsigned int n_segments2 = n_segments;
	for (int s = 0 ; s < flagMerges.size() ; s++) n_segments2 -= flagMerges[s];
	Diplotypes2.reserve(n_segments2);
	Lengths2.reserve(n_segments2); // length of each segment
	Lengths_bio2.reserve(n_segments2); // length of each segment (siblings excluded)

	//Step1: initialize cursors
	unsigned int prev_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned int curr_dipcount = countDiplotypes(Diplotypes[0]);
	unsigned char prev_dipcodes [64];
	makeDiplotypes(Diplotypes[0]);
	std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
	unsigned int toffset = prev_dipcount;
	unsigned int n_curr_transitions = 0;
	unsigned int aoffset = 0, voffset = 0;

	auto isAmbiguous = [&](unsigned int locus) -> bool {
		bool is_amb = VAR_GET_AMB(MOD2(locus), Variants[DIV2(locus)]);
		SuperSiteContext ctx = getSuperSiteContext(locus);
		if (ctx.is_member) {
			// Siblings never contribute separate ambiguous events; anchors are
			// ambiguous if and only if they were marked AMB at build time,
			// which is driven by heterozygosity only (supersite scaffolds are
			// ignored to keep Ambiguous indexing in sync with bial).
			if (!ctx.is_anchor) return false;
		}
		return is_amb;
	};

	for (int s = 1 ; s < flagMerges.size() -1 ; s ++) {
		//Step1: update cursors (1)
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		n_curr_transitions = prev_dipcount * curr_dipcount;
		makeDiplotypes(Diplotypes[s]);

		//case1: merge to be done
		if (flagMerges[s]) {

			// CRITICAL: Pre-compute ambiguous anchor counts for correct segment membership.
			// We need this because arel (ambiguous anchor index) advances at a different
			// rate than vrel (all-variant index) when siblings are present.
			unsigned int n_amb_first = 0, n_amb_second = 0;
			for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++)
				if (isAmbiguous(voffset + vrel)) n_amb_first++;
			for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++)
				if (isAmbiguous(voffset + Lengths[s-1] + vrel)) n_amb_second++;

			// Guard against silent overflow when merging segments (Lengths_bio is uint16_t).
			const unsigned int len_sum = static_cast<unsigned int>(Lengths[s-1]) +
			                             static_cast<unsigned int>(Lengths[s]);
			const unsigned int len_bio_sum = static_cast<unsigned int>(Lengths_bio[s-1]) +
			                                 static_cast<unsigned int>(Lengths_bio[s]);
			if (len_sum > std::numeric_limits<unsigned short>::max() ||
			    len_bio_sum > std::numeric_limits<unsigned short>::max()) {
				std::fprintf(stderr,
				             "[PRUNE_LEN_OVERFLOW] sample=%s seg=%d len=%u+%u=%u len_bio=%u+%u=%u\n",
				             name.c_str(), s,
				             Lengths[s-1], Lengths[s], len_sum,
				             Lengths_bio[s-1], Lengths_bio[s], len_bio_sum);
				std::fflush(stderr);
				std::abort();
			}
			Lengths2.push_back(static_cast<unsigned short>(len_sum));
			Lengths_bio2.push_back(static_cast<unsigned short>(len_bio_sum));
			Diplotypes2.push_back(0x0000000000000000UL);
			for (int t = 0 ; t < n_curr_transitions ; t ++) { vecTransitions[t].prob = currProbs[toffset + t]; vecTransitions[t].idx = t; }
			sort(vecTransitions.begin(), vecTransitions.begin() + n_curr_transitions);

			int n_haps = 0;
			int n_skipped = 0;
			vector < int > Mhaps = vector < int >(HAP_NUMBER * HAP_NUMBER, -1);
			for (int t = 0 ; t < n_curr_transitions ; t ++) {
				unsigned int prev_dip = prev_dipcodes[vecTransitions[t].idx/curr_dipcount];
				unsigned int next_dip = curr_dipcodes[vecTransitions[t].idx%curr_dipcount];
				unsigned int prev_h0 = DIP_HAP0(prev_dip);
				unsigned int prev_h1 = DIP_HAP1(prev_dip);
				unsigned int next_h0 = DIP_HAP0(next_dip);
				unsigned int next_h1 = DIP_HAP1(next_dip);
				unsigned int merged_h0 = prev_h0 * HAP_NUMBER + next_h0;
				unsigned int merged_h1 = prev_h1 * HAP_NUMBER + next_h1;
				bool new_h0 = (Mhaps[merged_h0] < 0);
				bool new_h1 = ((Mhaps[merged_h1] < 0) && (merged_h0 != merged_h1));
					if ((n_haps + new_h0 + new_h1) <= HAP_NUMBER) {
						if (new_h0) {
							Mhaps[merged_h0] = n_haps;
							// CRITICAL: Use vrel (all-variant index) to determine segment membership.
							// vrel < Lengths[s-1] means we're in segment s-1, otherwise segment s.
							// This works correctly for both biallelic and supersite because Lengths
							// includes siblings, while arel only counts ambiguous anchors.
							for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
								if (isAmbiguous(voffset+vrel)) {
									bool in_first_segment = (arel < n_amb_first);
									unsigned int old_hap = in_first_segment ? prev_h0 : next_h0;
									bool was_set = HAP_GET(Ambiguous[aoffset+arel], old_hap);
									if (was_set) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h0]);
									arel ++;
								}
							}
							n_haps ++;
						}
						if (new_h1) {
							Mhaps[merged_h1] = n_haps;
							// CRITICAL: Use vrel (all-variant index) to determine segment membership.
							// vrel < Lengths[s-1] means we're in segment s-1, otherwise segment s.
							// This works correctly for both biallelic and supersite because Lengths
							// includes siblings, while arel only counts ambiguous anchors.
							for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
								if (isAmbiguous(voffset+vrel)) {
									bool in_first_segment = (arel < n_amb_first);
									unsigned int old_hap = in_first_segment ? prev_h1 : next_h1;
									bool was_set = HAP_GET(Ambiguous[aoffset+arel], old_hap);
									if (was_set) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h1]);
									arel ++;
								}
							}
						n_haps ++;
					}
					DIP_SET(Diplotypes2.back(), Mhaps[merged_h0] * HAP_NUMBER + Mhaps[merged_h1]);
				} else {
					n_skipped++;
				}
			}
			if (n_haps != HAP_NUMBER) {
				std::fprintf(stderr, "[PRUNE_ERROR] %s seg=%d n_haps=%d HAP_NUMBER=%d n_curr_transitions=%d n_skipped=%d prev_dipcount=%d curr_dipcount=%d\n",
							 name.c_str(), s, n_haps, HAP_NUMBER, n_curr_transitions, n_skipped, prev_dipcount, curr_dipcount);
				std::fprintf(stderr, "[PRUNE_ERROR] Lengths[s-1]=%u Lengths[s]=%u n_amb_first=%u n_amb_second=%u\n",
							 Lengths[s-1], Lengths[s], n_amb_first, n_amb_second);
				std::fprintf(stderr, "[PRUNE_ERROR] Mhaps filled: ");
				for (int i = 0; i < HAP_NUMBER * HAP_NUMBER; i++) {
					if (Mhaps[i] >= 0) std::fprintf(stderr, "[%d]=%d ", i, Mhaps[i]);
				}
				std::fprintf(stderr, "\n");
			}
			assert(n_haps == HAP_NUMBER);
		//Case2: no merge to be done, push last segment
		} else if (!flagMerges[s-1]) {
			//cout << name << " C " << aoffset << endl;
			for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths[s-1] ; vrel ++) {
				if (isAmbiguous(voffset + vrel)) {
					Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
					arel ++;
				}
			}
			Lengths2.push_back(Lengths[s-1]);
			Lengths_bio2.push_back(Lengths_bio[s-1]);
			Diplotypes2.push_back(Diplotypes[s-1]);
		}

		//Update cursors
		for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++) if (isAmbiguous(voffset + vrel)) aoffset++;
		voffset += Lengths[s-1];
		std::copy(curr_dipcodes, curr_dipcodes+curr_dipcount, prev_dipcodes);
		prev_dipcount = curr_dipcount;
		toffset += n_curr_transitions;
	}
	if (!flagMerges[flagMerges.size()-2]) {
		for (unsigned int vrel = 0, arel = 0 ; vrel < Lengths.back() ; vrel ++) {
			if (isAmbiguous(voffset + vrel)) {
				Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
				arel ++;
			}
		}
		Lengths2.push_back(Lengths.back());
		Lengths_bio2.push_back(Lengths_bio.back());
		Diplotypes2.push_back(Diplotypes.back());
	}
	//free();
	Ambiguous = Ambiguous2;
	Diplotypes = Diplotypes2;
	Lengths = Lengths2;
	Lengths_bio = Lengths_bio2;
	n_segments = n_segments2;
	n_transitions = countTransitions();

	/*
	for (int s = 0 ; s < n_segments ; s++) {
		cout << name << " " << s << " " << isOrdered(Diplotypes[s]) << endl;
	}
	 */
}

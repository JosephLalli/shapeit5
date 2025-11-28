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
		// Enable with SHAPEIT5_DETERMINISTIC_SORT=1 environment variable.
		// This is needed because supersite processes anchor+siblings (multiple FP ops per locus)
		// while biallelic processes single variants (one FP op per locus), causing ULP-level
		// differences in transition probabilities that can affect sort order.
		static bool use_epsilon = (std::getenv("SHAPEIT5_DETERMINISTIC_SORT") != nullptr);
		static const double epsilon = 1e-12; // ~40 ULPs for values around 1e-9

		if (use_epsilon) {
			double diff = std::abs(prob - t.prob);
			double max_val = std::max(std::abs(prob), std::abs(t.prob));
			// If probabilities differ by less than epsilon (relative), use index as tie-breaker
			if (diff <= epsilon * max_val) {
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
	bool operator < (const TransStatistics & s) const { return entropy < s.entropy; }
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
		//Step3: check number of variants in merged segment
		unsigned int segment_length = Lengths[s-1] + Lengths[s];
		if (segment_length < std::numeric_limits< unsigned short >::max()) {
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

				// LOG TRANSITION PROBABILITIES BEFORE SORT (Scenario 8, prune iterations only)
				bool is_test_sample = (name.find("_sample") != std::string::npos);
				static int prune_iter_count = 0;
				if (is_test_sample && n_segments == 3 && s >= 1 && s <= 2) {
					std::fprintf(stderr, "[TRANS_PROBS_PRESORT] sample=%s seg=%d n_trans=%d prune_count=%d\n",
					             name.c_str(), s, n_curr_transitions, prune_iter_count);
					int max_print = (n_curr_transitions < 20) ? n_curr_transitions : 20;
					for (int t = 0; t < max_print; t++) {
						std::fprintf(stderr, "  [%d] idx=%u prob=%.20g (hex:%a)\n",
						             t, vecTransitions[t].idx, vecTransitions[t].prob, vecTransitions[t].prob);
					}
					if (s == 2) prune_iter_count++;
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
						std::fprintf(stderr, "[MAP_MERGE_DEBUG] %s seg=%d n_haps=%d cumSumProbs=%.6f thresholdProbMass=%.6f t=%d/%d\n",
									 name.c_str(), s, n_haps, cumSumProbs, thresholdProbMass, t, n_curr_transitions);
					}
					final_n_haps = n_haps;
				}
				if (final_n_haps != HAP_NUMBER && segment_length < std::numeric_limits< unsigned short >::max() && n_ambiguous_merged < MAX_AMB) {
					std::fprintf(stderr, "[MAP_MERGE_SKIP] %s seg=%d final_n_haps=%d (not 8), merged=%s\n",
								 name.c_str(), s, final_n_haps, vecTransStatistics[s-1].merged ? "true" : "false");
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
	bool enable_debug = (!debug::SUPERDEBUG_SAMPLENAME.empty() && name == debug::SUPERDEBUG_SAMPLENAME) ||
						(name.find("_sample") != std::string::npos);
	if (enable_debug) {
		std::fprintf(stderr, "[MAP_ENTROPY_SORT] %s sorted by entropy (size=%zu):\n", name.c_str(), vecTransStatistics.size());
		for (unsigned int s = 0; s < vecTransStatistics.size(); s++) {
			std::fprintf(stderr, "  [%u] idx=%u entropy=%.10f merged=%d\n",
						 s, vecTransStatistics[s].idx, (double)vecTransStatistics[s].entropy,
						 vecTransStatistics[s].merged ? 1 : 0);
			std::fflush(stderr);
		}
	}
	for (unsigned int s = 0 ; s < vecTransStatistics.size() ; s ++) {
		bool no_adjacent_merges = !flagMerges[vecTransStatistics[s].idx-1] && !flagMerges[vecTransStatistics[s].idx+1];
		bool can_be_merged = vecTransStatistics[s].merged;
		bool will_merge = no_adjacent_merges && can_be_merged;
		flagMerges[vecTransStatistics[s].idx] = will_merge;
		if (enable_debug) {
			std::fprintf(stderr, "[MAP_FINAL] idx=%u no_adj=%d can_merge=%d will_merge=%d\n",
						 vecTransStatistics[s].idx, no_adjacent_merges, can_be_merged, will_merge);
		}
	}
}

void genotype::performMerges(vector < double > & currProbs, vector < bool > & flagMerges) {
	const bool superdebug = (!debug::SUPERDEBUG_SAMPLENAME.empty() && name == debug::SUPERDEBUG_SAMPLENAME) ||
							(name.find("_sample") != std::string::npos);
	if (superdebug) {
		std::cout << "[PRUNE_DEBUG] " << name << " enter performMerges n_segments=" << n_segments
				  << " flagMerges.size=" << flagMerges.size()
				  << " Ambiguous.size=" << Ambiguous.size()
				  << " Lengths[0]=" << (Lengths.empty()?0:Lengths[0])
				  << " Lengths_bio[0]=" << (Lengths_bio.empty()?0:Lengths_bio[0])
				  << " flagM[0]=" << (flagMerges.empty()?0:flagMerges[0])
				  << std::endl;
	}
	auto logAmbiguousSite = [&](const char* ctx, int seg_idx, unsigned int vabs, unsigned int vrel,
								unsigned int arel, unsigned char amb_code, unsigned int l_bio_left,
								unsigned int l_bio_right) {
		if (!superdebug) return;
		std::cout << "[PRUNE_DEBUG] " << name
				  << " ctx=" << ctx
				  << " seg=" << seg_idx
				  << " vabs=" << vabs
				  << " vrel=" << vrel
				  << " arel=" << arel
				  << " amb_code=0x" << std::hex << static_cast<int>(amb_code) << std::dec
				  << " len_bio_left=" << l_bio_left
				  << " len_bio_right=" << l_bio_right
				  << std::endl;
	};
	auto logFlagVector = [&](const char* ctx, const std::vector<bool>& flags) {
		if (!superdebug) return;
		std::cout << "[PRUNE_DEBUG] " << name << " " << ctx << " flagMerges:";
		for (bool f : flags) std::cout << (f ? "1" : "0");
		std::cout << std::endl;
	};

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
			if (!ctx.is_anchor) return false;
			is_amb = ctx.has_het || ctx.has_sca;
		}
		return is_amb;
	};

	logFlagVector("pre-loop", flagMerges);
	for (int s = 1 ; s < flagMerges.size() -1 ; s ++) {
		//Step1: update cursors (1)
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		n_curr_transitions = prev_dipcount * curr_dipcount;
		makeDiplotypes(Diplotypes[s]);

		//case1: merge to be done
		if (flagMerges[s]) {
			logFlagVector("merge-at", flagMerges);

			// CRITICAL: Pre-compute ambiguous anchor counts for correct segment membership.
			// We need this because arel (ambiguous anchor index) advances at a different
			// rate than vrel (all-variant index) when siblings are present.
			unsigned int n_amb_first = 0, n_amb_second = 0;
			for (unsigned int vrel = 0 ; vrel < Lengths[s-1] ; vrel ++)
				if (isAmbiguous(voffset + vrel)) n_amb_first++;
			for (unsigned int vrel = 0 ; vrel < Lengths[s] ; vrel ++)
				if (isAmbiguous(voffset + Lengths[s-1] + vrel)) n_amb_second++;

			if (superdebug) {
				std::fprintf(stderr, "[PRUNE_MERGE] seg=%d Lengths[s-1]=%u Lengths[s]=%u "
							 "Lengths_bio[s-1]=%u Lengths_bio[s]=%u "
							 "n_amb_first=%u n_amb_second=%u\n",
							 s, Lengths[s-1], Lengths[s], Lengths_bio[s-1], Lengths_bio[s],
							 n_amb_first, n_amb_second);
			}

			Lengths2.push_back(Lengths[s-1]+Lengths[s]);
			Lengths_bio2.push_back(Lengths_bio[s-1]+Lengths_bio[s]);
			Diplotypes2.push_back(0x0000000000000000UL);
			for (int t = 0 ; t < n_curr_transitions ; t ++) { vecTransitions[t].prob = currProbs[toffset + t]; vecTransitions[t].idx = t; }
			sort(vecTransitions.begin(), vecTransitions.begin() + n_curr_transitions);
			if (superdebug) {
				std::fprintf(stderr, "[PRUNE_DEBUG] %s seg=%d Top 10 Transitions:\n", name.c_str(), s);
				for (int t = 0; t < 10 && t < n_curr_transitions; ++t) {
					std::fprintf(stderr, "  t=%d prob=%.20g idx=%u\n", t, vecTransitions[t].prob, vecTransitions[t].idx);
				}
			}

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
						if (superdebug) {
							std::fprintf(stderr, "[MERGE_AMB_H0] %s seg=%d merged_h0=%u->%d prev_h0=%u next_h0=%u Lengths[s-1]=%u Lengths[s]=%u\n",
										 name.c_str(), s, merged_h0, n_haps, prev_h0, next_h0, Lengths[s-1], Lengths[s]);
						}
						for (unsigned int vrel = 0, arel = 0 ; vrel < (Lengths[s-1]+Lengths[s]) ; vrel ++) {
							if (isAmbiguous(voffset+vrel)) {
								bool in_first_segment = (vrel < Lengths[s-1]);
								unsigned int old_hap = in_first_segment ? prev_h0 : next_h0;
								bool was_set = HAP_GET(Ambiguous[aoffset+arel], old_hap);
								if (was_set) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h0]);
								if (superdebug) {
									std::fprintf(stderr, "  [H0] vrel=%u arel=%u locus=%u in_first=%d old_hap=%u was_set=%d new_val=%02x\n",
												 vrel, arel, voffset+vrel, in_first_segment, old_hap, was_set, Ambiguous2[aoffset+arel]);
								}
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
								bool in_first_segment = (vrel < Lengths[s-1]);
								unsigned int old_hap = in_first_segment ? prev_h1 : next_h1;
								if (HAP_GET(Ambiguous[aoffset+arel], old_hap)) HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h1]);
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
					logAmbiguousSite("copy", s, voffset + vrel, vrel, arel,
									 Ambiguous[aoffset+arel], Lengths_bio[s-1], 0);
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
	if (superdebug) {
		std::cout << "[PRUNE_DEBUG] " << name
				  << " exit performMerges n_segments=" << n_segments2
				  << " Ambiguous.size=" << Ambiguous2.size()
				  << " copied_amb=" << Ambiguous.size()
				  << " Lengths_bio[0]=" << (Lengths_bio2.empty()?0:Lengths_bio2[0])
				  << std::endl;
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

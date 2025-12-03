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
#include <models/supersite_trace_utils.h>
#include <cstdlib>

using namespace std;

void genotype::sample(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	const bool use_forward = (sample_rng.getDouble() < 0.5f);
	const bool trace_sample = []() {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
		return env && env[0] != '\0' && env[0] != '0';
	}();
	if (trace_sample) {
		std::fprintf(stderr, "[SAMPLE_PATH] sample=%s direction=%s n_segments=%u\n",
		             name.c_str(), use_forward ? "forward" : "backward", n_segments);
	}
	if (use_forward) sampleForward(CurrentTransProbabilities, CurrentMissingProbabilities);
	else sampleBackward(CurrentTransProbabilities, CurrentMissingProbabilities);
}

void genotype::sampleForward(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	if (revert_buffer_fix) {
		double sumProbs = 0.0;
		unsigned int prev_sampled = 0;
		unsigned int curr_dipcount = 0, prev_dipcount = 1;
		vector < double > currProbs = vector < double > (64, 0.0);
		vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
		for (unsigned int s = 0, toffset = 0 ; s < n_segments ; s ++) {
			sumProbs = 0.0;
			curr_dipcount = countDiplotypes(Diplotypes[s]);
			for (unsigned int tabs = toffset + prev_sampled*curr_dipcount, trel = 0 ; trel < curr_dipcount ; ++trel, ++tabs)
				sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
			prev_sampled = sample_rng.sample(currProbs, sumProbs);
			makeDiplotypes(Diplotypes[s]);
			DipSampled[s] = curr_dipcodes[prev_sampled];
			toffset += prev_dipcount * curr_dipcount;
			prev_dipcount = curr_dipcount;
		}
		make(DipSampled, CurrentMissingProbabilities);
		return;
	}

	double sumProbs = 0.0;
	unsigned int prev_sampled = 0;
	unsigned int curr_dipcount = 0, prev_dipcount = 1;
	vector < double > currProbs;
	currProbs.reserve(64);
	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
	const bool trace_sample = []() {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
		return env && env[0] != '\0' && env[0] != '0';
	}();
	const size_t trans_buf_size = CurrentTransProbabilities.size();
	for (unsigned int s = 0, toffset = 0 ; s < n_segments ; s ++) {
		sumProbs = 0.0;
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		currProbs.resize(curr_dipcount);
		for (unsigned int tabs = toffset + prev_sampled*curr_dipcount, trel = 0 ; trel < curr_dipcount ; ++trel, ++tabs) {
			if (tabs >= trans_buf_size) {
				std::fprintf(stderr,
				             "[TRANS_READ_OOB][sampleForward] sample=%s segment=%u tabs=%u buf_size=%zu toffset=%u prev_sampled=%u curr_dipcount=%u prev_dipcount=%u n_segments=%u n_transitions=%u\n",
				             name.c_str(), s, tabs, trans_buf_size, toffset, prev_sampled, curr_dipcount, prev_dipcount, n_segments, n_transitions);
				std::fflush(stderr);
				std::abort();
			}
			sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
		}

		// SAMPLE TRACE: Log sampling details before drawing
		if (supersite_trace_enabled() && s == 0) {
			std::fprintf(stderr, "[SAMPLE_FWD] sample=%s segment=%u sumProbs=%.15f\n",
			             name.c_str(), s, sumProbs);
			double cumulative = 0.0;
			for (unsigned int i = 0; i < curr_dipcount; i++) {
				cumulative += currProbs[i];
				std::fprintf(stderr, "  dip[%u] prob=%.15f norm=%.15f cum=%.15f\n",
				             i, currProbs[i], currProbs[i]/sumProbs, cumulative/sumProbs);
			}
		}
		// Targeted transition dump for segment 3 divergence
		if (trace_sample && n_segments > 1 && s == 3) {
			std::array<unsigned char, 64> prev_codes{};
			if (s >= 1) {
				unsigned int idx = 0;
				for (unsigned int d = 0; d < 64; ++d) {
					if (DIP_GET(Diplotypes[s-1], d)) prev_codes[idx++] = static_cast<unsigned char>(d);
				}
			}
			unsigned char prev_dipcode = (s>=1 && prev_sampled < prev_dipcount) ? prev_codes[prev_sampled] : 0;
			std::fprintf(stderr, "[SAMPLE_DEBUG_TRANS] sample=%s dir=forward seg=%u prev_seg_code=%u prev_sampled=%u prev_dipcount=%u toffset=%u dip_mask=0x%016llx\n",
			             name.c_str(), s, (unsigned)prev_dipcode, prev_sampled, prev_dipcount, toffset,
			             static_cast<unsigned long long>(Diplotypes[s]));
			for (unsigned int i = 0; i < curr_dipcount; ++i) {
				unsigned int tabs = toffset + prev_sampled * curr_dipcount + i;
				std::fprintf(stderr, "  trans[%u] code=%u prob=%.15f\n",
				             tabs, static_cast<unsigned>(curr_dipcodes[i]), CurrentTransProbabilities[tabs]);
			}
		}

			prev_sampled = sample_rng.sample(currProbs, sumProbs);

		makeDiplotypes(Diplotypes[s]);
		unsigned char selected_dipcode = curr_dipcodes[prev_sampled];

		// Targeted debug for dipcode divergence (burn3, repeat_factor=8, segment 3)
		if (trace_sample && n_segments > 1 && s == 3) {
			std::fprintf(stderr, "[SAMPLE_DEBUG_DIP] sample=%s dir=forward seg=%u dipcount=%u toffset=%u dip_mask=0x%016llx\n",
			             name.c_str(), s, curr_dipcount, toffset, static_cast<unsigned long long>(Diplotypes[s]));
			double cumulative = 0.0;
			for (unsigned int i = 0; i < curr_dipcount; ++i) {
				cumulative += currProbs[i];
				std::fprintf(stderr, "  dip[%u] code=%u prob=%.15f norm=%.15f cum=%.15f\n",
				             i, static_cast<unsigned>(curr_dipcodes[i]), currProbs[i], currProbs[i]/sumProbs, cumulative/sumProbs);
			}
			unsigned char hap0 = (selected_dipcode >> 3);
			unsigned char hap1 = (selected_dipcode & 7);
			std::fprintf(stderr, "  selected_idx=%u dipcode=%u hap0=%u hap1=%u\n",
			             prev_sampled, (unsigned)selected_dipcode, (unsigned)hap0, (unsigned)hap1);
		}

		// SAMPLE TRACE: Log selected index after drawing (all segments when tracing)
		if (trace_sample) {
			unsigned char hap0 = (selected_dipcode >> 3);
			unsigned char hap1 = (selected_dipcode & 7);
			std::fprintf(stderr, "[SAMPLE_FWD_PICK] sample=%s seg=%u selected_idx=%u dipcode=%u hap0=%u hap1=%u\n",
			             name.c_str(), s, prev_sampled, static_cast<unsigned>(selected_dipcode),
			             static_cast<unsigned>(hap0), static_cast<unsigned>(hap1));
		}

		DipSampled[s] = selected_dipcode;
		if (!DIP_GET(Diplotypes[s], selected_dipcode)) {
			std::fprintf(stderr,
			             "[SAMPLE_DIPC_INVALID][forward] sample=%s seg=%u dipcode=%u dipcount=%u mask=0x%016llx\n",
			             name.c_str(), s, static_cast<unsigned>(selected_dipcode), countDiplotypes(Diplotypes[s]),
			             static_cast<unsigned long long>(Diplotypes[s]));
			std::fflush(stderr);
			std::abort();
		}
		toffset += prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}
	if (trace_sample) {
		std::fprintf(stderr, "[SAMPLE_DEBUG] sample=%s direction=forward n_segments=%u dip_mask=0x%016llx DipSampled0=%u\n",
		             name.c_str(),
		             n_segments,
		             static_cast<unsigned long long>(n_segments ? Diplotypes[0] : 0),
		             n_segments ? static_cast<unsigned>(DipSampled[0]) : 0);
	}
	make(DipSampled, CurrentMissingProbabilities);

	// HYPOTHESIS 2 DEBUGGING
	if (super_sites && !debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME) {
		for (const auto& ss : *super_sites) {
			if ((int)ss.global_site_id == debug::SUPERDEBUG_BP) {
				debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After make() in sampleForward");
				break;
			}
		}
	}

	// CRITICAL: Post-HMM supersite projection
	// This is the ONLY place where projectSupersites() is called in the entire codebase.
	// It runs after every sample() operation in all MCMC stages (burn-in, pruning, main).
	// The projection maps anchor phasing results to member split variants.
	// No additional projection is needed during finalization - this call ensures it's complete.
	if (super_sites) projectSupersites();

	// HYPOTHESIS 2 DEBUGGING
	if (super_sites && !debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME) {
		for (const auto& ss : *super_sites) {
			if ((int)ss.global_site_id == debug::SUPERDEBUG_BP) {
				debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After projectSupersites() in sampleForward");
				break;
			}
		}
	}
} // <-- Missing brace added here

void genotype::sampleBackward(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {

	if (revert_buffer_fix) {
		double sumProbs = 0.0;
		int next_sampled = -1;
		unsigned int curr_dipcount = 0, next_dipcount = countDiplotypes(Diplotypes[n_segments - 1]);
		vector < double > currProbs = vector < double > (64 * 64, 0.0);
		vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);

		for (int s = n_segments - 2, toffset = n_transitions ; s >= 0 ; s --) {
			sumProbs = 0.0;
			curr_dipcount = countDiplotypes(Diplotypes[s]);
			toffset -= next_dipcount * curr_dipcount;

			if (next_sampled >= 0) {
				currProbs.resize(64);
				for (unsigned int tabs = toffset+next_sampled, trel = 0 ; trel < curr_dipcount ; ++trel, tabs += next_dipcount)
					sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
				next_sampled = sample_rng.sample(currProbs, sumProbs);
				makeDiplotypes(Diplotypes[s]);
				DipSampled[s] = curr_dipcodes[next_sampled];
			} else {
				for (unsigned int tabs = toffset, trel = 0 ; tabs < n_transitions ; ++trel, ++tabs)
					sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
				next_sampled = sample_rng.sample(currProbs, sumProbs);
				makeDiplotypes(Diplotypes[s+1]);
				DipSampled[s+1] = curr_dipcodes[next_sampled % next_dipcount];
				makeDiplotypes(Diplotypes[s]);
				next_sampled = next_sampled / next_dipcount;
				DipSampled[s] = curr_dipcodes[next_sampled];
			}
			next_dipcount = curr_dipcount;
		}
		make(DipSampled, CurrentMissingProbabilities);
		return;
	}

	double sumProbs = 0.0;
	int next_sampled = -1;
	unsigned int curr_dipcount = 0, next_dipcount = countDiplotypes(Diplotypes[n_segments - 1]);
	vector < double > currProbs;
	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
	const size_t trans_buf_size = CurrentTransProbabilities.size();

	const bool trace_sample = []() {
		const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
        return env && env[0] != '\0' && env[0] != '0';
	}();
	if (trace_sample) {
		std::fprintf(stdout, "[SAMPLE_BWD_ENTRY] sample=%s n_segments=%u\n", name.c_str(), n_segments);
		std::fflush(stdout);
	}

	if (n_segments == 1) {
		// Single-segment genotypes skip the backward loop entirely. Sample directly
		// from the first transition block so make() receives a valid dipcode instead
		// of the default zero (which is often disallowed by the diplotype mask).
		curr_dipcount = countDiplotypes(Diplotypes[0]);
		currProbs.assign(curr_dipcount, 0.0);
		for (unsigned int tabs = 0; tabs < curr_dipcount; ++tabs) {
			if (tabs >= trans_buf_size) {
				std::fprintf(stderr,
				             "[TRANS_READ_OOB][sampleBackward] sample=%s segment=0 tabs=%u buf_size=%zu curr_dipcount=%u n_segments=%u n_transitions=%u\n",
				             name.c_str(), tabs, trans_buf_size, curr_dipcount, n_segments, n_transitions);
				std::fflush(stderr);
				std::abort();
			}
			sumProbs += (currProbs[tabs] = CurrentTransProbabilities[tabs]);
		}
			next_sampled = sample_rng.sample(currProbs, sumProbs);
		if (next_sampled < 0 || next_sampled >= static_cast<int>(curr_dipcount)) {
			std::fprintf(stderr,
			             "[SAMPLE_IDX_OOB][sampleBackward-single] sample=%s sampled_idx=%d curr_dipcount=%u buf_size=%zu\n",
			             name.c_str(), next_sampled, curr_dipcount, trans_buf_size);
			std::fflush(stderr);
			std::abort();
		}
		makeDiplotypes(Diplotypes[0]);
		DipSampled[0] = curr_dipcodes[next_sampled];
		if (trace_sample) {
			std::fprintf(stderr, "[SAMPLE_DEBUG] sample=%s direction=backward-single n_segments=%u dip_mask=0x%016llx DipSampled0=%u\n",
			             name.c_str(),
			             n_segments,
			             static_cast<unsigned long long>(n_segments ? Diplotypes[0] : 0),
			             n_segments ? static_cast<unsigned>(DipSampled[0]) : 0);
		}
		make(DipSampled, CurrentMissingProbabilities);
		return;
	}

	for (int s = n_segments - 2, toffset = n_transitions ; s >= 0 ; s --) {
		sumProbs = 0.0;
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		toffset -= next_dipcount * curr_dipcount;

		if (next_sampled >= 0) {
			currProbs.assign(curr_dipcount, 0.0);
			for (unsigned int tabs = toffset+next_sampled, trel = 0 ; trel < curr_dipcount && tabs < trans_buf_size ; ++trel, tabs += next_dipcount) {
				if (tabs >= trans_buf_size) {
					std::fprintf(stderr,
					             "[TRANS_READ_OOB][sampleBackward] sample=%s segment=%d tabs=%u buf_size=%zu toffset=%d next_sampled=%d curr_dipcount=%u next_dipcount=%u n_segments=%u n_transitions=%u\n",
					             name.c_str(), s, tabs, trans_buf_size, toffset, next_sampled, curr_dipcount, next_dipcount, n_segments, n_transitions);
					std::fflush(stderr);
					std::abort();
				}
				sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
			}

			next_sampled = sample_rng.sample(currProbs, sumProbs);
			if (next_sampled >= (int)curr_dipcount) {
				std::fprintf(stderr,
				             "[SAMPLE_IDX_OOB][sampleBackward] sample=%s segment=%d sampled_idx=%d curr_dipcount=%u next_dipcount=%u toffset=%d trans_buf_size=%zu\n",
				             name.c_str(), s, next_sampled, curr_dipcount, next_dipcount, toffset, trans_buf_size);
				std::fflush(stderr);
				std::abort();
			}

			// SAMPLE TRACE: Log sampling details before drawing
			if (supersite_trace_enabled() && s == 0) {
				std::fprintf(stderr, "[SAMPLE_BWD_A] sample=%s segment=%d sumProbs=%.15f\n",
				             name.c_str(), s, sumProbs);
				double cumulative = 0.0;
				for (unsigned int i = 0; i < curr_dipcount; i++) {
					cumulative += currProbs[i];
					std::fprintf(stderr, "  dip[%u] prob=%.15f norm=%.15f cum=%.15f\n",
					             i, currProbs[i], currProbs[i]/sumProbs, cumulative/sumProbs);
				}
			}

			makeDiplotypes(Diplotypes[s]);
			unsigned char selected_dipcode = curr_dipcodes[next_sampled];

			// SAMPLE TRACE: Log selected index after drawing (all segments when tracing)
			if (trace_sample) {
				unsigned char hap0 = (selected_dipcode >> 3);
				unsigned char hap1 = (selected_dipcode & 7);
				std::fprintf(stderr, "[SAMPLE_BWD_PICK] sample=%s seg=%d selected_idx=%d dipcode=%u hap0=%u hap1=%u\n",
				             name.c_str(), s, next_sampled, static_cast<unsigned>(selected_dipcode),
				             static_cast<unsigned>(hap0), static_cast<unsigned>(hap1));
			}

			DipSampled[s] = selected_dipcode;
		} else {
			const unsigned int n_joint = curr_dipcount * next_dipcount;
			currProbs.assign(n_joint, 0.0);
			for (unsigned int tabs = toffset, trel = 0 ; trel < n_joint && tabs < n_transitions ; ++trel, ++tabs) {
				if (tabs >= trans_buf_size) {
					std::fprintf(stderr,
					             "[TRANS_READ_OOB][sampleBackward] sample=%s segment=%d tabs=%u buf_size=%zu toffset=%d curr_dipcount=%u next_dipcount=%u n_segments=%u n_transitions=%u\n",
					             name.c_str(), s, tabs, trans_buf_size, toffset, curr_dipcount, next_dipcount, n_segments, n_transitions);
					std::fflush(stderr);
					std::abort();
				}
				sumProbs += (currProbs[trel] = CurrentTransProbabilities[tabs]);
			}

			// SAMPLE TRACE: Log sampling details before drawing (initial backward step)
			if (supersite_trace_enabled()) {
				std::fprintf(stderr, "[SAMPLE_BWD_INIT] sample=%s segment=%d sumProbs=%.15f\n",
				             name.c_str(), s, sumProbs);
				double cumulative = 0.0;
				unsigned int n_probs = (n_transitions - toffset);
				for (unsigned int i = 0; i < n_probs && i < 64; i++) {
					cumulative += currProbs[i];
					std::fprintf(stderr, "  joint[%u] prob=%.15f norm=%.15f cum=%.15f\n",
					             i, currProbs[i], currProbs[i]/sumProbs, cumulative/sumProbs);
				}
			}

				next_sampled = sample_rng.sample(currProbs, sumProbs);
			if (next_sampled >= static_cast<int>(n_joint)) {
				std::fprintf(stderr,
				             "[SAMPLE_IDX_OOB][sampleBackward-init] sample=%s segment=%d sampled_idx=%d n_joint=%u curr_dipcount=%u next_dipcount=%u toffset=%d trans_buf_size=%zu\n",
				             name.c_str(), s, next_sampled, n_joint, curr_dipcount, next_dipcount, toffset, trans_buf_size);
				std::fflush(stderr);
				std::abort();
			}

			// SAMPLE TRACE: Log selected index after drawing
			if (supersite_trace_enabled()) {
				std::fprintf(stderr, "  selected_joint_idx=%d -> seg[%d]_idx=%d seg[%d]_idx=%d\n",
				             next_sampled, s+1, next_sampled % next_dipcount, s, next_sampled / next_dipcount);
			}

			makeDiplotypes(Diplotypes[s+1]);
			DipSampled[s+1] = curr_dipcodes[next_sampled % next_dipcount];
			makeDiplotypes(Diplotypes[s]);
			next_sampled = next_sampled / next_dipcount;
			if (next_sampled >= (int)curr_dipcount) {
				std::fprintf(stderr,
				             "[SAMPLE_IDX_OOB][sampleBackward] sample=%s segment=%d sampled_idx=%d curr_dipcount=%u next_dipcount=%u toffset=%d trans_buf_size=%zu (post-divide)\n",
				             name.c_str(), s, next_sampled, curr_dipcount, next_dipcount, toffset, trans_buf_size);
				std::fflush(stderr);
				std::abort();
			}
			DipSampled[s] = curr_dipcodes[next_sampled];
			if (!DIP_GET(Diplotypes[s+1], DipSampled[s+1]) || !DIP_GET(Diplotypes[s], DipSampled[s])) {
				std::fprintf(stderr,
				             "[SAMPLE_DIPC_INVALID][backward] sample=%s seg=%d dip_s=%u dip_s1=%u mask_s=0x%016llx mask_s1=0x%016llx\n",
				             name.c_str(), s, static_cast<unsigned>(DipSampled[s]), static_cast<unsigned>(DipSampled[s+1]),
				             static_cast<unsigned long long>(Diplotypes[s]), static_cast<unsigned long long>(Diplotypes[s+1]));
				std::fflush(stderr);
				std::abort();
			}
		}
		next_dipcount = curr_dipcount;
	}
	if (trace_sample) {
		std::fprintf(stderr, "[SAMPLE_DEBUG] sample=%s direction=backward n_segments=%u dip_mask=0x%016llx DipSampled0=%u\n",
		             name.c_str(),
		             n_segments,
		             static_cast<unsigned long long>(n_segments ? Diplotypes[0] : 0),
		             n_segments ? static_cast<unsigned>(DipSampled[0]) : 0);
	}
	make(DipSampled, CurrentMissingProbabilities);

	// HYPOTHESIS 2 DEBUGGING
	if (super_sites && !debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME) {
		for (const auto& ss : *super_sites) {
			if ((int)ss.global_site_id == debug::SUPERDEBUG_BP) {
				debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After make() in sampleBackward");
				break;
			}
		}
	}

	// CRITICAL: Post-HMM supersite projection (same as in sampleForward)
	// See detailed comment in sampleForward() above - this is where projection happens.
	if (super_sites) projectSupersites();

	// HYPOTHESIS 2 DEBUGGING
	if (super_sites && !debug::SUPERDEBUG_SAMPLENAME.empty() && this->name == debug::SUPERDEBUG_SAMPLENAME) {
		for (const auto& ss : *super_sites) {
			if ((int)ss.global_site_id == debug::SUPERDEBUG_BP) {
				debug::print_supersite_state(this, ss, *super_site_var_index, "Hypo2: After projectSupersites() in sampleBackward");
				break;
			}
		}
	}
}

void genotype::solve() {
	unsigned int curr_dipcount = 0, prev_dipcount = 1;
	vector < vector < double > > maxProbs = vector < vector < double > > (n_segments, vector < double > ());
	vector < vector < int > > maxIndexes = vector < vector < int > > (n_segments, vector < int > ());

	for (int s = 0, toffset = 0, trel = 0 ; s < n_segments ; s ++) {
		curr_dipcount = countDiplotypes(Diplotypes[s]);
		maxProbs[s] = vector < double > (curr_dipcount, 0.0);
		maxIndexes[s] = vector < int > (curr_dipcount, 0);
		for (int t = 0 ; t < prev_dipcount * curr_dipcount ; t++) {
			int prev_dip = t/curr_dipcount;
			int next_dip = t%curr_dipcount;
			//double currProb = (s?maxProbs[s-1][prev_dip]:1.0) * StoredProbs[t+toffset];
			double currProb = (s?maxProbs[s-1][prev_dip]:1.0) * (ProbMask[t+toffset]?ProbStored[trel++]:1e-6);
			if (currProb > maxProbs[s][next_dip]) {
				maxProbs[s][next_dip] = currProb;
				maxIndexes[s][next_dip] = prev_dip;
			}
		}
		double sumProb = 0.0;
		for (int d = 0 ; d < curr_dipcount ; d ++) sumProb += maxProbs[s][d];
		for (int d = 0 ; d < curr_dipcount ; d ++) maxProbs[s][d] /= sumProb;
		toffset += prev_dipcount * curr_dipcount;
		prev_dipcount = curr_dipcount;
	}

	vector < unsigned char > DipSampled = vector < unsigned char >(n_segments, 0);
	unsigned int bestDip = alg.imax(maxProbs.back());
	makeDiplotypes(Diplotypes.back());
	DipSampled.back() = curr_dipcodes[bestDip];
	for (int s = DipSampled.size() - 2 ; s >= 0 ; s --) {
		bestDip = maxIndexes[s+1][bestDip];
		makeDiplotypes(Diplotypes[s]);
		DipSampled[s] = curr_dipcodes[bestDip];
	}
	make(DipSampled);
}

void genotype::store(vector < double > & CurrentTransProbabilities, vector < float > & CurrentMissingProbabilities) {
	if (ProbMask.size() == 0) {
		n_stored_transitionProbs = 0;
		ProbMask = vector < bool > (n_transitions, false);
		for (unsigned int t = 0 ; t < n_transitions ; t ++) if (CurrentTransProbabilities[t] >= 1e-6) {
			n_stored_transitionProbs ++;
			ProbMask[t] = true;
		}
		ProbStored = vector  < float > (n_stored_transitionProbs, 0.0);
		ProbMissing = vector < float > (n_missing * HAP_NUMBER, 0.0);
	}
	for (unsigned int t = 0, trel = 0 ; t < n_transitions ; t ++) {
		if (ProbMask[t]) ProbStored[trel++] += CurrentTransProbabilities[t];
	}
	for (unsigned int m = 0 ; m < (n_missing * HAP_NUMBER) ; m ++) ProbMissing[m] += CurrentMissingProbabilities[m];
	n_storage_events ++;
}

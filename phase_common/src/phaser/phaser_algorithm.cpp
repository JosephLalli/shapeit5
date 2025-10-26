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

#include <phaser/phaser_header.h>

using namespace std;

void * phaseWindow_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_worker, id_job;
	pthread_mutex_lock(&S->mutex_workers);
	id_worker = S->i_workers ++;
	pthread_mutex_unlock(&S->mutex_workers);
	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_ind) vrb.progress("  * HMM computations", id_job*1.0/S->G.n_ind);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_ind) S->phaseWindow(id_worker, id_job);
		else pthread_exit(NULL);
	}
}

void phaser::phaseWindow(int id_worker, int id_job) {
	threadData[id_worker].make(id_job, options["hmm-window"].as < double > ());

	//HMM compute in windows
	for (int w = 0 ; w < threadData[id_worker].size() ; w ++) {
		if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
		statH.push(threadData[id_worker].Kstates[w].size()*1.0);
		statS.push(threadData[id_worker].Windows.W[w].lengthBP(V) * 1.0e-6);
		if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);

		int outcome = 0;

		if (G.vecG[id_job]->double_precision) {
			//Run using double precision as underflow happened previously
			haplotype_segment_double HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M);
			HS.forward();
			outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M);
		} else {
			//Try single precision as this is faster
			haplotype_segment_single HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M);
			HS.forward();
			outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M);

			//Underflow happening with single precision, rerun using double precision
			if (outcome != 0) {
				haplotype_segment_double HS(G.vecG[id_job], H.H_opt_hap, threadData[id_worker].Kstates[w], threadData[id_worker].Windows.W[w], M);
				HS.forward();
				outcome = HS.backward(threadData[id_worker].T, threadData[id_worker].M);
				G.vecG[id_job]->double_precision = true;
				n_underflow_recovered_precision++;
			}
		}

		//
		switch (outcome) {
		case -2: vrb.error("Diploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		case -1: vrb.error("Haploid underflow impossible to recover for [" + G.vecG[id_job]->name + "]");
		}
		n_underflow_recovered_summing += outcome;
	}

	//Copy over new IBD2 constraints into H
	if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
	H.Kbanned.pushIBD2(id_job, threadData[id_worker].Kbanned);
	if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);

	//Sampling / Merging / Storing
	vector < bool > flagMerges;
	switch (iteration_types[iteration_stage]) {
	case STAGE_BURN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						break;
	case STAGE_PRUN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						G.vecG[id_job]->mapMerges(threadData[id_worker].T, options["mcmc-prune"].as < double > (), flagMerges);
						G.vecG[id_job]->performMerges(threadData[id_worker].T, flagMerges);
						break;
	case STAGE_MAIN:	G.vecG[id_job]->sample(threadData[id_worker].T, threadData[id_worker].M);
						G.vecG[id_job]->store(threadData[id_worker].T, threadData[id_worker].M);
						break;
	}

	// UNIFIED: Per-sample multiallelic enforcement for ALL modes immediately after sampling
	// This ensures violations are corrected while donor context is available and before
	// any downstream operations (IBD2 collapse, H update) that could be affected
	if (oneallele_enforcer.enabled() && multiallelic_map.size() > 0) {
		// Reset per-sample epoch stats before enforcement to prevent cumulative counting
		oneallele_enforcer.reset_sample_epoch_stats();
		
		// Enforce constraints immediately after sampling
		// - MICRO-DONOR: Uses threadData[id_worker].Kstates for donor-weighted scoring
		// - MICRO: Uses threadData[id_worker].Kstates if available, otherwise donor-agnostic
		// - TRANSITION: Kstates available but not used (transition-only scoring)
		oneallele_enforcer.enforce_sample(multiallelic_map, *G.vecG[id_job], V, threadData[id_worker].Kstates,
										  current_iteration_context, id_job);
		
		// Get stats from this sample's enforcement
		shapeit5::modules::OneAlleleEpochStats sample_stats = oneallele_enforcer.sample_epoch_stats();
		
		// Thread-safe accumulation into global epoch stats
		if (options["thread"].as < int > () > 1) pthread_mutex_lock(&mutex_workers);
		oneallele_enforcer.accumulate_sample_stats(sample_stats);
		if (options["thread"].as < int > () > 1) pthread_mutex_unlock(&mutex_workers);
	}
}

void phaser::phaseWindow() {
	tac.clock();
	int n_thread = options["thread"].as < int > ();
	n_underflow_recovered_summing = 0;
	n_underflow_recovered_precision = 0;
	i_workers = 0; i_jobs = 0;
	
	statH.clear(); statS.clear();
	storedKsizes.clear();
	if (n_thread > 1) {
		for (int t = 0 ; t < n_thread ; t++) pthread_create( &id_workers[t] , NULL, phaseWindow_callback, static_cast<void *>(this));
		for (int t = 0 ; t < n_thread ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_ind ; i ++) {
		phaseWindow(0, i);
		vrb.progress("  * HMM computations", (i+1)*1.0/G.n_ind);
	}
	vrb.bullet("HMM computations [K=" + stb.str(statH.mean(), 1) + "+/-" + stb.str(statH.sd(), 1) + " / W=" + stb.str(statS.mean(), 2) + "Mb / US=" + stb.str(n_underflow_recovered_summing) + " / UP=" + stb.str(n_underflow_recovered_precision) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
}

std::string phaser::get_iteration_string(int iter, int stage) {
	std::string stage_name;
	switch (iteration_types[stage]) {
		case STAGE_BURN: stage_name = "burn-in"; break;
		case STAGE_PRUN: stage_name = "pruning"; break;
		case STAGE_MAIN: stage_name = "main"; break;
		default: stage_name = "unknown"; break;
	}
	return stage_name + "-" + stb.str(iter+1) + "/" + stb.str(iteration_counts[stage]);
}

void phaser::phase() {
	unsigned long n_old_segments = G.numberOfSegments(), n_new_segments = 0, current_iteration = 0;
	for (iteration_stage = 0 ; iteration_stage < iteration_counts.size() ; iteration_stage ++) {
		for (int iter = 0 ; iter < iteration_counts[iteration_stage] ; iter ++) {
			// Store current iteration context for debug logging
			current_iteration_context = get_iteration_string(iter, iteration_stage);
			
			//VERBOSE
			switch (iteration_types[iteration_stage]) {
			case STAGE_BURN:	vrb.title("Burn-in iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_PRUN:	vrb.title("Pruning iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			case STAGE_MAIN:	vrb.title("Main iteration [" + stb.str(iter+1) + "/" + stb.str(iteration_counts[iteration_stage]) + "]"); break;
			}
			
			// Reset epoch statistics for this iteration
			if (oneallele_enforcer.enabled()) {
				oneallele_enforcer.reset_epoch_stats();
			}
			
			//SELECT NEW STATES WITH PBWT
			H.select();
			//PHASE DATA
			phaseWindow();
			//MERGE IBD2 PAIRS
			H.Kbanned.collapse();
			
			// Report accumulated per-sample enforcement statistics
			// Note: Enforcement already happened per-sample in phaseWindow(), we just report here
			if (oneallele_enforcer.enabled() && multiallelic_map.size() > 0) {
				const auto& epoch_stats = oneallele_enforcer.epoch_stats();
				if (epoch_stats.violations_found > 0 || epoch_stats.flips_applied > 0) {
					std::string mode_str;
					switch (oneallele_enforcer.mode()) {
						case shapeit5::modules::OneAlleleMode::TRANSITION: mode_str = "transition"; break;
						case shapeit5::modules::OneAlleleMode::MICRO: mode_str = "micro"; break;
						case shapeit5::modules::OneAlleleMode::MICRO_DONOR: mode_str = "micro-donor"; break;
					}
					
					std::string base_stats = "[violations=" + stb.str(epoch_stats.violations_found) + 
											" / flipped=" + stb.str(epoch_stats.flips_applied) + "]";
					
					// Add micro-donor specific details if available
					if (oneallele_enforcer.mode() == shapeit5::modules::OneAlleleMode::MICRO_DONOR && 
						(epoch_stats.emission_dominated_decisions > 0 || epoch_stats.transition_dominated_decisions > 0)) {
						base_stats += " [emission=" + stb.str(epoch_stats.emission_dominated_decisions) + 
									 " / transition=" + stb.str(epoch_stats.transition_dominated_decisions) + 
									 " / genotype-changes=" + stb.str(epoch_stats.genotype_changes) + "]";
					}
					
					vrb.bullet("Multiallelic correction (" + mode_str + ") " + base_stats);
				}
			}
			
			//UPDATE H with new sampled haplotypes
			H.updateHaplotypes(G);
			//if (options.count("pedigree")) H.checkScaffoldPedigrees(G, options["pedigree"].as < string > ());
			//TRANSPOSE H from Hfirst to Vfirst (for next PBWT compute)
			H.transposeHaplotypes_H2V(false);
			//UPDATE PS after prunning
			if (iteration_types[iteration_stage] == STAGE_PRUN) {
				n_new_segments = G.numberOfSegments();
				vrb.bullet("Trimming [pc=" + stb.str((1-n_new_segments*1.0/n_old_segments)*100, 2) + "%]");
			}
		}
	}
}

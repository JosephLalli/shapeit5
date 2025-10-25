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
#include <unordered_map>
#include <string>
#include <cmath>

using namespace std;

void * hmmcompute_callback(void * ptr) {
	phaser * S = static_cast< phaser * >( ptr );
	int id_job, id_thread;

	pthread_mutex_lock(&S->mutex_workers);
	id_thread = S->i_threads ++;
	pthread_mutex_unlock(&S->mutex_workers);

	for(;;) {
		pthread_mutex_lock(&S->mutex_workers);
		id_job = S->i_jobs ++;
		if (id_job <= S->G.n_samples) vrb.progress("  * Processing", (id_job+1)*1.0/S->G.n_samples);
		pthread_mutex_unlock(&S->mutex_workers);
		if (id_job < S->G.n_samples) S->hmmcompute(id_job, id_thread);
		else pthread_exit(NULL);
	}
}

void phaser::hmmcompute(int id_job, int id_thread) {
	//Mapping storage events
	vector < vector < unsigned int > > cevents;
	G.mapUnphasedOntoScaffold(id_job, cevents);

	//Viterbi paths
	vector < int > path0, path1;

	//Forward-Backward-Viterbi passes for hap0
	thread_hmms[id_thread]->setup(2*id_job+0);
	thread_hmms[id_thread]->viterbi(path0);
	double pf0 = thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, path0);


	//Forward-Backward-Viterbi passes for hap1
	thread_hmms[id_thread]->setup(2*id_job+1);
	thread_hmms[id_thread]->viterbi(path1);
	double pf1 = thread_hmms[id_thread]->forward();
	thread_hmms[id_thread]->backward(cevents, path1);

	//Phase remaining unphased using viterbi [singletons, etc ...]
	G.phaseCoalescentViterbi(id_job, path0, path1, M, options.count("score-singletons"));
}

void phaser::phase() {
	//STEP1: haplotype selection
	vrb.title("PBWT pass");
	H.initialize(V,	options["pbwt-modulo"].as < double > (),
			options["pbwt-mdr"].as < double > (),
			options["pbwt-depth-common"].as < int > (),
			options["pbwt-depth-rare"].as < int > (),
			options["pbwt-mac"].as < int > ());
	H.scanIBD2(V);
	H.select(V, G);

	//STEP2: HMM computations
	vrb.title("HMM computations");
	thread_hmms = vector < hmm_scaffold * > (nthreads);
	for(int t = 0; t < nthreads ; t ++) thread_hmms[t] = new hmm_scaffold(V, G, H, M);
	if (nthreads > 1) {
		i_jobs = i_threads = 0;
		for (int t = 0 ; t < nthreads ; t++) pthread_create( &id_workers[t] , NULL, hmmcompute_callback, static_cast<void *>(this));
		for (int t = 0 ; t < nthreads ; t++) pthread_join( id_workers[t] , NULL);
	} else for (int i = 0 ; i < G.n_samples ; i ++) {
		hmmcompute(i, 0);
		vrb.progress("  * Processing", (i+1)*1.0/G.n_samples);
	}
	for(int t = 0; t < nthreads ; t ++) delete thread_hmms[t];
	vrb.bullet("Processing (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");

	//STEP3: MERGE BACK ALL TOGETHER
    G.merge_by_transpose_I2V();

    //STEP3.5: Enforce one-allele-per-haplotype at multiallelic rare positions
    if (enforce_oneallele_rare) {
        // Build groups of rare variants by (chr,bp)
        std::unordered_map<long long, std::vector<int> > pos_groups; // key: (chrom_id<<32)|bp
        pos_groups.reserve(V.sizeRare());
        // Build mapping chrom name to small id
        std::unordered_map<std::string,int> chr2id;
        int next_chr_id = 0;
        for (int vt = 0; vt < V.sizeFull(); ++vt) {
            if (V.vec_full[vt]->type != VARTYPE_RARE) continue;
            const std::string &chr = V.vec_full[vt]->chr;
            auto it = chr2id.find(chr);
            if (it == chr2id.end()) it = chr2id.emplace(chr, next_chr_id++).first;
            long long key = (static_cast<long long>(it->second) << 32) | static_cast<unsigned int>(V.vec_full[vt]->bp);
            pos_groups[key].push_back(V.vec_full[vt]->idx_rare);
        }

        oneallele_rare_stats.positions_checked = pos_groups.size();

        auto clamp01 = [](double x){ return std::max(1e-6, std::min(1.0-1e-6, x)); };

        // For each multiallelic site
        for (const auto &kv : pos_groups) {
            const std::vector<int> &vrs = kv.second;
            if (vrs.size() < 2) continue;

            // Bucket entries per sample
            std::vector< std::vector< sparse_genotype* > > buckets(G.n_samples);
            for (int vr : vrs) {
                auto &entries = G.GRvar_genotypes[vr];
                for (auto &e : entries) {
                    buckets[e.idx].push_back(&e);
                }
            }

            // Per-sample resolution
            for (int s = 0; s < (int)G.n_samples; ++s) {
                auto &vec = buckets[s];
                if (vec.size() < 2) continue;

                // Count ALT per hap and collect phased het contributors
                int alt_h0 = 0, alt_h1 = 0;
                for (auto *e : vec) {
                    if (e->het && e->pha && !e->mis) { alt_h0 += e->al0; alt_h1 += e->al1; }
                }
                bool violation = (alt_h0 > 1) || (alt_h1 > 1);
                if (!violation) continue;
                oneallele_rare_stats.sample_violations_found++;
                
                // Count extreme violations (>2 ALT alleles at same position in same sample)
                int total_alts = alt_h0 + alt_h1;
                if (total_alts > 2) {
                    oneallele_rare_stats.extreme_violations_found++;
                }

                // Resolve by flipping lowest-PP among ALT on offending hap; may require >1 flip if >2 contributors
                for (int iter = 0; iter < 4 && ((alt_h0 > 1) || (alt_h1 > 1)); ++iter) {
                    int offending = (alt_h0 > 1) ? 0 : 1;
                    // Collect candidates on offending hap
                    sparse_genotype *anchor = nullptr;
                    sparse_genotype *to_flip = nullptr;
                    double best_anchor_pp = -1.0;
                    double worst_pp = 2.0;
                    for (auto *e : vec) {
                        if (!(e->het && e->pha && !e->mis)) continue;
                        int hap_alt = offending == 0 ? e->al0 : e->al1;
                        if (hap_alt == 1) {
                            double pp = clamp01(e->prob);
                            // track anchor (highest PP) and flip target (lowest PP)
                            if (pp > best_anchor_pp) { best_anchor_pp = pp; anchor = e; }
                            if (pp < worst_pp) { worst_pp = pp; to_flip = e; }
                        }
                    }
                    if (!to_flip || !anchor || to_flip == anchor) break;

                    // Compute joint PP if aligning to anchor
                    double p1 = clamp01(to_flip->prob);
                    double p2 = clamp01(anchor->prob);
                    // Numerically stable: newPP = sigmoid(logit(p2) - logit(p1))
                    auto logit = [](double x){ return std::log(x) - std::log1p(-x); };
                    double delta = logit(p1) - logit(p2);
                    double newPP = 1.0 / (1.0 + std::exp(delta));
                    newPP = clamp01(newPP);

                    // Flip the selected entry
                    {
                        unsigned int t_al0 = to_flip->al0;
                        to_flip->al0 = to_flip->al1;
                        to_flip->al1 = t_al0;
                    }
                    to_flip->prob = newPP;
                    anchor->prob = newPP; // set both to joint confidence
                    oneallele_rare_stats.flips_applied++;

                    // Recompute counts quickly for this pair
                    alt_h0 = alt_h1 = 0;
                    for (auto *e : vec) {
                        if (e->het && e->pha && !e->mis) { alt_h0 += e->al0; alt_h1 += e->al1; }
                    }
                }
            }
        }
    }

    //VERBOSE
    vrb.bullet("Solving summary:");
	G.nmiss_total = G.nmiss_imputation + G.nmiss_families + G.nmiss_monomorphic;
	G.nhets_total = G.nhets_families + G.nhets_imputation + G.nhets_coalescent;
	vrb.bullet2("#Hets=" + stb.str(G.nhets_total) + " / Phased by HMM=" + stb.str(G.nhets_imputation) + ", PED=" + stb.str(G.nhets_families) + ", SING=" + stb.str(G.nhets_coalescent));
	vrb.bullet2("#Miss=" + stb.str(G.nmiss_total) + " / Imputed by HMM=" + stb.str(G.nmiss_imputation) + ", PED=" + stb.str(G.nmiss_families) + ", MONO=" + stb.str(G.nmiss_monomorphic));

}

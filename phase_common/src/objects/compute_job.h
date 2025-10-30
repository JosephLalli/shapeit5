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

#ifndef _COMPUTE_THREAD_H
#define _COMPUTE_THREAD_H

#include <utils/otools.h>

#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/genotype_set.h>
#include <containers/variant_map.h>
#include <containers/window_set.h>
#include <models/super_site_accessor.h>  // Phase 3: For SuperSite struct

class compute_job {
public:

	//DATA
	variant_map & V;
	genotype_set & G;
	conditioning_set & H;
	
	// Phase 3: Supersite references (set by phaser, read-only)
	const std::vector<SuperSite>* super_sites;
	const std::vector<int>* locus_to_super_idx;
	const std::vector<int>* super_site_var_index;

	//Probabilities
	std::vector < double > T;
	std::vector < float > M;
	
	// Phase 3: Supersite multi-class posteriors (per window)
	std::vector < float > SC;  // CurrentSuperClassPosteriors: layout [ss0: HAP_NUMBER*C, ss1: HAP_NUMBER*C, ...]
	std::vector < bool > anchor_has_missing;  // Per-supersite flag: true if all members missing for this sample

	//Windows
	window_set Windows;

	//States
	std::vector < track > Kbanned;
	std::vector < std::vector < unsigned int > > Kstates;

	//Random states
	std::vector < unsigned int > Ordering;
	int Oiterator;

	compute_job(variant_map & , genotype_set & , conditioning_set & , unsigned int n_max_transitions , unsigned int n_max_missing,
	           const std::vector<SuperSite>* ss = nullptr,
	           const std::vector<int>* locus_ss_idx = nullptr,
	           const std::vector<int>* ss_var_idx = nullptr);
	~compute_job();

	void free();
	void make(unsigned int, double);
	unsigned int size();
};

inline
unsigned int compute_job::size() {
	 return Windows.size();
}

#endif

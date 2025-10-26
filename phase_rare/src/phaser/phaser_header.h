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

#ifndef _PHASER_H
#define _PHASER_H

#include <utils/otools.h>
#include <objects/hmm_parameters.h>

#include <containers/genotype_set/genotype_set_header.h>
#include <containers/state_set.h>
#include <containers/conditioning_set/conditioning_set_header.h>
#include <containers/variant_map.h>

#include <models/hmm_scaffold/hmm_scaffold_header.h>


class phaser {
public:
	//COMMAND LINE OPTIONS
	bpo::options_description descriptions;
	bpo::variables_map options;

	//INTERNAL DATA
	conditioning_set H;
	genotype_set G;
	hmm_parameters M;
    variant_map V;
    state_set P;

    // One-allele constraint (rare)
    enum class OneAlleleRareMode {
        PP_BASIC,           // Basic PP-based enforcement (existing)
        PP_ENHANCED,        // PP-based with Li-Stephens enhancement
        SPARSE_TRANSITION,  // Sparse transition scoring with donor context
        SPARSE_MICRO        // Sparse enumeration for complex cases
    };
    bool enforce_oneallele_rare;
    OneAlleleRareMode oneallele_rare_mode;
    std::string oneallele_rare_stats_path;
    struct OneAlleleRareStats {
        // Base metrics (existing)
        unsigned long long positions_checked = 0;
        unsigned long long sample_violations_found = 0;
        unsigned long long flips_applied = 0;
        unsigned long long extreme_violations_found = 0; // >2 ALT alleles at same position in same sample
        
        // Enhanced metrics (new - for improved rare variant enforcement)
        unsigned long long sparse_donor_resolutions = 0;    // Resolutions using PBWT donor context
        unsigned long long pp_only_resolutions = 0;         // Resolutions using PP-only scoring
        unsigned long long genotype_changes = 0;            // REF<->ALT changes (not just phase)
        unsigned long long phase_only_changes = 0;          // Only haplotype assignment changes
        unsigned long long complex_enumeration_cases = 0;   // Cases requiring full enumeration
        unsigned long long li_stephens_enhanced = 0;        // PP resolutions enhanced with Li-Stephens
    } oneallele_rare_stats;

	//MULTI-THREADING
	int i_jobs, i_threads, nthreads;
	std::vector < pthread_t > id_workers;
	pthread_mutex_t mutex_workers;
	std::vector < std::vector < std::pair < int, float > > > thread_data;
	std::vector < hmm_scaffold * > thread_hmms;

	//STATS
	int totalSite, doneSite;
	unsigned long int n_common_yphased;
	unsigned long int n_common_nphased;
	unsigned long int n_rare_yphased;
	unsigned long int n_rare_nphased;
	stats1D statCS;

	//GENOMIC REGION
	std::string chrid;
	int input_start;
	int input_stop;
	int scaffold_start;
	int scaffold_stop;
	std::string input_gregion;
	std::string scaffold_gregion;

	//CONSTRUCTOR
	phaser();
	~phaser();

	//METHODS
	void hmmcompute(int, int);
	void phase();


	//PARAMETERS
	void buildCoordinates();
	void declare_options();
	void parse_command_line(std::vector < std::string > &);
	void check_options();
	void verbose_options();
	void verbose_files();

	//
	void read_files_and_initialise();
	void phase(std::vector < std::string > &);
	void write_files_and_finalise();
};


#endif


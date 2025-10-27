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

#ifndef _SNP_SET_H
#define _SNP_SET_H

#include <utils/otools.h>
#include <objects/variant.h>
#include <io/gmap_reader.h>
#include <containers/bitmatrix.h>

#include <cstdint>

struct supersite_summary {
	uint32_t total_sites = 0;
	uint32_t total_super_sites = 0;
	uint32_t total_variants = 0;
	uint32_t collapsed_variants = 0;
	uint32_t haplotype_conflicts = 0;
};

struct supersite_desc {
	uint32_t first_variant_index = 0;
	uint16_t variant_span = 0;
	uint16_t n_alt = 0;
	uint8_t bitwidth = 0;
	uint32_t panel_offset = 0;
	uint32_t codes_count = 0;
	uint32_t conflicting_haps = 0;
	bool is_super_site = false;
};

class variant_map {
public :
	//DATA
	std::vector < variant * > vec_pos;			//vector of variants ordered by position in bp
	std::multimap < int, variant * > map_pos;	//associative container of variant with position in bp

	//SUPERSITE DATA
	std::vector<supersite_desc> supersites;
	std::vector<uint32_t> variant_to_site;
	std::vector<uint8_t> variant_alt_code;
	std::vector<uint8_t> variant_is_anchor;
	std::vector<uint32_t> supersite_alt_variant_index;
	std::vector<uint32_t> supersite_alt_variant_offset;
	std::vector<uint8_t> supersite_codes;

	// SUPERSITE POSTERIOR OFFSETS (v2 multi-code impute)
	std::vector<uint32_t> supersite_posterior_offset; // per supersite, start offset in ProbMissingMulti
	uint32_t total_posterior_size = 0;               // sum_s ( (n_alt+1) * HAP_NUMBER )

	//CONSTRUCTOR/DESTRUCTOR
	variant_map();
	~variant_map();

	//METHODS
	int size();
	std::vector < variant * > getByPos(int);
	//vector < variant * > getByRef(int, string &, string &);
	//variant * getByIndex(int);
	void push(variant *);
	void setGeneticMap(gmap_reader&);
	void setGeneticMap();
	int setCentiMorgan(std::vector < int > & pos_bp, std::vector < double > & pos_cM);
	int interpolateCentiMorgan(std::vector < int > & pos_bp, std::vector < double > & pos_cM);
	unsigned int length();
	double lengthcM();
	supersite_summary buildSupersites(const bitmatrix & hap_matrix, unsigned long n_hap);
 
 	// compute posterior offsets for v2 per-code accumulators
	void buildPosteriorOffsets();
};

#endif

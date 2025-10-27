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

#include <containers/variant_map.h>

#include <algorithm>
#include <limits>

using std::vector;
using std::pair;
using std::multimap;
using std::sort;

namespace {

inline uint8_t required_bitwidth(uint16_t n_alt) {
	uint32_t states = static_cast<uint32_t>(n_alt) + 1u;
	uint8_t bits = 0u;
	while ((1u << bits) < states) ++bits;
	return bits == 0u ? 1u : bits;
}

} // namespace

variant_map::variant_map() {
}

variant_map::~variant_map() {
	for (int s = 0 ; s < vec_pos.size() ; s++) delete vec_pos[s];
	vec_pos.clear();
	map_pos.clear();
}

int variant_map::size() {
	return vec_pos.size();
}
/*
variant * variant_map::getByIndex (int i) {
	return vec_pos[i];
}
*/
vector < variant * > variant_map::getByPos (int pos) {
	vector < variant * > vecS;
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) vecS.push_back(it->second);
	return vecS;
}
/*
vector < variant * > variant_map::getByRef(int pos, string & ref, string & alt) {
	vector < variant * > vecS = vector < variant * >();
	pair < multimap < int , variant * >::iterator , multimap < int , variant * >::iterator > ret = map_pos.equal_range(pos);
	for (multimap < int , variant * >::iterator it = ret.first ; it != ret.second ; ++it) {
		if (it->second->ref == ref && it->second->alt == alt) vecS.push_back(it->second);
	}
	return vecS;
}
*/
void variant_map::push(variant * v) {
	vec_pos.push_back(v);
	map_pos.insert(pair < int , variant * > (v->bp, v));
}

int variant_map::setCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int cpt = 0;
	for (int l = 0 ; l < pos_cM.size() ; l ++) {
		vector  < variant * > vecS = getByPos(pos_bp[l]);
		for (int si = 0 ; si < vecS.size() ; si ++) {
			vecS[si]->cm = pos_cM[l];
			cpt++;
		}
	}
	return cpt;
}
/*
int variant_map::interpolateCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int cpt = 0;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);
	for (int s = 0 ; s < vec_pos.size() ; s ++) {
		if (vec_pos[s]->cm < 0) {
			if (vec_pos[s]->bp < pos_bp[0]) vec_pos[s]->cm = pos_cM[0] - mean_rate * (pos_bp[0] - vec_pos[s]->bp);
			else if (vec_pos[s]->bp > pos_bp.back()) vec_pos[s]->cm = pos_cM.back() + mean_rate * (vec_pos[s]->bp - pos_bp.back());
			else {
				int index_from, index_to;
				for ( index_from = 0 ; index_from < pos_bp.size() && pos_bp[index_from] < vec_pos[s]->bp; ) index_from ++ ;
				for ( index_to = pos_cM.size() - 1 ; index_to >= 0 && pos_bp[index_to] > vec_pos[s]->bp; ) index_to --;
				index_from--;
				index_to++;
				vec_pos[s]->cm = pos_cM[index_from] + (vec_pos[s]->bp - pos_bp[index_from]) * (pos_cM[index_to] - pos_cM[index_from]) / (pos_bp[index_to] - pos_bp[index_from]);
				}
			cpt++;
		}
		vrb.progress("  * cM interpolation", (s+1)*1.0/vec_pos.size());
	}
	return cpt;
}
*/
int variant_map::interpolateCentiMorgan(vector < int > & pos_bp, vector < double > & pos_cM) {
	int n_interpolated = 0, i_locus = 0;
	double base, rate, dist;
	double mean_rate = (pos_cM.back() - pos_cM[0]) / (pos_bp.back() - pos_bp[0]);

	//Set up first positions to be mean rate
	while (i_locus<vec_pos.size() && vec_pos[i_locus]->bp < pos_bp[0]) {
		base = pos_cM[0];
		dist = (pos_bp[0] - vec_pos[i_locus]->bp);
		vec_pos[i_locus]->cm = base - mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}

	//Set up middle positions using interpolation
	int closest_pos = 1;
	for (; i_locus < vec_pos.size() ; ) {
		if (vec_pos[i_locus]->cm == -1) {

			//Find suitable interpolation interval
			while (vec_pos[i_locus]->bp > pos_bp[closest_pos] && closest_pos < pos_bp.size()) closest_pos++;

			//Interpolate
			if (closest_pos < pos_bp.size()) {
				assert(vec_pos[i_locus]->bp < pos_bp[closest_pos]);
				assert(vec_pos[i_locus]->bp > pos_bp[closest_pos-1]);
				base = pos_cM[closest_pos-1];
				rate = (pos_cM[closest_pos] - pos_cM[closest_pos-1]) / (pos_bp[closest_pos] - pos_bp[closest_pos-1]);
				dist = (vec_pos[i_locus]->bp - pos_bp[closest_pos-1]);
				vec_pos[i_locus]->cm = base + rate * dist;
				n_interpolated ++;
				i_locus ++;
			} else break;
		} else i_locus ++;
	}

	//Set up last positions to be mean rate
	while (i_locus < vec_pos.size()) {
		base = pos_cM.back();
		dist = (vec_pos[i_locus]->bp - pos_bp.back());
		vec_pos[i_locus]->cm = base + mean_rate * dist;
		n_interpolated ++;
		i_locus ++;
	}
	return n_interpolated;
}

unsigned int variant_map::length() {
	return vec_pos.back()->bp - vec_pos[0]->bp + 1;
}

double variant_map::lengthcM() {
	return vec_pos.back()->cm - vec_pos[0]->cm;
}

supersite_summary variant_map::buildSupersites(const bitmatrix & hap_matrix, unsigned long n_hap) {
	supersite_summary summary;
	summary.total_variants = static_cast<uint32_t>(vec_pos.size());

	variant_to_site.assign(vec_pos.size(), 0u);
	variant_alt_code.assign(vec_pos.size(), 0u);
	variant_is_anchor.assign(vec_pos.size(), 0u);
	supersites.clear();
	supersites.reserve(vec_pos.size());
	supersite_alt_variant_index.clear();
	supersite_alt_variant_index.reserve(vec_pos.size());
	supersite_alt_variant_offset.clear();
	supersite_alt_variant_offset.reserve(vec_pos.size() + 1);
	supersite_codes.clear();

	if (vec_pos.empty()) {
		supersite_alt_variant_offset.push_back(0u);
		return summary;
	}

	if (hap_matrix.n_cols < vec_pos.size()) vrb.error("Supersite build: haplotype matrix has fewer columns than variants");
	if (hap_matrix.n_rows < n_hap) vrb.error("Supersite build: haplotype matrix has fewer rows than haplotypes");

	for (size_t i = 0; i < vec_pos.size();) {
		const std::string & chr = vec_pos[i]->chr;
		const int bp = vec_pos[i]->bp;
		size_t j = i + 1;
		while (j < vec_pos.size() && vec_pos[j]->chr == chr && vec_pos[j]->bp == bp) ++j;
		const size_t span = j - i;

		supersite_alt_variant_offset.push_back(static_cast<uint32_t>(supersite_alt_variant_index.size()));

		vector<uint32_t> block(span);
		for (size_t k = 0; k < span; ++k) block[k] = static_cast<uint32_t>(i + k);

		if (span > 1) {
			auto cmp = [&](uint32_t lhs, uint32_t rhs) {
				const std::string & alt_lhs = vec_pos[lhs]->alt;
				const std::string & alt_rhs = vec_pos[rhs]->alt;
				if (alt_lhs.size() != alt_rhs.size()) return alt_lhs.size() > alt_rhs.size();
				if (alt_lhs != alt_rhs) return alt_lhs < alt_rhs;
				return lhs < rhs;
			};
			sort(block.begin(), block.end(), cmp);
		}

		for (uint32_t idx : block) supersite_alt_variant_index.push_back(idx);
		if (span > 0) variant_is_anchor[block[0]] = 1u;

		supersite_desc desc;
		desc.first_variant_index = static_cast<uint32_t>(i);
		desc.variant_span = static_cast<uint16_t>(span);
		desc.n_alt = static_cast<uint16_t>(span);
		desc.is_super_site = span > 1;
		desc.bitwidth = desc.is_super_site ? required_bitwidth(desc.n_alt) : 1u;
		desc.panel_offset = desc.is_super_site ? static_cast<uint32_t>(supersite_codes.size()) : std::numeric_limits<uint32_t>::max();
		desc.codes_count = desc.is_super_site ? static_cast<uint32_t>(n_hap) : 0u;
		desc.conflicting_haps = 0u;

		const uint32_t site_index = static_cast<uint32_t>(supersites.size());
		for (size_t k = 0; k < span; ++k) {
			const uint32_t variant_idx = static_cast<uint32_t>(i + k);
			variant_to_site[variant_idx] = site_index;
		}

		if (desc.is_super_site) {
			size_t offset = supersite_codes.size();
			supersite_codes.resize(offset + static_cast<size_t>(n_hap), 0u);
			desc.panel_offset = static_cast<uint32_t>(offset);

			for (size_t local = 0; local < span; ++local) {
				variant_alt_code[block[local]] = static_cast<uint8_t>(local + 1);
			}

			for (unsigned long h = 0; h < n_hap; ++h) {
				uint8_t code = 0u;
				bool conflict = false;
				for (size_t local = 0; local < span; ++local) {
					uint32_t variant_idx = block[local];
					if (hap_matrix.get(static_cast<unsigned int>(h), variant_idx)) {
						uint8_t alt_code = static_cast<uint8_t>(local + 1);
						if (code == 0u) code = alt_code;
						else if (code != alt_code) conflict = true;
					}
				}
				supersite_codes[desc.panel_offset + h] = code;
				if (conflict) desc.conflicting_haps++;
			}

			summary.total_super_sites++;
			summary.collapsed_variants += static_cast<uint32_t>(span - 1);
			summary.haplotype_conflicts += desc.conflicting_haps;
		} else {
			variant_alt_code[i] = 1u;
		}

		supersites.push_back(desc);
		summary.total_sites++;
		i = j;
	}

	supersite_alt_variant_offset.push_back(static_cast<uint32_t>(supersite_alt_variant_index.size()));

	return summary;
}

void variant_map::buildPosteriorOffsets() {
	supersite_posterior_offset.clear();
	supersite_posterior_offset.reserve(supersites.size());
	uint32_t acc = 0u;
	for (size_t s = 0; s < supersites.size(); ++s) {
		supersite_posterior_offset.push_back(acc);
		uint32_t n_codes = static_cast<uint32_t>(supersites[s].n_alt) + 1u;
		acc += n_codes * HAP_NUMBER;
	}
	total_posterior_size = acc;
}

void variant_map::setGeneticMap(gmap_reader & readerGM) {
	tac.clock();
	int n_set = setCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	int n_interpolated = interpolateCentiMorgan(readerGM.pos_bp, readerGM.pos_cm);
	double baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
	vrb.bullet("cM interpolation [s=" + stb.str(n_set) + " / i=" + stb.str(n_interpolated) + "] (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet("Region length [" + stb.str(vec_pos.back()->bp-vec_pos[0]->bp+1) + " bp / " + stb.str(vec_pos.back()->cm-vec_pos[0]->cm, 1) + " cM]");
}

void variant_map::setGeneticMap() {
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm = vec_pos[l]->bp * 1.0 / 1e6;
	double baseline = vec_pos[0]->cm;
	for (int l = 0 ; l < vec_pos.size() ; l ++) vec_pos[l]->cm -= baseline;
	vrb.bullet("Region length [" + stb.str(vec_pos.back()->bp-vec_pos[0]->bp+1) + " bp / " + stb.str(vec_pos.back()->cm-vec_pos[0]->cm, 1) + " cM (assuming 1cM per Mb)]");
}

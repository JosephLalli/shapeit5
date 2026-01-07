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

#include <io/haplotype_reader.h>

using namespace std;

namespace {
std::string safe_string(const char *value) {
	return value ? std::string(value) : std::string(".");
}

std::string to_lower_ascii(const std::string &value) {
	std::string out = value;
	for (char &ch : out) {
		ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
	}
	return out;
}

std::string join_alleles(const std::vector<std::string> &alleles) {
	std::string out;
	for (size_t i = 0; i < alleles.size(); ++i) {
		if (i) out += ",";
		out += alleles[i];
	}
	return out;
}

std::string join_alt_alleles(const bcf1_t *line) {
	std::string alt_str;
	for (int ai = 1; ai < line->n_allele; ++ai) {
		if (ai > 1) alt_str += ",";
		alt_str += safe_string(line->d.allele[ai]);
	}
	return alt_str;
}

std::string join_alt_alleles(const std::vector<std::string> &alleles) {
	std::string out;
	for (size_t i = 1; i < alleles.size(); ++i) {
		if (i > 1) out += ",";
		out += alleles[i];
	}
	return out;
}

std::vector<std::string> extract_alleles(const bcf1_t *line) {
	std::vector<std::string> alleles;
	alleles.reserve(line->n_allele);
	for (int ai = 0; ai < line->n_allele; ++ai) {
		alleles.push_back(safe_string(line->d.allele[ai]));
	}
	return alleles;
}

struct allele_order {
	std::vector<std::string> alleles;
	std::vector<std::string> alleles_lower;
	std::vector<int> mapping;
	bool reordered;
};

allele_order lexicographic_order(const bcf1_t *line) {
	std::vector<std::string> original = extract_alleles(line);
	const int n = static_cast<int>(original.size());
	if (n == 0) {
		vrb.error("Empty allele list in VCF/BCF record");
	}
	std::vector<std::string> original_lower;
	original_lower.reserve(n);
	for (const auto &allele : original) {
		original_lower.push_back(to_lower_ascii(allele));
	}

	std::vector<int> alt_indices;
	alt_indices.reserve(n > 0 ? n - 1 : 0);
	for (int i = 1; i < n; ++i) alt_indices.push_back(i);
	std::stable_sort(alt_indices.begin(), alt_indices.end(),
	                 [&original_lower](int a, int b) { return original_lower[a] < original_lower[b]; });

	std::vector<std::string> sorted;
	std::vector<std::string> sorted_lower;
	sorted.reserve(n);
	sorted_lower.reserve(n);
	sorted.push_back(original[0]);
	sorted_lower.push_back(original_lower[0]);
	for (int idx : alt_indices) sorted.push_back(original[idx]);
	for (int idx : alt_indices) sorted_lower.push_back(original_lower[idx]);

	for (int i = 1; i < n; ++i) {
		if (sorted_lower[i] == sorted_lower[i - 1]) {
			vrb.error("Duplicate allele after case-insensitive normalization: " +
			          sorted[i - 1] + " and " + sorted[i]);
		}
	}

	std::vector<int> mapping(n, 0);
	if (n > 1) {
		for (int i = 1; i < n; ++i) {
			mapping[alt_indices[i - 1]] = i;
		}
	}

	bool reordered = false;
	for (int i = 1; i < n; ++i) {
		if (mapping[i] != i) {
			reordered = true;
			break;
		}
	}
	return {sorted, sorted_lower, mapping, reordered};
}
} // namespace

haplotype_reader::haplotype_reader(haplotype_set & _H, string _region, double _minPP,  int _nthreads, std::string _site_log_path) : H(_H) {
	nthreads = _nthreads;
	region = _region;
	minPP = _minPP;
	site_log_path = _site_log_path;
}

haplotype_reader::~haplotype_reader() {
	region = "";
}

void haplotype_reader::readHaplotypes(string ftruth, string festi, string ffreq, bool dupid) {
	tac.clock();
	vrb.title("Reading VCF/BCF input files");
	vrb.bullet("Validation ["  + ftruth + "]");
	vrb.bullet("Estimation ["  + festi + "]");
	vrb.bullet("Frequency  ["  + ffreq + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, ftruth.c_str()))) vrb.error("Problem opening index file for [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, festi.c_str()))) vrb.error("Problem opening index file for [" + festi + "]");
	if(!(bcf_sr_add_reader (sr, ffreq.c_str()))) vrb.error("Problem opening index file for [" + ffreq + "]");

	//IDs in truth
	int n_samples_truth = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_truth ; i ++) {
		string sample_id = string(sr->readers[0].header->samples[i]);
		if (dupid) sample_id = sample_id + "_" + sample_id;
		H.push(sample_id);
	}
	vrb.bullet("#Validation samples = " + stb.str(n_samples_truth));

	//IDs in estimation
	int n_samples_estimated = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mapping = vector < int > (n_samples_estimated, -1);
	for (int i = 0 ; i < n_samples_estimated ; i ++) {
		map < string, int > :: iterator itM = H.mapSamples.find(string(sr->readers[1].header->samples[i]));
		if (itM != H.mapSamples.end()) {
			mapping[i] = itM->second;
			H.IDXesti.push_back(itM->second);
		}
	}
	sort(H.IDXesti.begin(), H.IDXesti.end());
	vrb.bullet("#Estimation samples = " + stb.str(n_samples_estimated));
	vrb.bullet("#Overlapping samples = " + stb.str(H.IDXesti.size()));

	std::ofstream site_log;
	if (!site_log_path.empty()) {
		site_log.open(site_log_path.c_str());
		if (site_log.fail()) vrb.error("Impossible to create site log file [" + site_log_path + "]");
		site_log << "#site_index\ttruth_pos\ttruth_id\ttruth_n_allele\ttruth_ref\ttruth_alt\ttruth_missing\t"
		         << "est_pos\test_id\test_n_allele\test_ref\test_alt\test_missing\test_estimated\t"
		         << "freq_pos\tfreq_id\tfreq_n_allele\tfreq_ref\tfreq_alt\n";
	}
	const bool log_sites = site_log.is_open();
	int n_reordered_truth = 0;
	int n_reordered_est = 0;
	int n_reordered_freq = 0;

	//Read GT data
	int n_variant_tot = 0;
	float * vPP = NULL;
	float * vAF = NULL;
	int rAC=0, nAC=0, *vAC = NULL, rAN=0, nAN=0, *vAN = NULL, rAF = 0, nAF = 0, nPP = 0, rPP = 0;
	int nset = 0, *gt_arr_t = NULL, *gt_arr_e = NULL, ngt_arr_t = 0, ngt_arr_e = 0;
	bcf1_t * line_t, * line_e, * line_f;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 3) {
			line_t =  bcf_sr_get_line(sr, 0);
			line_e =  bcf_sr_get_line(sr, 1);
			line_f =  bcf_sr_get_line(sr, 2);
			if (line_f->n_allele >= 2) {
				//1. Unpack variant infos
				bcf_unpack(line_f, BCF_UN_ALL);
				bcf_unpack(line_t, BCF_UN_STR);
				bcf_unpack(line_e, BCF_UN_STR);
				if (line_t->n_allele != line_f->n_allele || line_e->n_allele != line_t->n_allele) {
					vrb.error("Allele count mismatch between validation/estimation/frequency records");
				}
				const allele_order order_truth = lexicographic_order(line_t);
				const allele_order order_est = lexicographic_order(line_e);
				const allele_order order_freq = lexicographic_order(line_f);
				if (order_truth.alleles_lower != order_est.alleles_lower) {
					vrb.error("Allele mismatch between validation and estimation after case-insensitive lexicographic sorting "
					          "(validation=" + join_alleles(order_truth.alleles) +
					          " estimation=" + join_alleles(order_est.alleles) + ")");
				}
				if (order_truth.alleles_lower != order_freq.alleles_lower) {
					vrb.error("Allele mismatch between validation and frequency after case-insensitive lexicographic sorting "
					          "(validation=" + join_alleles(order_truth.alleles) +
					          " frequency=" + join_alleles(order_freq.alleles) + ")");
				}
				if (order_truth.reordered) n_reordered_truth++;
				if (order_est.reordered) n_reordered_est++;
				if (order_freq.reordered) n_reordered_freq++;
				const int n_alts = static_cast<int>(order_truth.alleles.size()) - 1;
				int missing_truth = 0;
				int missing_est = 0;
				int estimated_count = 0;
				H.Positions.push_back(line_t->pos + 1);
				H.RSIDs.push_back(string(line_t->d.id));
				H.REFs.push_back(string(line_t->d.allele[0]));
				std::string alt_str = join_alt_alleles(order_truth.alleles);
				H.ALTs.push_back(alt_str);
				H.NAlts.push_back(static_cast<uint16_t>(n_alts));

				rAN = bcf_get_info_int32(sr->readers[2].header, line_f, "AN", &vAN, &nAN);
				rAC = bcf_get_info_int32(sr->readers[2].header, line_f, "AC", &vAC, &nAC);
				if (rAN != 1) vrb.error("AN field is needed");
				int alt_count = 0;
				if (rAC > 0 && nAC > 0) {
					for (int ai = 0; ai < rAC; ++ai) alt_count += vAC[ai];
				} else {
					rAF = bcf_get_info_float(sr->readers[2].header, line_f, "AF", &vAF, &nAF);
					if (rAF <= 0 || nAF <= 0) vrb.error("AC or AF field is needed");
					float af_sum = 0.0f;
					for (int ai = 0; ai < rAF; ++ai) af_sum += vAF[ai];
					alt_count = static_cast<int>(std::round(af_sum * static_cast<float>(vAN[0])));
				}
				int ref_count = vAN[0] - alt_count;
				if (ref_count < 0) ref_count = 0;
				H.MAC.push_back(std::min(alt_count, ref_count));
				H.MinorAlleles.push_back(alt_count < ref_count);

				//2. Validation
				bcf_get_genotypes(sr->readers[0].header, line_t, &gt_arr_t, &ngt_arr_t);
				for(int h = 0 ; h < 2 * n_samples_truth ; h += 2) {
					const int32_t raw0 = gt_arr_t[h+0];
					const int32_t raw1 = gt_arr_t[h+1];
					bool mi = (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing);
					int allele0 = bcf_gt_allele(raw0);
					int allele1 = bcf_gt_allele(raw1);
					if (allele0 < 0 || allele0 >= line_t->n_allele) { allele0 = 0; mi = true; }
					else allele0 = order_truth.mapping[allele0];
					if (allele1 < 0 || allele1 >= line_t->n_allele) { allele1 = 0; mi = true; }
					else allele1 = order_truth.mapping[allele1];
					if (allele0 < 0 || allele0 > n_alts) { allele0 = 0; mi = true; }
					if (allele1 < 0 || allele1 > n_alts) { allele1 = 0; mi = true; }
					H.Htrue[h+0].push_back(static_cast<uint16_t>(allele0));
					H.Htrue[h+1].push_back(static_cast<uint16_t>(allele1));
					H.Missing[h/2].push_back(mi);
					if (mi) missing_truth++;
				}

				//3. Estimation
				bcf_get_genotypes(sr->readers[1].header, line_e, &gt_arr_e, &ngt_arr_e);
				for(int h = 0 ; h < 2 * n_samples_estimated ; h += 2) {
					int index = mapping[h/2];
					if (index >= 0) {
						const int32_t raw0 = gt_arr_e[h+0];
						const int32_t raw1 = gt_arr_e[h+1];
						bool mi = (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing);
						int allele0 = bcf_gt_allele(raw0);
						int allele1 = bcf_gt_allele(raw1);
						if (allele0 < 0 || allele0 >= line_e->n_allele) { allele0 = 0; mi = true; }
						else allele0 = order_est.mapping[allele0];
						if (allele1 < 0 || allele1 >= line_e->n_allele) { allele1 = 0; mi = true; }
						else allele1 = order_est.mapping[allele1];
						if (allele0 < 0 || allele0 > n_alts) { allele0 = 0; mi = true; }
						if (allele1 < 0 || allele1 > n_alts) { allele1 = 0; mi = true; }
						H.Hesti[2*index+0].push_back(static_cast<uint16_t>(allele0));
						H.Hesti[2*index+1].push_back(static_cast<uint16_t>(allele1));
						H.MissingEst[index].push_back(mi);
						if (mi) missing_est++;
					}
				}

				//4. Probabilities
				rPP = bcf_get_format_float(sr->readers[1].header, line_e, "PP", &vPP, &nPP);
				const bool has_pp = (rPP == n_samples_estimated);
				for(int i = 0 ; i < n_samples_estimated ; i ++) {
					int index = mapping[i];
					if (index >= 0) {
						const bool pp_present = has_pp && !bcf_float_is_missing(vPP[i]);
						H.Hprob[index].push_back(pp_present);
						bool estimated_here = !H.MissingEst[index].back();
						if (pp_present) {
							string key = stb.str(H.Hprob[index].size() - 1) + "_" + stb.str(index);
							H.Vprob.insert(pair < string, float > ( key, vPP[i]));
							if (vPP[i] <= minPP) estimated_here = false;
						}
						H.Estimated[index].push_back(estimated_here);
						if (estimated_here) estimated_count++;
					}
				}

				if (log_sites) {
					const int site_index = H.n_variants + 1;
					const std::string alt_str_t = join_alt_alleles(order_truth.alleles);
					const std::string alt_str_e = join_alt_alleles(order_est.alleles);
					site_log << site_index << "\t" << (line_t->pos + 1) << "\t" << safe_string(line_t->d.id) << "\t"
					         << line_t->n_allele << "\t" << safe_string(line_t->d.allele[0]) << "\t" << alt_str_t
					         << "\t" << missing_truth << "\t" << (line_e->pos + 1) << "\t" << safe_string(line_e->d.id)
					         << "\t" << line_e->n_allele << "\t" << safe_string(line_e->d.allele[0]) << "\t" << alt_str_e
					         << "\t" << missing_est << "\t" << estimated_count << "\t" << (line_f->pos + 1) << "\t"
					         << safe_string(line_f->d.id) << "\t" << line_f->n_allele << "\t"
					         << safe_string(line_f->d.allele[0]) << "\t" << join_alt_alleles(order_freq.alleles) << "\n";
				}

				H.n_variants ++;

			}
		}
		n_variant_tot ++;
		if (n_variant_tot % 10000 == 0) vrb.bullet (stb.str(n_variant_tot) + " lines processed");
	}
	vrb.bullet("#Total variants = " + stb.str(n_variant_tot));
	vrb.bullet("#Overlapping variants = " + stb.str(H.n_variants));
	vrb.bullet("#Prob stored [PP field] = " + stb.str(H.Vprob.size()));
	if (n_reordered_truth > 0) {
		vrb.bullet("#Lexicographically reordered validation alleles at " + stb.str(n_reordered_truth) + " sites");
	}
	if (n_reordered_est > 0) {
		vrb.bullet("#Lexicographically reordered estimation alleles at " + stb.str(n_reordered_est) + " sites");
	}
	if (n_reordered_freq > 0) {
		vrb.bullet("#Lexicographically reordered frequency alleles at " + stb.str(n_reordered_freq) + " sites");
	}
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
	free(gt_arr_t); free(gt_arr_e);
	bcf_sr_destroy(sr);
}


void haplotype_reader::readHaplotypes(string ftruth, string festi, bool dupid) {
	tac.clock();
	vrb.title("Reading VCF/BCF input files");
	vrb.bullet("Validation ["  + ftruth + "]");
	vrb.bullet("Estimation ["  + festi + "]");

	bcf_srs_t * sr =  bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);
	if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1) vrb.error("Impossible to jump to region [" + region + "] in [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, ftruth.c_str()))) vrb.error("Problem opening index file for [" + ftruth + "]");
	if(!(bcf_sr_add_reader (sr, festi.c_str()))) vrb.error("Problem opening index file for [" + festi + "]");

	//IDs in truth
	int n_samples_truth = bcf_hdr_nsamples(sr->readers[0].header);
	for (int i = 0 ; i < n_samples_truth ; i ++) {
		string sample_id = string(sr->readers[0].header->samples[i]);
		if (dupid) sample_id = sample_id + "_" + sample_id;
		H.push(sample_id);
	}
	vrb.bullet("#Validation samples = " + stb.str(n_samples_truth));

	//IDs in estimation
	int n_samples_estimated = bcf_hdr_nsamples(sr->readers[1].header);
	vector < int > mapping = vector < int > (n_samples_estimated, -1);
	for (int i = 0 ; i < n_samples_estimated ; i ++) {
		map < string, int > :: iterator itM = H.mapSamples.find(string(sr->readers[1].header->samples[i]));
		if (itM != H.mapSamples.end()) {
			mapping[i] = itM->second;
			H.IDXesti.push_back(itM->second);
		}
	}
	sort(H.IDXesti.begin(), H.IDXesti.end());
	vrb.bullet("#Estimation samples = " + stb.str(n_samples_estimated));
	vrb.bullet("#Overlapping samples = " + stb.str(H.IDXesti.size()));

	std::ofstream site_log;
	if (!site_log_path.empty()) {
		site_log.open(site_log_path.c_str());
		if (site_log.fail()) vrb.error("Impossible to create site log file [" + site_log_path + "]");
		site_log << "#site_index\ttruth_pos\ttruth_id\ttruth_n_allele\ttruth_ref\ttruth_alt\ttruth_missing\t"
		         << "est_pos\test_id\test_n_allele\test_ref\test_alt\test_missing\test_estimated\n";
	}
	const bool log_sites = site_log.is_open();
	int n_reordered_truth = 0;
	int n_reordered_est = 0;

	//Read GT data
	int n_variant_tot = 0;
	float * vPP = NULL;
	int nPP = 0, rPP = 0;
	int nset = 0, *gt_arr_t = NULL, *gt_arr_e = NULL, ngt_arr_t = 0, ngt_arr_e = 0;
	bcf1_t * line_t, * line_e;
	while ((nset = bcf_sr_next_line (sr))) {
		if (nset == 2) {
			line_t =  bcf_sr_get_line(sr, 0);
			line_e =  bcf_sr_get_line(sr, 1);
			if (line_t->n_allele >= 2) {
				//1. Unpack variant infos
				bcf_unpack(line_t, BCF_UN_ALL);
				bcf_unpack(line_e, BCF_UN_STR);
				if (line_e->n_allele != line_t->n_allele) {
					vrb.error("Allele count mismatch between validation/estimation records");
				}
				const allele_order order_truth = lexicographic_order(line_t);
				const allele_order order_est = lexicographic_order(line_e);
				if (order_truth.alleles_lower != order_est.alleles_lower) {
					vrb.error("Allele mismatch between validation and estimation after case-insensitive lexicographic sorting "
					          "(validation=" + join_alleles(order_truth.alleles) +
					          " estimation=" + join_alleles(order_est.alleles) + ")");
				}
				if (order_truth.reordered) n_reordered_truth++;
				if (order_est.reordered) n_reordered_est++;
				const int n_alts = static_cast<int>(order_truth.alleles.size()) - 1;
				int missing_truth = 0;
				int missing_est = 0;
				int estimated_count = 0;
				H.Positions.push_back(line_t->pos + 1);
				H.RSIDs.push_back(string(line_t->d.id));
				H.REFs.push_back(string(line_t->d.allele[0]));
				std::string alt_str = join_alt_alleles(order_truth.alleles);
				H.ALTs.push_back(alt_str);
				H.NAlts.push_back(static_cast<uint16_t>(n_alts));

				//2. Validation
				int vAC = 0, vAN = 0;
				bcf_get_genotypes(sr->readers[0].header, line_t, &gt_arr_t, &ngt_arr_t);
				for(int h = 0 ; h < 2 * n_samples_truth ; h += 2) {
					const int32_t raw0 = gt_arr_t[h+0];
					const int32_t raw1 = gt_arr_t[h+1];
					bool mi = (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing);
					int allele0 = bcf_gt_allele(raw0);
					int allele1 = bcf_gt_allele(raw1);
					if (allele0 < 0 || allele0 >= line_t->n_allele) { allele0 = 0; mi = true; }
					else allele0 = order_truth.mapping[allele0];
					if (allele1 < 0 || allele1 >= line_t->n_allele) { allele1 = 0; mi = true; }
					else allele1 = order_truth.mapping[allele1];
					if (allele0 < 0 || allele0 > n_alts) { allele0 = 0; mi = true; }
					if (allele1 < 0 || allele1 > n_alts) { allele1 = 0; mi = true; }
					H.Htrue[h+0].push_back(static_cast<uint16_t>(allele0));
					H.Htrue[h+1].push_back(static_cast<uint16_t>(allele1));
					H.Missing[h/2].push_back(mi);
					if (!mi) {
						vAC += (allele0 != 0) + (allele1 != 0);
						vAN += 2;
					} else {
						missing_truth++;
					}
				}
				H.MAC.push_back(min(vAC, (vAN - vAC)));
				H.MinorAlleles.push_back(vAC < (vAN - vAC));

				//3. Estimation
				bcf_get_genotypes(sr->readers[1].header, line_e, &gt_arr_e, &ngt_arr_e);
				for(int h = 0 ; h < 2 * n_samples_estimated ; h += 2) {
					int index = mapping[h/2];
					if (index >= 0) {
						const int32_t raw0 = gt_arr_e[h+0];
						const int32_t raw1 = gt_arr_e[h+1];
						bool mi = (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing);
						int allele0 = bcf_gt_allele(raw0);
						int allele1 = bcf_gt_allele(raw1);
						if (allele0 < 0 || allele0 >= line_e->n_allele) { allele0 = 0; mi = true; }
						else allele0 = order_est.mapping[allele0];
						if (allele1 < 0 || allele1 >= line_e->n_allele) { allele1 = 0; mi = true; }
						else allele1 = order_est.mapping[allele1];
						if (allele0 < 0 || allele0 > n_alts) { allele0 = 0; mi = true; }
						if (allele1 < 0 || allele1 > n_alts) { allele1 = 0; mi = true; }
						H.Hesti[2*index+0].push_back(static_cast<uint16_t>(allele0));
						H.Hesti[2*index+1].push_back(static_cast<uint16_t>(allele1));
						H.MissingEst[index].push_back(mi);
						if (mi) missing_est++;
					}
				}

				//4. Probabilities
				rPP = bcf_get_format_float(sr->readers[1].header, line_e, "PP", &vPP, &nPP);
				const bool has_pp = (rPP == n_samples_estimated);
				for(int i = 0 ; i < n_samples_estimated ; i ++) {
					int index = mapping[i];
					if (index >= 0) {
						const bool pp_present = has_pp && !bcf_float_is_missing(vPP[i]);
						H.Hprob[index].push_back(pp_present);
						bool estimated_here = !H.MissingEst[index].back();
						if (pp_present) {
							string key = stb.str(H.Hprob[index].size() - 1) + "_" + stb.str(index);
							H.Vprob.insert(pair < string, float > ( key, vPP[i]));
							if (vPP[i] <= minPP) estimated_here = false;
						}
						H.Estimated[index].push_back(estimated_here);
						if (estimated_here) estimated_count++;
					}
				}

				if (log_sites) {
					const int site_index = H.n_variants + 1;
					const std::string alt_str_e = join_alt_alleles(order_est.alleles);
					site_log << site_index << "\t" << (line_t->pos + 1) << "\t" << safe_string(line_t->d.id) << "\t"
					         << line_t->n_allele << "\t" << safe_string(line_t->d.allele[0]) << "\t" << alt_str
					         << "\t" << missing_truth << "\t" << (line_e->pos + 1) << "\t" << safe_string(line_e->d.id)
					         << "\t" << line_e->n_allele << "\t" << safe_string(line_e->d.allele[0]) << "\t" << alt_str_e
					         << "\t" << missing_est << "\t" << estimated_count << "\n";
				}

				H.n_variants ++;
			}
		}
		n_variant_tot ++;
		if (n_variant_tot % 10000 == 0) vrb.bullet (stb.str(n_variant_tot) + " lines processed");
	}
	vrb.bullet("#Total variants = " + stb.str(n_variant_tot));
	vrb.bullet("#Overlapping variants = " + stb.str(H.n_variants));
	vrb.bullet("#Prob stored [PP field] = " + stb.str(H.Vprob.size()));
	if (n_reordered_truth > 0) {
		vrb.bullet("#Lexicographically reordered validation alleles at " + stb.str(n_reordered_truth) + " sites");
	}
	if (n_reordered_est > 0) {
		vrb.bullet("#Lexicographically reordered estimation alleles at " + stb.str(n_reordered_est) + " sites");
	}
	vrb.bullet("Timing: " + stb.str(tac.rel_time()*1.0/1000, 2) + "s");
	free(gt_arr_t); free(gt_arr_e);
	bcf_sr_destroy(sr);
}

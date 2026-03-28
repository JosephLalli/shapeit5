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

#include <io/genotype_reader/genotype_reader_header.h>

namespace {
inline void rare_mask_set(std::array<uint64_t, 4>& mask, uint16_t code) {
	const uint16_t idx = code >> 6;
	const uint16_t bit = code & 63u;
	mask[idx] |= (1ULL << bit);
}
}

void genotype_reader::scanGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF scanning");

	//File idx
	int32_t idx_file_main = 0;
	int32_t idx_file_ref = panels[1] ? 1 : -1;
	int32_t idx_file_scaf = panels[2] ? (panels[1] + panels[2]) : -1;

	//Initialize synced reader
	bcf_srs_t *sr = bcf_sr_init();
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if (nthreads > 1) bcf_sr_set_threads(sr, nthreads);

	//Set region filter
	if (!region.empty()) {
		if (bcf_sr_set_regions(sr, region.c_str(), 0) == -1)
			vrb.error("Impossible to jump to region [" + region + "]");
		if (bcf_sr_set_targets(sr, region.c_str(), 0, 0) == -1)
			vrb.error("Impossible to constrain to region [" + region + "]");
	}

	//Add file readers
	for (int f = 0; f < 3; f++) {
		if (panels[f]) {
			if (!bcf_sr_add_reader(sr, filenames[f].c_str())) {
				switch (sr->errnum) {
				case not_bgzf:
					vrb.error("Opening [" + filenames[f] + "]: not compressed with bgzip");
					break;
				case idx_load_failed:
					vrb.error("Opening [" + filenames[f] + "]: impossible to load index file");
					break;
				case file_type_error:
					vrb.error("Opening [" + filenames[f] + "]: file format not supported by HTSlib");
					break;
				default:
					vrb.error("Opening [" + filenames[f] + "]: unknown error");
					break;
				}
			}
		}
	}

	//Sample processing
	n_main_samples = bcf_hdr_nsamples(sr->readers[idx_file_main].header);
	n_ref_samples = panels[1] ? bcf_hdr_nsamples(sr->readers[idx_file_ref].header) : 0;
	if ((n_main_samples + n_ref_samples) < 50) {
		if (n_ref_samples) vrb.error("Less than 50 samples is not enough to get reliable phasing");
		else vrb.error("Less than 50 samples is not enough to get reliable phasing, consider using a reference panel to increase sample size");
	}

	uint32_t n_variants_with_scaffold = 0;
	uint32_t n_variants_noverlap = 0, n_variants_multiallelic = 0, n_variants_notsnp = 0, n_variants_rare = 0, n_variants_nscaf = 0, n_variants_nref = 0;
	uint32_t n_rare_alt_masked = 0;
	uint32_t n_variants_main_format = 0, n_variants_ref_format = 0, n_variants_scaf_format = 0;
	has_multiallelic_records = false;
	has_binary_haplotype = false;
	rare_alt_masks.clear();

	//Buffers for AC/AN extraction (reused across iterations)
	int32_t *vAC = nullptr, *vAN = nullptr;
	int nAC = 0, nAN = 0;

	while (bcf_sr_next_line(sr)) {

		//By defaults, we do not want the variants
		variant_mask.push_back(false);

		//See which file has a record
		bool has_main = bcf_sr_has_line(sr, idx_file_main);
		bool has_ref = (panels[1] && bcf_sr_has_line(sr, idx_file_ref));
		bool has_scaf = (panels[2] && bcf_sr_has_line(sr, idx_file_scaf));

		//Extract variant metadata from first available record
		bcf1_t *rec = nullptr;
		int rec_file_idx = -1;
		if (has_main) {
			rec = bcf_sr_get_line(sr, idx_file_main);
			rec_file_idx = idx_file_main;
		} else if (has_ref) {
			rec = bcf_sr_get_line(sr, idx_file_ref);
			rec_file_idx = idx_file_ref;
		} else if (has_scaf) {
			rec = bcf_sr_get_line(sr, idx_file_scaf);
			rec_file_idx = idx_file_scaf;
		}

		if (rec == nullptr || rec->n_allele < 2) continue;

		//Extract metadata
		bcf_hdr_t *hdr = sr->readers[rec_file_idx].header;
		std::string chr = std::string(bcf_hdr_id2name(hdr, rec->rid));
		uint32_t pos = rec->pos + 1;
		std::string rsid = std::string(rec->d.id);
		std::string ref_allele = std::string(rec->d.allele[0]);
		int n_allele = rec->n_allele;

		//Build alt string and alts vector
		std::vector<std::string> alts;
		for (int ai = 1; ai < n_allele; ++ai) {
			alts.push_back(std::string(rec->d.allele[ai]));
		}
		std::string alt_str;
		if (alts.size() == 1) {
			alt_str = alts[0];
		} else {
			for (size_t ai = 0; ai < alts.size(); ++ai) {
				if (ai) alt_str += ",";
				alt_str += alts[ai];
			}
		}

		const bool record_multiallelic = (n_allele > 2);

		//Not in reference panel
		if (panels[1] && has_main && !has_ref) { n_variants_noverlap++; continue; }

		//In reference, but not in main panel
		if (panels[1] && !has_main && has_ref) { n_variants_nref++; continue; }

		//In scaffold, but not in main panel
		if (panels[2] && !has_main && has_scaf) { n_variants_nscaf++; continue; }

		//Main panel format check (must have valid record)
		if (!has_main) { n_variants_main_format++; continue; }

		//Scaffold panel - count sites with scaffold data
		if (panels[2] && has_scaf) {
			n_variants_with_scaffold++;
		}

		//Keep SNPs only
		if (filter_snp_only) {
			auto is_snp_base = [](const std::string& base) {
				return (base == "A") || (base == "T") || (base == "G") || (base == "C");
			};
			bool bref = is_snp_base(ref_allele);
			bool balt = !alts.empty();
			for (const auto& a : alts) {
				balt &= is_snp_base(a);
			}
			n_variants_notsnp += (!bref || !balt);
			if (!bref || !balt) continue;
		}

		//Keep common only
		std::array<uint64_t, 4> rare_mask = {0u, 0u, 0u, 0u};
		if (filter_min_maf > 0) {
			int rAC = bcf_get_info_int32(hdr, rec, "AC", &vAC, &nAC);
			int rAN = bcf_get_info_int32(hdr, rec, "AN", &vAN, &nAN);

			if (record_multiallelic) {
				const int n_alts = n_allele - 1;
				if (rAC >= 1 && rAN == 1 && vAN[0] > 0) {
					for (int ai = 0; ai < n_alts; ++ai) {
						const int32_t ac = (ai < rAC) ? vAC[ai] : 0;
						const float af = static_cast<float>(ac) / static_cast<float>(vAN[0]);
						const float maf = std::min(1.0f - af, af);
						if (maf < filter_min_maf) {
							rare_mask_set(rare_mask, static_cast<uint16_t>(ai + 1));
							n_rare_alt_masked++;
						}
					}
				} else {
					for (int ai = 0; ai < n_alts; ++ai) {
						rare_mask_set(rare_mask, static_cast<uint16_t>(ai + 1));
						n_rare_alt_masked++;
					}
				}
			} else {
				float af = 0.0f;
				if (rAC >= 1 && rAN == 1 && vAN[0] > 0) {
					uint32_t ac_sum = 0;
					for (int ai = 0; ai < rAC; ++ai) {
						if (vAC[ai] > 0) ac_sum += static_cast<uint32_t>(vAC[ai]);
					}
					af = static_cast<float>(ac_sum) / static_cast<float>(vAN[0]);
				}

				float maf = std::min(1.0f - af, af);
				n_variants_rare += (maf < filter_min_maf);
				if (maf < filter_min_maf) continue;
			}
		}

		//Push variant information
		const uint16_t n_alts_count = (n_allele > 0) ? static_cast<uint16_t>(n_allele - 1) : 0u;
		V.push(new variant(chr, pos, rsid, ref_allele, alt_str, n_alts_count, V.size()));
		rare_alt_masks.push_back(rare_mask);

		//Flag it!
		variant_mask.back() = true;
		n_variants++;
		if (record_multiallelic) {
			++n_variants_multiallelic;
			has_multiallelic_records = true;
		}
	}

	//Cleanup
	free(vAC);
	free(vAN);
	bcf_sr_destroy(sr);

	n_supersites = n_variants_multiallelic;
	if (has_multiallelic_records && has_binary_haplotype) {
		vrb.error("Multiallelic records detected but reference/scaffold panel uses binary haplotype format; use BCF/VCF for multiallelic inputs");
	}

	vrb.bullet("VCF/BCF scanning done (" + stb.str(tac.rel_time()*1.0/1000, 2) + "s)");
	vrb.bullet2("Variants [#sites=" + stb.str(n_variants) + " / region=" + region + "]");
	if (n_variants_with_scaffold) vrb.bullet3(stb.str(n_variants_with_scaffold) + " sites with scaffold data");
	if (n_variants_noverlap) vrb.bullet3(stb.str(n_variants_noverlap) + " sites removed in main panel [not in reference panel]");
	if (n_variants_multiallelic) vrb.bullet3(stb.str(n_variants_multiallelic) + " multiallelic sites retained");
	if (n_variants_notsnp) vrb.bullet3(stb.str(n_variants_notsnp) + " sites removed in main panel [not SNPs]");
	if (n_variants_rare) vrb.bullet3(stb.str(n_variants_rare) + " sites removed in main panel [below MAF threshold]");
	if (n_rare_alt_masked) vrb.bullet3(stb.str(n_rare_alt_masked) + " ALT alleles masked in multiallelic sites [below MAF threshold]");
	if (n_variants_nref) vrb.bullet3(stb.str(n_variants_nref) + " sites removed in reference panel [not in main panel]");
	if (n_variants_nscaf) vrb.bullet3(stb.str(n_variants_nscaf) + " sites removed in scaffold panel [not in main panel]");
	if (n_variants_main_format) vrb.bullet3(stb.str(n_variants_main_format) + " sites removed [record in main panel not in a supported format]");
	if (n_variants_ref_format) vrb.bullet3(stb.str(n_variants_ref_format) + " sites removed [record in reference panel not in a supported format]");
	if (n_variants_scaf_format) vrb.bullet3(stb.str(n_variants_scaf_format) + " sites removed [record in scaffold panel not in a supported format]");
	if (n_variants == 0) vrb.error("No variants to be phased!");
	vrb.bullet2("Samples [#target=" + stb.str(n_main_samples) + " / #reference=" + stb.str(n_ref_samples) + "]");
}

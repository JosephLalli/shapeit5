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
inline bool rare_mask_contains(const std::array<uint64_t, 4>& mask, uint8_t code) {
	if (code == 0u) return false;
	const uint16_t idx = static_cast<uint16_t>(code) >> 6;
	const uint16_t bit = static_cast<uint16_t>(code) & 63u;
	return (mask[idx] & (1ULL << bit)) != 0u;
}
}

void genotype_reader::readGenotypes() {
	tac.clock();
	vrb.wait("  * VCF/BCF parsing");
	if (has_multiallelic_records && has_binary_haplotype) {
		vrb.error("Multiallelic records detected but reference/scaffold panel uses binary haplotype format; use BCF/VCF for multiallelic inputs");
	}

	//File idx
	int32_t idx_file_main = 0;
	int32_t idx_file_ref = panels[1]?1:-1;
	int32_t idx_file_scaf = panels[2]?(panels[1]+panels[2]):-1;

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

	//Opening file(s)
	for (int32_t f = 0 ; f < 3 ; f ++) {
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

	//Main sample IDs processing
	std::vector < std::string > main_names;
	std::map < std::string, uint32_t > map_names;
	bcf_hdr_t *hdr_main = sr->readers[idx_file_main].header;
	int32_t header_main_samples = bcf_hdr_nsamples(hdr_main);
	if (header_main_samples != static_cast<int32_t>(n_main_samples)) {
		vrb.error("Main panel sample count mismatch between scan and parse");
	}
	main_names.reserve(static_cast<size_t>(header_main_samples));
	for (int32_t i = 0 ; i < header_main_samples ; i ++) {
		main_names.push_back(std::string(hdr_main->samples[i]));
	}
	for (int32_t i = 0 ; i < static_cast<int32_t>(n_main_samples) ; i ++) {
		G.vecG[i]->name = main_names[i];
		map_names.insert(std::pair < std::string, int32_t > (G.vecG[i]->name, i));
	}

	//Scaffold sample IDs processing
	int32_t n_scaf_samples = 0, n_with_scaffold = 0;
	std::vector < uint32_t > mappingS2M;
	if (panels[2]) {
		std::vector < std::string > scaf_names;
		bcf_hdr_t *hdr_scaf = sr->readers[idx_file_scaf].header;
		n_scaf_samples = bcf_hdr_nsamples(hdr_scaf);
		scaf_names.reserve(static_cast<size_t>(n_scaf_samples));
		for (int32_t i = 0 ; i < n_scaf_samples ; i ++) {
			scaf_names.push_back(std::string(hdr_scaf->samples[i]));
		}
		mappingS2M = std::vector < uint32_t > (n_scaf_samples, -1);
		for (int32_t i = 0 ; i < n_scaf_samples ; i ++) {
			std::map < std::string, uint32_t > :: iterator it = map_names.find(scaf_names[i]);
			if (it != map_names.end()) { mappingS2M[i] = it->second; n_with_scaffold++; }
		}
		vrb.bullet3(stb.str(n_with_scaffold) + " samples with scaffold data");
	}

	//Genotype buffers (allocated by bcf_get_genotypes)
	int32_t * main_buffer = nullptr, * ref_buffer = nullptr, * scaf_buffer = nullptr;
	int n_main_buffer = 0, n_ref_buffer = 0, n_scaf_buffer = 0;

	//Parsing VCF/BCF
	uint32_t i_variant_total = 0, i_variant_kept = 0, i_supersite_kept = 0;
	while (bcf_sr_next_line(sr)) {
		if (i_variant_total >= variant_mask.size()) {
			vrb.error("Variant mask shorter than number of scanned records");
		}
		if (variant_mask[i_variant_total]) {
			const bool has_main = bcf_sr_has_line(sr, idx_file_main);
			const bool has_ref = panels[1] && bcf_sr_has_line(sr, idx_file_ref);
			const bool has_scaf = panels[2] && bcf_sr_has_line(sr, idx_file_scaf);
			if (!has_main) {
				vrb.error("Main panel record missing during VCF/BCF parsing");
			}
			bcf1_t *rec_main = bcf_sr_get_line(sr, idx_file_main);
			if (!rec_main) {
				vrb.error("Main panel record missing during VCF/BCF parsing");
			}
			const bool record_multiallelic = (rec_main->n_allele > 2);
			const int n_alts = rec_main->n_allele > 0 ? (rec_main->n_allele - 1) : 0;
			const int ss_idx = record_multiallelic ? static_cast<int>(i_supersite_kept) : -1;
			const size_t ss_offset = record_multiallelic ? static_cast<size_t>(ss_idx) * H.n_hap : 0u;
			const std::array<uint64_t, 4> empty_rare_mask = {0u, 0u, 0u, 0u};
			const std::array<uint64_t, 4>& rare_mask = (record_multiallelic && i_variant_kept < rare_alt_masks.size())
				? rare_alt_masks[i_variant_kept]
				: empty_rare_mask;

			// ====== Retrieve MAIN data ============== //
			int nGT_main = bcf_get_genotypes(hdr_main, rec_main, &main_buffer, &n_main_buffer);
			if (nGT_main <= 0) {
				vrb.error("Missing genotypes in main panel [" + filenames[0] + "]");
			}
			if (n_main_samples == 0 || (nGT_main % static_cast<int>(n_main_samples)) != 0) {
				vrb.error("Unexpected ploidy in main panel [" + filenames[0] + "]");
			}
			int ploidy_main = nGT_main / static_cast<int>(n_main_samples);
			if (ploidy_main != 2) {
				vrb.error("Only diploid genotypes are supported in main panel [" + filenames[0] + "]");
			}
			for(int32_t i = 0 ; i < 2 * static_cast<int32_t>(n_main_samples) ; i += 2) {
				genotype* g = G.vecG[DIV2(i)];
				if (record_multiallelic) {
					const int32_t raw0 = main_buffer[i+0];
					const int32_t raw1 = main_buffer[i+1];
					const bool mi = (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing);
					unsigned char &vbyte = g->Variants[DIV2(i_variant_kept)];
					if (mi) {
						g->setSupersiteObservedGt(ss_idx, 0u, 0u, 0x3u);
						VAR_SET_MIS(MOD2(i_variant_kept), vbyte);
						VAR_CLR_HAP0(MOD2(i_variant_kept), vbyte);
						VAR_CLR_HAP1(MOD2(i_variant_kept), vbyte);
						V.vec_pos[i_variant_kept]->cmis++;
						n_genotypes[3] ++;
						continue;
					}
					int allele0 = bcf_gt_allele(raw0);
					int allele1 = bcf_gt_allele(raw1);
					if (allele0 < 0 || allele0 > n_alts) allele0 = 0;
					if (allele1 < 0 || allele1 > n_alts) allele1 = 0;
					const uint8_t code0 = static_cast<uint8_t>(allele0);
					const uint8_t code1 = static_cast<uint8_t>(allele1);
					if (rare_mask_contains(rare_mask, code0) || rare_mask_contains(rare_mask, code1)) {
						g->setSupersiteObservedGt(ss_idx, 0u, 0u, 0x3u);
						VAR_SET_MIS(MOD2(i_variant_kept), vbyte);
						VAR_CLR_HAP0(MOD2(i_variant_kept), vbyte);
						VAR_CLR_HAP1(MOD2(i_variant_kept), vbyte);
						V.vec_pos[i_variant_kept]->cmis++;
						n_genotypes[3] ++;
						continue;
					}
					g->setSupersiteObservedGt(ss_idx, code0, code1, 0u);
					const unsigned sample_idx = static_cast<unsigned>(DIV2(i));
					const size_t hap0 = ss_offset + sample_idx * 2u;
					const size_t hap1 = hap0 + 1u;
					if (hap1 < H.H_supersite_codes.size()) {
						H.H_supersite_codes[hap0] = code0;
						H.H_supersite_codes[hap1] = code1;
					}
					const bool a0 = (code0 != 0);
					const bool a1 = (code1 != 0);
					a0 ? VAR_SET_HAP0(MOD2(i_variant_kept), vbyte) : VAR_CLR_HAP0(MOD2(i_variant_kept), vbyte);
					a1 ? VAR_SET_HAP1(MOD2(i_variant_kept), vbyte) : VAR_CLR_HAP1(MOD2(i_variant_kept), vbyte);
					if (code0 == code1) VAR_SET_HOM(MOD2(i_variant_kept), vbyte);
					else VAR_SET_HET(MOD2(i_variant_kept), vbyte);
					V.vec_pos[i_variant_kept]->cref += (code0 == 0) + (code1 == 0);
					V.vec_pos[i_variant_kept]->calt += (code0 > 0) + (code1 > 0);
					n_genotypes[(code0 > 0) + (code1 > 0)] ++;
				} else {
					bool a0 = (bcf_gt_allele(main_buffer[i+0])==1);
					bool a1 = (bcf_gt_allele(main_buffer[i+1])==1);
					bool mi = (main_buffer[i+0] == bcf_gt_missing || main_buffer[i+1] == bcf_gt_missing);
					bool he = !mi && a0 != a1;
					bool ho = !mi && a0 == a1;
					if (a0) VAR_SET_HAP0(MOD2(i_variant_kept), g->Variants[DIV2(i_variant_kept)]);
					if (a1) VAR_SET_HAP1(MOD2(i_variant_kept), g->Variants[DIV2(i_variant_kept)]);
					if (mi) VAR_SET_MIS(MOD2(i_variant_kept), g->Variants[DIV2(i_variant_kept)]);
					if (he) VAR_SET_HET(MOD2(i_variant_kept), g->Variants[DIV2(i_variant_kept)]);
					if (mi) { V.vec_pos[i_variant_kept]->cmis++; n_genotypes[3] ++; }
					else { V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1); V.vec_pos[i_variant_kept]->calt += a0+a1; n_genotypes[a0+a1] ++; }
				}
			}


			// ====== Retrieve REFERENCE data ============== //
			if (panels[1]) {
				if (!has_ref) {
					vrb.error("Reference panel record missing during VCF/BCF parsing");
				}
				bcf_hdr_t *hdr_ref = sr->readers[idx_file_ref].header;
				bcf1_t *rec_ref = bcf_sr_get_line(sr, idx_file_ref);
				int nGT_ref = bcf_get_genotypes(hdr_ref, rec_ref, &ref_buffer, &n_ref_buffer);
				if (nGT_ref <= 0) {
					vrb.error("Missing genotypes in reference panel [" + filenames[idx_file_ref] + "]");
				}
				if (n_ref_samples == 0 || (nGT_ref % static_cast<int>(n_ref_samples)) != 0) {
					vrb.error("Unexpected ploidy in reference panel [" + filenames[idx_file_ref] + "]");
				}
				int ploidy_ref = nGT_ref / static_cast<int>(n_ref_samples);
				if (ploidy_ref != 2) {
					vrb.error("Only diploid genotypes are supported in reference panel [" + filenames[idx_file_ref] + "]");
				}
				for(int32_t i = 0 ; i < 2 * static_cast<int32_t>(n_ref_samples) ; i += 2) {
					if (record_multiallelic) {
						const int32_t raw0 = ref_buffer[i+0];
						const int32_t raw1 = ref_buffer[i+1];
						if (raw0 == bcf_gt_missing || raw1 == bcf_gt_missing) vrb.error("Missing genotype(s) in reference panel");
						if (!bcf_gt_is_phased(ref_buffer[i+1])) vrb.error("Unphased genotype(s) in reference panel");
						int allele0 = bcf_gt_allele(raw0);
						int allele1 = bcf_gt_allele(raw1);
						if (allele0 < 0 || allele0 > n_alts) allele0 = 0;
						if (allele1 < 0 || allele1 > n_alts) allele1 = 0;
						const uint8_t code0 = static_cast<uint8_t>(allele0);
						const uint8_t code1 = static_cast<uint8_t>(allele1);
						const unsigned ref_sample = static_cast<unsigned>(DIV2(i));
						const size_t hap0 = ss_offset + (2u * n_main_samples) + ref_sample * 2u;
						const size_t hap1 = hap0 + 1u;
						if (hap1 < H.H_supersite_codes.size()) {
							H.H_supersite_codes[hap0] = code0;
							H.H_supersite_codes[hap1] = code1;
						}
						const bool a0 = (code0 != 0);
						const bool a1 = (code1 != 0);
						H.H_opt_hap.set(i+2*n_main_samples+0, i_variant_kept, a0);
						H.H_opt_hap.set(i+2*n_main_samples+1, i_variant_kept, a1);
						V.vec_pos[i_variant_kept]->cref += (code0 == 0) + (code1 == 0);
						V.vec_pos[i_variant_kept]->calt += (code0 > 0) + (code1 > 0);
						n_alleles[a0]++; n_alleles[a1]++;
					} else {
						bool a0 = (bcf_gt_allele(ref_buffer[i+0])==1);
						bool a1 = (bcf_gt_allele(ref_buffer[i+1])==1);
						if (ref_buffer[i+0] == bcf_gt_missing || ref_buffer[i+1] == bcf_gt_missing) vrb.error("Missing genotype(s) in reference panel");
						if (!bcf_gt_is_phased(ref_buffer[i+1])) vrb.error("Unphased genotype(s) in reference panel");
						H.H_opt_hap.set(i+2*n_main_samples+0, i_variant_kept, a0);
						H.H_opt_hap.set(i+2*n_main_samples+1, i_variant_kept, a1);
						V.vec_pos[i_variant_kept]->cref += (1-a0)+(1-a1);
						V.vec_pos[i_variant_kept]->calt += a0+a1;
						n_alleles[a0]++; n_alleles[a1]++;
					}
				}
			}

			// ====== Retrieve SCAFFOLD data ============== //
			if (panels[2] && has_scaf && !record_multiallelic) {
				bcf_hdr_t *hdr_scaf = sr->readers[idx_file_scaf].header;
				bcf1_t *rec_scaf = bcf_sr_get_line(sr, idx_file_scaf);
				int nGT_scaf = bcf_get_genotypes(hdr_scaf, rec_scaf, &scaf_buffer, &n_scaf_buffer);
				if (nGT_scaf <= 0) {
					vrb.error("Missing genotypes in scaffold panel [" + filenames[idx_file_scaf] + "]");
				}
				if (n_scaf_samples == 0 || (nGT_scaf % static_cast<int>(n_scaf_samples)) != 0) {
					vrb.error("Unexpected ploidy in scaffold panel [" + filenames[idx_file_scaf] + "]");
				}
				int ploidy_scaf = nGT_scaf / static_cast<int>(n_scaf_samples);
				if (ploidy_scaf != 2) {
					vrb.error("Only diploid genotypes are supported in scaffold panel [" + filenames[idx_file_scaf] + "]");
				}
				for(int32_t i = 0 ; i < 2 * n_scaf_samples ; i += 2) {
					int32_t ind = mappingS2M[DIV2(i)];
					if (ind >= 0) {
						bool sa0 = (bcf_gt_allele(scaf_buffer[i+0])==1);
						bool sa1 = (bcf_gt_allele(scaf_buffer[i+1])==1);
						bool sph = (bcf_gt_is_phased(scaf_buffer[i+0]) || bcf_gt_is_phased(scaf_buffer[i+1]));
						bool smi = (scaf_buffer[i+0] == bcf_gt_missing || scaf_buffer[i+1] == bcf_gt_missing);
						if ((sa0 != sa1) && !smi && sph && VAR_GET_HET(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)])) {
							VAR_SET_SCA(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							sa0?VAR_SET_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP0(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							sa1?VAR_SET_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]):VAR_CLR_HAP1(MOD2(i_variant_kept), G.vecG[ind]->Variants[DIV2(i_variant_kept)]);
							n_genotypes[4] ++;
						}
					}
				}
			}
			vrb.progress("  * VCF/BCF parsing", i_variant_kept * 1.0 / n_variants);
			if (record_multiallelic) i_supersite_kept++;
			i_variant_kept ++;
		}
		i_variant_total++;
	}
	free(main_buffer);
	free(ref_buffer);
	free(scaf_buffer);
	bcf_sr_destroy(sr);

	// Report
	uint64_t n_genotypes_total = accumulate(n_genotypes.begin(), n_genotypes.begin() + 4, 0UL);
	std::string str0 = "0/0=" + stb.str(n_genotypes[0]*100.0/n_genotypes_total, 3) + "%";
	std::string str1 = "0/1=" + stb.str(n_genotypes[1]*100.0/n_genotypes_total, 3) + "%";
	std::string str2 = "1/1=" + stb.str(n_genotypes[2]*100.0/n_genotypes_total, 3) + "%";
	std::string str3 = "./.=" + stb.str(n_genotypes[3]*100.0/n_genotypes_total, 3) + "%";
	std::string str4 = "0|1=" + stb.str(n_genotypes[4]*100.0/n_genotypes[1], 3) + "%";
	std::string str5 = stb.str(tac.rel_time()*1.0/1000, 2) + "s";
	vrb.bullet("VCF/BCF parsing done ("+str5+")");
	vrb.bullet2("Genotypes ["+str0+", "+str1+", "+str2+", "+str3+", "+str4+"]");

	if (panels[1]) {
		uint64_t n_alleles_total = n_alleles[0] + n_alleles[1];
		str0 = "0=" + stb.str(n_alleles[0]*100.0/n_alleles_total, 3) + "%";
		str1 = "1=" + stb.str(n_alleles[1]*100.0/n_alleles_total, 3) + "%";
		vrb.bullet2("Reference haplotypes ["+str0+", "+str1+"]");
	}
}

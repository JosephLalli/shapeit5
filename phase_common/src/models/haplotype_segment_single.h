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

#ifndef _HAPLOTYPE_SEGMENT_SINGLE_H
#define _HAPLOTYPE_SEGMENT_SINGLE_H

#include <utils/otools.h>
#include <objects/compute_job.h>
#include <objects/hmm_parameters.h>
#include <models/super_site_macros.h>

#include <immintrin.h>
#include <boost/align/aligned_allocator.hpp>

class haplotype_segment_single {
private:
	//EXTERNAL DATA
	hmm_parameters & M;
	genotype * G;
	bitmatrix Hhap, Hvar;

	//COORDINATES & CONSTANTS
	int segment_first;
	int segment_last;
	int locus_first;
	int locus_last;
	int ambiguous_first;
	int ambiguous_last;
	int missing_first;
	int missing_last;
	int transition_first;
	int transition_last;
	unsigned int n_cond_haps;
	unsigned int n_missing;

	//CURSORS
	int curr_segment_index;
	int curr_segment_locus;
	int curr_abs_locus;
	int prev_abs_locus;
	int curr_rel_locus;
	int curr_rel_locus_offset;
	int curr_abs_ambiguous;
	int curr_abs_transition;
	int curr_abs_missing;
	int curr_rel_missing;


	//DYNAMIC ARRAYS
	float probSumT;
	aligned_vector32 < float > prob;
	aligned_vector32 < float > probSumK;
	aligned_vector32 < float > probSumH;
	std::vector < aligned_vector32 < float > > Alpha;
	std::vector < aligned_vector32 < float > > AlphaSum;
	std::vector < int > AlphaLocus;
	aligned_vector32 < float > AlphaSumSum;
	std::vector < aligned_vector32 < float > > AlphaMissing;
	std::vector < aligned_vector32 < float > > AlphaSumMissing;
	float HProbs [HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));
	double DProbs [HAP_NUMBER * HAP_NUMBER * HAP_NUMBER * HAP_NUMBER] __attribute__ ((aligned(32)));

	//STATIC ARRAYS
	float sumHProbs;
	double sumDProbs;
	float g0[HAP_NUMBER], g1[HAP_NUMBER];
	float nt, yt;

	//SUPER-SITE SUPPORT
	const std::vector<SuperSite>* super_sites;
	const std::vector<bool>* is_super_site;
	const std::vector<int>* locus_to_super_idx;
	const uint8_t* panel_codes;
	const std::vector<int>* super_site_var_index;
	const std::vector<unsigned int>* cond_idx;
	aligned_vector32<uint8_t> ss_cond_codes;
	aligned_vector32<float> ss_emissions;

	//INLINED AND UNROLLED ROUTINES
	void INIT_HOM();
	void INIT_AMB();
	void INIT_MIS();
	bool RUN_HOM(char);
	void RUN_AMB();
	void RUN_MIS();
	void COLLAPSE_HOM();
	void COLLAPSE_AMB();
	void COLLAPSE_MIS();
	void SUMK();
	void IMPUTE(std::vector < float > & );
	bool TRANS_HAP();
	bool TRANS_DIP_MULT();
	bool TRANS_DIP_ADD();
	void SET_FIRST_TRANS(std::vector < double > & );
	int SET_OTHER_TRANS(std::vector < double > & );

public:
	//CONSTRUCTOR/DESTRUCTOR
	haplotype_segment_single(genotype *, bitmatrix &, std::vector < unsigned int > &, window &, hmm_parameters &,
		const std::vector<SuperSite>* _super_sites = nullptr,
		const std::vector<bool>* _is_super_site = nullptr,
		const std::vector<int>* _locus_to_super_idx = nullptr,
		const uint8_t* _panel_codes = nullptr,
		const std::vector<int>* _super_site_var_index = nullptr);
	~haplotype_segment_single();

	//void fetch();
	void forward();
	int backward(std::vector < double > &, std::vector < float > &);
};

/*******************************************************************************/
/*****************			HOMOZYGOUS GENOTYPE			************************/
/*******************************************************************************/


inline
void haplotype_segment_single::INIT_HOM() {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { INIT_MIS(); return; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _emit;
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return;
    }
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _mm256_set1_ps((ag==ah)?1.0f:M.ed/M.ee);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
bool haplotype_segment_single::RUN_HOM(char rare_allele) {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { RUN_MIS(); return true; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
        __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
        _tFreq = _mm256_mul_ps(_tFreq, _factor);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            _prob = _mm256_mul_ps(_prob, _emit);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return true;
    }
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    if (rare_allele < 0 || ag == rare_allele) {
        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
        __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
        _tFreq = _mm256_mul_ps(_tFreq, _factor);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
		__m256 _mismatch = _mm256_set1_ps(M.ed/M.ee);
		for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
			bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
			__m256 _prob = _mm256_load_ps(&prob[i]);
			_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
			if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
			_sum = _mm256_add_ps(_sum, _prob);
			_mm256_store_ps(&prob[i], _prob);
		}
		_mm256_store_ps(&probSumH[0], _sum);
		probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
		return true;
	}
	return false;
}

inline
void haplotype_segment_single::COLLAPSE_HOM() {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { COLLAPSE_MIS(); return; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            _prob = _mm256_mul_ps(_prob, _emit);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return;
    }
    bool ag = VAR_GET_HAP0(MOD2(curr_abs_locus), G->Variants[DIV2(curr_abs_locus)]);
    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);					//Check divide by probSumT here!
    __m256 _nt = _mm256_set1_ps(nt / probSumT);
    __m256 _mismatch = _mm256_set1_ps(M.ed/M.ee);
    for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
        bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
        __m256 _prob = _mm256_set1_ps(probSumK[k]);
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
        if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
        _sum = _mm256_add_ps(_sum, _prob);
        _mm256_store_ps(&prob[i], _prob);
    }
    _mm256_store_ps(&probSumH[0], _sum);
    probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			HETEROZYGOUS GENOTYPE			********************/
/*******************************************************************************/

inline
void haplotype_segment_single::INIT_AMB() {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { INIT_MIS(); return; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _emit;
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return;
    }
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
    }
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _emit[ah];
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_single::RUN_AMB() {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { RUN_MIS(); return; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
        __m256 _tFreq = _mm256_load_ps(&probSumH[0]);
        _tFreq = _mm256_mul_ps(_tFreq, _factor);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            _prob = _mm256_mul_ps(_prob, _emit);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return;
    }
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
    }
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
	__m256 _tFreq = _mm256_load_ps(&probSumH[0]);
	_tFreq = _mm256_mul_ps(_tFreq, _factor);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	__m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_load_ps(&prob[i]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_prob = _mm256_mul_ps(_prob, _emit[ah]);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*
inline
void haplotype_segment_single::RUN_AMB() {
	unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
	for (int h = 0 ; h < HAP_NUMBER ; h ++) {
		g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
		g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
	}
	__m256d _sum0 = _mm256_set1_pd(0.0f);
	__m256d _sum1 = _mm256_set1_pd(0.0f);
	__m256d _factor = _mm256_set1_pd(yt / (n_cond_haps * probSumT));
	__m256d _tFreq0 = _mm256_load_pd(&probSumH[0]);
	__m256d _tFreq1 = _mm256_load_pd(&probSumH[4]);
	_tFreq0 = _mm256_mul_pd(_tFreq0, _factor);
	_tFreq1 = _mm256_mul_pd(_tFreq1, _factor);
	__m256d _nt = _mm256_set1_pd(nt / probSumT);
	__m256d _emit0[2], _emit1[2];
	_emit0[0] = _mm256_loadu_pd(&g0[0]);
	_emit0[1] = _mm256_loadu_pd(&g1[0]);
	_emit1[0] = _mm256_loadu_pd(&g0[4]);
	_emit1[1] = _mm256_loadu_pd(&g1[4]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256d _prob0 = _mm256_load_pd(&prob[i+0]);
		__m256d _prob1 = _mm256_load_pd(&prob[i+4]);
		_prob0 = _mm256_fmadd_pd(_prob0, _nt, _tFreq0);
		_prob1 = _mm256_fmadd_pd(_prob1, _nt, _tFreq1);
		_prob0 = _mm256_mul_pd(_prob0, _emit0[ah]);
		_prob1 = _mm256_mul_pd(_prob1, _emit1[ah]);
		_sum0 = _mm256_add_pd(_sum0, _prob0);
		_sum1 = _mm256_add_pd(_sum1, _prob1);
		_mm256_store_pd(&prob[i+0], _prob0);
		_mm256_store_pd(&prob[i+4], _prob1);
	}
	_mm256_store_pd(&probSumH[0], _sum0);
	_mm256_store_pd(&probSumH[4], _sum1);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}
*/

inline
void haplotype_segment_single::COLLAPSE_AMB() {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        uint8_t sample_code = getSampleSuperSiteAlleleCode(G, ss, *super_site_var_index, 0);
        if (sample_code == SUPERSITE_CODE_MISSING) { COLLAPSE_MIS(); return; }
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        precomputeSuperSiteEmissions_FloatScalar(ss_cond_codes.data(), n_cond_haps, sample_code, 1.0f, (float)(M.ed/M.ee), ss_emissions);
        __m256 _sum = _mm256_set1_ps(0.0f);
        __m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
        __m256 _nt = _mm256_set1_ps(nt / probSumT);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            __m256 _emit = _mm256_set1_ps(ss_emissions[k]);
            __m256 _prob = _mm256_set1_ps(probSumK[k]);
            _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
            _prob = _mm256_mul_ps(_prob, _emit);
            _sum = _mm256_add_ps(_sum, _prob);
            _mm256_store_ps(&prob[i], _prob);
        }
        _mm256_store_ps(&probSumH[0], _sum);
        probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
        return;
    }
    unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
    for (int h = 0 ; h < HAP_NUMBER ; h ++) {
        g0[h] = HAP_GET(amb_code,h)?M.ed/M.ee:1.0f;
        g1[h] = HAP_GET(amb_code,h)?1.0f:M.ed/M.ee;
    }
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	__m256 _emit[2]; _emit[0] = _mm256_loadu_ps(&g0[0]); _emit[1] = _mm256_loadu_ps(&g1[0]);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
		__m256 _prob = _mm256_set1_ps(probSumK[k]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_prob = _mm256_mul_ps(_prob, _emit[ah]);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************			MISSING GENOTYPE			************************/
/*******************************************************************************/

inline
void haplotype_segment_single::INIT_MIS() {
	fill(prob.begin(), prob.end(), 1.0f/(HAP_NUMBER * n_cond_haps));
	fill(probSumH.begin(), probSumH.end(), 1.0f/HAP_NUMBER);
	probSumT = 1.0f;
}

inline
void haplotype_segment_single::RUN_MIS() {
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _factor = _mm256_set1_ps(yt / (n_cond_haps * probSumT));
	__m256 _tFreq = _mm256_load_ps(&probSumH[0]);
	_tFreq = _mm256_mul_ps(_tFreq, _factor);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256 _prob = _mm256_load_ps(&prob[i]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

inline
void haplotype_segment_single::COLLAPSE_MIS() {
	__m256 _sum = _mm256_set1_ps(0.0f);
	__m256 _tFreq = _mm256_set1_ps(yt / n_cond_haps);
	__m256 _nt = _mm256_set1_ps(nt / probSumT);
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		__m256 _prob = _mm256_set1_ps(probSumK[k]);
		_prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
		_sum = _mm256_add_ps(_sum, _prob);
		_mm256_store_ps(&prob[i], _prob);
	}
	_mm256_store_ps(&probSumH[0], _sum);
	probSumT = probSumH[0] + probSumH[1] + probSumH[2] + probSumH[3] + probSumH[4] + probSumH[5] + probSumH[6] + probSumH[7];
}

/*******************************************************************************/
/*****************					SUM Ks				************************/
/*******************************************************************************/

inline
void haplotype_segment_single::SUMK() {
	for(int k = 0, i = 0 ; k != n_cond_haps ; ++k, i += HAP_NUMBER) {
		probSumK[k] = prob[i+0] + prob[i+1] + prob[i+2] + prob[i+3] + prob[i+4] + prob[i+5] + prob[i+6] + prob[i+7];
	}
}

/*******************************************************************************/
/*****************		TRANSITION COMPUTATIONS			************************/
/*******************************************************************************/

inline
bool haplotype_segment_single::TRANS_HAP() {
	sumHProbs = 0.0f;
	unsigned int  curr_rel_segment_index = curr_segment_index-segment_first;
	yt = M.getForwardTransProb(AlphaLocus[curr_rel_segment_index - 1], prev_abs_locus);
	nt = 1.0f - yt;
	float fact1 = nt / AlphaSumSum[curr_rel_segment_index - 1];
	for (int h1 = 0 ; h1 < HAP_NUMBER ; h1++) {
		__m256 _sum = _mm256_set1_ps(0.0f);
		float fact2 = (AlphaSum[curr_rel_segment_index-1][h1]/AlphaSumSum[curr_rel_segment_index-1]) * yt / n_cond_haps;
		for (int k = 0 ; k < n_cond_haps ; k ++) {
			__m256 _alpha = _mm256_set1_ps(Alpha[curr_rel_segment_index-1][k*HAP_NUMBER + h1] * fact1 + fact2);
			__m256 _beta = _mm256_load_ps(&prob[k*HAP_NUMBER]);
			_sum = _mm256_add_ps(_sum, _mm256_mul_ps(_alpha, _beta));
		}
		_mm256_store_ps(&HProbs[h1*HAP_NUMBER], _sum);
		sumHProbs += HProbs[h1*HAP_NUMBER+0]+HProbs[h1*HAP_NUMBER+1]+HProbs[h1*HAP_NUMBER+2]+HProbs[h1*HAP_NUMBER+3]+HProbs[h1*HAP_NUMBER+4]+HProbs[h1*HAP_NUMBER+5]+HProbs[h1*HAP_NUMBER+6]+HProbs[h1*HAP_NUMBER+7];
	}
	return (std::isnan(sumHProbs) || std::isinf(sumHProbs) || sumHProbs < std::numeric_limits<float>::min());
}

inline
bool haplotype_segment_single::TRANS_DIP_MULT() {
	sumDProbs= 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) * ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

inline
bool haplotype_segment_single::TRANS_DIP_ADD() {
	sumDProbs = 0.0f;
	double scaling = 1.0 / sumHProbs;
	for (int pd = 0, t = 0 ; pd < 64 ; ++pd) {
		if (DIP_GET(G->Diplotypes[curr_segment_index-1], pd)) {
			for (int nd = 0 ; nd < 64 ; ++nd) {
				if (DIP_GET(G->Diplotypes[curr_segment_index], nd)) {
					DProbs[t] = DProbs[t] = (((double)HProbs[DIP_HAP0(pd)*HAP_NUMBER+DIP_HAP0(nd)]) * scaling) + ((double)(HProbs[DIP_HAP1(pd)*HAP_NUMBER+DIP_HAP1(nd)]) * scaling);
					sumDProbs += DProbs[t];
					t++;
				}
			}
		}
	}
	return (std::isnan(sumDProbs) || std::isinf(sumDProbs) || sumDProbs < std::numeric_limits<double>::min());
}

inline
void haplotype_segment_single::IMPUTE(std::vector < float > & missing_probabilities) {
    int ss_idx = -1;
    if (super_sites && locus_to_super_idx) ss_idx = (*locus_to_super_idx)[curr_abs_locus];

    __m256 _sum = _mm256_set1_ps(0.0f);
    __m256 _alphaSum = _mm256_load_ps(&AlphaSumMissing[curr_rel_missing][0]);
    __m256 _ones = _mm256_set1_ps(1.0f);
    _alphaSum = _mm256_div_ps(_ones, _alphaSum);

    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Determine class for current split variant (0=REF, 1..n_alts)
        int target_class = 0;
        for (uint8_t ai = 0; ai < ss.var_count; ++ai) {
            if ((*super_site_var_index)[ss.var_start + ai] == curr_abs_locus) { target_class = (int)ai + 1; break; }
        }
        // Unpack cond hap codes
        for (int k = 0; k < (int)n_cond_haps; ++k) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        // Prepare accumulators per class
        __m256 sum_classes[SUPERSITE_MAX_ALTS + 1];
        for (int c = 0; c <= (int)ss.var_count; ++c) sum_classes[c] = _mm256_set1_ps(0.0f);
        __m256 denom = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            int code = (int)ss_cond_codes[k]; if (code > (int)ss.var_count) code = 0;
            __m256 _prob = _mm256_load_ps(&prob[i]);
            __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
            _sum = _mm256_mul_ps(_mm256_mul_ps(_alpha, _alphaSum), _prob);
            sum_classes[code] = _mm256_add_ps(sum_classes[code], _sum);
            denom = _mm256_add_ps(denom, _sum);
        }
        __m256 p = sum_classes[target_class];
        alignas(32) float pv[8], dv[8];
        _mm256_store_ps(pv, p);
        _mm256_store_ps(dv, denom);
        for (int h = 0; h < HAP_NUMBER; ++h) {
            float d = dv[h];
            missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = (d > 0.0f) ? (pv[h] / d) : 0.0f;
        }
    } else {
        __m256 _sumA[2]; _sumA[0] = _mm256_set1_ps(0.0f); _sumA[1] = _mm256_set1_ps(0.0f);
        for (int k = 0, i = 0; k != (int)n_cond_haps; ++k, i += HAP_NUMBER) {
            bool ah = Hvar.get(curr_rel_locus+curr_rel_locus_offset, k);
            __m256 _prob = _mm256_load_ps(&prob[i]);
            __m256 _alpha = _mm256_load_ps(&AlphaMissing[curr_rel_missing][i]);
            _sum = _mm256_mul_ps(_mm256_mul_ps(_alpha, _alphaSum), _prob);
            _sumA[ah] = _mm256_add_ps(_sumA[ah], _sum);
        }
        float* prob0 = (float*)&_sumA[0];
        float* prob1 = (float*)&_sumA[1];
        for (int h = 0; h < HAP_NUMBER; h++) missing_probabilities[curr_abs_missing * HAP_NUMBER + h] = prob1[h] / (prob0[h] + prob1[h]);
    }
}

#endif

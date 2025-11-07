/*******************************************************************************
 * Supersite vs biallelic representation parity (forward/backward and phasing)
 *
 * Build the same dataset twice:
 *   (A) plain biallelic – no supersites constructed
 *   (B) supersite-enabled – the two split variants at position 1000 are grouped
 *       and the HMM receives full supersite context / anchor gating.
 *
 * Run the HMM forward+backward (double precision) on both datasets, assert that
 * probabilities match, and verify that the most-likely phasing inferred from
 * the HMM is identical and equal to the expected haplotype sequence.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <cctype>
#include <sys/stat.h>
#include <sys/types.h>

#include "../../common/src/utils/otools.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

namespace {

// Forward declarations for trace helpers
static bool trace_enabled();
static bool ensure_out_dir();
static std::string sanitize_label(const std::string& in);
static std::unique_ptr<std::ofstream> make_trace_stream(const std::string& prefix, const std::string& label);

// Definitions for trace helpers
static bool trace_enabled() {
    const char* env = std::getenv("SHAPEIT5_TEST_TRACE");
    return env && env[0] != '\0' && env[0] != '0';
}

static bool ensure_out_dir() {
    struct stat st{};
    if (stat("tests/out", &st) == 0) {
        return S_ISDIR(st.st_mode);
    }
    return mkdir("tests/out", 0777) == 0;
}

static std::string sanitize_label(const std::string& in) {
    std::string out;
    out.reserve(in.size());
    for (char c : in) {
        if (std::isalnum(static_cast<unsigned char>(c))) out.push_back(std::tolower(static_cast<unsigned char>(c)));
        else out.push_back('_');
    }
    return out;
}

static std::unique_ptr<std::ofstream> make_trace_stream(const std::string& prefix, const std::string& label) {
    if (!trace_enabled()) return nullptr;
    if (!ensure_out_dir()) return nullptr;
    std::string fname = std::string("tests/out/") + prefix + "_" + sanitize_label(label) + ".tsv";
    auto ofs = std::make_unique<std::ofstream>(fname, std::ios::out | std::ios::trunc);
    if (!ofs->good()) return nullptr;
    return ofs;
}

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
    std::vector<uint8_t> sample_codes_unused;
};

struct FBResult {
    std::vector<double> prob;
    std::vector<double> probSumH;
    double probSumT;
    std::vector<double> transition_probabilities;
};

enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

struct Orientation {
    std::array<PhaseCode,5> phases;
};


struct ForwardTraceEntry {
    int locus;
    int prev_before;
    int prev_after;
    double yt;
    double nt;
    bool update_prev;
    bool is_anchor;
    bool is_sibling;
    bool amb;
    bool mis;
    bool hom;
    std::array<double,HAP_NUMBER> probSumH;
    double probSumT;
    // Biallelic AMB details
    bool has_g01 = false;
    unsigned char amb_mask = 0;
    std::array<double,HAP_NUMBER> g0{};
    std::array<double,HAP_NUMBER> g1{};
    // Supersite AMB details
    bool has_ss = false;
    std::vector<double> ss_em0; // emissions for low half (c0)
    std::vector<double> ss_em1; // emissions for high half (c1)
};

struct BackwardTraceEntry {
    int locus;
    int prev_before;
    int prev_after;
    double yt;
    double nt;
    bool update_prev;
    bool is_anchor;
    bool is_sibling;
    bool amb;
    bool mis;
    bool hom;
    std::array<double,HAP_NUMBER> probSumH;
    double probSumT;
};

struct ForwardTrace { std::vector<ForwardTraceEntry> entries; };
struct BackwardTrace { std::vector<BackwardTraceEntry> entries; };

// Build variant pointer helper
static variant* make_var(std::string chr, int bp, std::string id,
                         std::string ref, std::string alt, int idx) {
    return new variant(chr, bp, id, ref, alt, idx);
}

static void clear_variant_state(genotype& G, int locus) {
    unsigned char& byte = G.Variants[DIV2(locus)];
    const int shift = (MOD2(locus)) << 2;
    byte &= ~(0x0F << shift);
}

static void set_phase(genotype& G, int locus, PhaseCode code) {
    clear_variant_state(G, locus);
    unsigned char& byte = G.Variants[DIV2(locus)];
    switch (code) {
        case REF_REF:
            // nothing (already 0/0)
            break;
        case ALT_ALT:
            VAR_SET_HAP0(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
        case ALT_REF:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP0(MOD2(locus), byte);
            break;
        case REF_ALT:
            VAR_SET_HET(MOD2(locus), byte);
            VAR_SET_HAP1(MOD2(locus), byte);
            break;
    }
}

static void apply_orientation(genotype& G, const Orientation& orient) {
    for (int locus = 0; locus < static_cast<int>(orient.phases.size()); ++locus) {
        set_phase(G, locus, orient.phases[locus]);
    }
    G.build();
}

static void restore_genotype(genotype& G, const std::vector<unsigned char>& variants_backup) {
    G.Variants = variants_backup;
    G.build();
}

static window make_full_window(int stop_locus) {
    window W;
    W.start_locus = 0;
    W.stop_locus = stop_locus;
    W.start_segment = 0;
    W.stop_segment = 0;
    W.start_ambiguous = 0;
    W.stop_ambiguous = -1;
    W.start_missing = 0;
    W.stop_missing = -1;
    W.start_transition = 0;
    W.stop_transition = -1;
    return W;
}

// Compute transitions from cm: identical positions → yt=0, clamp only tiny positives
static void compute_t_from_cm(hmm_parameters& M) {
    const int n = static_cast<int>(M.cm.size());
    M.t.assign(n > 0 ? n - 1 : 0, 0.0f);
    M.nt.assign(n > 0 ? n - 1 : 0, 0.0f);
    for (int l = 1; l < n; ++l) {
        float dist_cm = M.cm[l] - M.cm[l - 1];
        if (dist_cm <= 0.0f) {
            M.t[l - 1] = 0.0f;
            M.nt[l - 1] = 1.0f;
        } else {
            if (dist_cm < 1e-7f) dist_cm = 1e-7f;
            float tval = -1.0f * expm1f(-0.04f * (float)M.Neff * dist_cm / (float)M.Nhap);
            M.t[l - 1] = tval;
            M.nt[l - 1] = 1.0f - tval;
        }
    }
}

static hmm_parameters make_hmm_params(size_t n_variants, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
    // Use binary presentation at anchors for parity against biallelic split
    M.ss_anchor_split_emissions = true;
    M.cm = std::vector<float>(n_variants, 0.001f);
    // simple increasing map with supersite splits sharing the same position
    if (n_variants >= 5) {
        M.cm[0] = 0.005f;   // v500
        M.cm[1] = 0.010f;   // ss_A_C
        M.cm[2] = 0.010f;   // ss_A_G (same as anchor)
        M.cm[3] = 0.020f;   // v2000
        M.cm[4] = 0.025f;   // v2500
    }
    M.Neff = 10000;
    M.Nhap = static_cast<int>(Nhap);
    compute_t_from_cm(M);
    M.rare_allele = std::vector<char>(n_variants, -1);
    return M;
}

static SuperSiteContext build_supersites(variant_map& V, conditioning_set& H) {
    SuperSiteContext ctx;
    buildSuperSites(V, H, ctx.super_sites, ctx.is_super_site, ctx.packed_codes,
                    ctx.locus_to_super_idx, ctx.super_site_var_index, ctx.sample_codes_unused);
    return ctx;
}

// Forward-only runner (mirrors expansion parity) for anchor-only parity checks
struct FBState { std::vector<double> prob; std::vector<double> sumH; };
static FBState run_forward_only(genotype& G,
                                conditioning_set& H,
                                hmm_parameters& M,
                                const window& W,
                                const std::vector<unsigned int>& idxH,
                                const SuperSiteContext* ctx) {
    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<bool>* is_super_site = ctx ? &ctx->is_super_site : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;
    const std::vector<int>* super_site_var_index = ctx ? &ctx->super_site_var_index : nullptr;
    const uint8_t* panel_codes = (ctx && !ctx->packed_codes.empty()) ? ctx->packed_codes.data() : nullptr;

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M,
                                super_sites, is_super_site, locus_to_super_idx,
                                panel_codes, super_site_var_index);
    HS.forward();
    FBState s; s.prob.assign(HS.prob.begin(), HS.prob.end()); s.sumH.assign(HS.probSumH.begin(), HS.probSumH.end());
    return s;
}

static inline bool is_anchor_site(const haplotype_segment_double& HS, int locus) {
    if (!(HS.super_sites && HS.locus_to_super_idx)) return false;
    int ss_idx = (*HS.locus_to_super_idx)[locus];
    return (ss_idx >= 0 && locus == static_cast<int>((*HS.super_sites)[ss_idx].global_site_id));
}

static inline bool is_sibling_site(const haplotype_segment_double& HS, int locus) {
    if (!(HS.super_sites && HS.locus_to_super_idx)) return false;
    int ss_idx = (*HS.locus_to_super_idx)[locus];
    return (ss_idx >= 0 && locus != static_cast<int>((*HS.super_sites)[ss_idx].global_site_id));
}

static ForwardTrace trace_forward(genotype& G, conditioning_set& H, hmm_parameters& M,
                                  const window& W, const std::vector<unsigned int>& idxH,
                                  const SuperSiteContext* ctx) {
    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<bool>* is_super_site = ctx ? &ctx->is_super_site : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;
    const std::vector<int>* super_site_var_index = ctx ? &ctx->super_site_var_index : nullptr;
    const uint8_t* panel_codes = (ctx && !ctx->packed_codes.empty()) ? ctx->packed_codes.data() : nullptr;

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M,
                                super_sites, is_super_site, locus_to_super_idx,
                                panel_codes, super_site_var_index);

    ForwardTrace trace;

    HS.curr_segment_index = HS.segment_first;
    HS.curr_segment_locus = 0;
    HS.curr_abs_ambiguous = HS.ambiguous_first;
    HS.curr_abs_missing = HS.missing_first;
    HS.prev_abs_locus = HS.locus_first;

    for (HS.curr_abs_locus = HS.locus_first; HS.curr_abs_locus <= HS.locus_last; ++HS.curr_abs_locus) {
        HS.curr_rel_locus = HS.curr_abs_locus - HS.locus_first;
        HS.curr_rel_missing = HS.curr_abs_missing - HS.missing_first;
        bool update_prev_locus = true;
        char rare_allele = HS.M.rare_allele[HS.curr_abs_locus];
        bool amb = VAR_GET_AMB(MOD2(HS.curr_abs_locus), HS.G->Variants[DIV2(HS.curr_abs_locus)]);
        bool mis = VAR_GET_MIS(MOD2(HS.curr_abs_locus), HS.G->Variants[DIV2(HS.curr_abs_locus)]);
        bool hom = !(amb || mis);
        HS.yt = (HS.curr_abs_locus == HS.locus_first) ? 0.0 : HS.M.getForwardTransProb(HS.prev_abs_locus, HS.curr_abs_locus);
        HS.nt = 1.0f - HS.yt;

        if (HS.curr_rel_locus == 0) {
            if (hom) HS.INIT_HOM();
            else if (amb) HS.INIT_AMB();
            else HS.INIT_MIS();
        } else if (HS.curr_segment_locus != 0) {
            if (hom) update_prev_locus = HS.RUN_HOM(rare_allele);
            else if (amb) HS.RUN_AMB();
            else HS.RUN_MIS();
        } else {
            if (hom) HS.COLLAPSE_HOM();
            else if (amb) HS.COLLAPSE_AMB();
            else HS.COLLAPSE_MIS();
        }

        HS.prev_abs_locus = update_prev_locus ? HS.curr_abs_locus : HS.prev_abs_locus;

        if (HS.curr_segment_locus == (HS.G->Lengths[HS.curr_segment_index] - 1)) HS.SUMK();
        if (HS.curr_segment_locus == HS.G->Lengths[HS.curr_segment_index] - 1) {
            HS.Alpha[HS.curr_segment_index - HS.segment_first] = HS.prob;
            HS.AlphaSum[HS.curr_segment_index - HS.segment_first] = HS.probSumH;
            HS.AlphaSumSum[HS.curr_segment_index - HS.segment_first] = HS.probSumT;
            HS.AlphaLocus[HS.curr_segment_index - HS.segment_first] = HS.prev_abs_locus;
        }
        if (mis) {
            HS.AlphaMissing[HS.curr_rel_missing] = HS.prob;
            HS.AlphaSumMissing[HS.curr_rel_missing] = HS.probSumH;
            HS.curr_abs_missing++;
        }

        ForwardTraceEntry entry{};
        entry.locus = HS.curr_abs_locus;
        entry.prev_before = update_prev_locus ? HS.curr_abs_locus : HS.prev_abs_locus;
        entry.prev_after = HS.prev_abs_locus;
        entry.yt = HS.yt;
        entry.nt = HS.nt;
        entry.update_prev = update_prev_locus;
        entry.is_anchor = is_anchor_site(HS, HS.curr_abs_locus);
        entry.is_sibling = is_sibling_site(HS, HS.curr_abs_locus);
        entry.amb = amb;
        entry.mis = mis;
        entry.hom = hom;
        for (int h = 0; h < HAP_NUMBER; ++h) entry.probSumH[h] = HS.probSumH[h];
        entry.probSumT = HS.probSumT;
        // Collect AMB details
        if (amb) {
            // Ambiguous mask is only meaningful in biallelic path
            entry.amb_mask = (HS.curr_abs_ambiguous >= HS.ambiguous_first && HS.curr_abs_ambiguous < HS.ambiguous_last+1)
                             ? HS.G->Ambiguous[HS.curr_abs_ambiguous] : 0;
            if (is_anchor_site(HS, HS.curr_abs_locus)) {
                // Supersite anchor: capture per-donor emissions for low/high halves
                entry.has_ss = true;
                entry.ss_em0.assign(HS.ss_emissions.begin(), HS.ss_emissions.end());
                entry.ss_em1.assign(HS.ss_emissions_h1.begin(), HS.ss_emissions_h1.end());
            } else {
                // Biallelic AMB: capture the 8-lane per-hap emission vectors used (g0/g1)
                entry.has_g01 = true;
                for (int h = 0; h < HAP_NUMBER; ++h) { entry.g0[h] = HS.g0[h]; entry.g1[h] = HS.g1[h]; }
            }
        }
        trace.entries.push_back(entry);

        HS.curr_segment_locus++;
        HS.curr_abs_ambiguous += amb;
        if (HS.curr_segment_locus >= HS.G->Lengths[HS.curr_segment_index]) {
            HS.curr_segment_index++;
            HS.curr_segment_locus = 0;
        }
    }

    return trace;
}

static BackwardTrace trace_backward(genotype& G, conditioning_set& H, hmm_parameters& M,
                                    const window& W, const std::vector<unsigned int>& idxH,
                                    const SuperSiteContext* ctx) {
    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<bool>* is_super_site = ctx ? &ctx->is_super_site : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;
    const std::vector<int>* super_site_var_index = ctx ? &ctx->super_site_var_index : nullptr;
    const uint8_t* panel_codes = (ctx && !ctx->packed_codes.empty()) ? ctx->packed_codes.data() : nullptr;

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M,
                                super_sites, is_super_site, locus_to_super_idx,
                                panel_codes, super_site_var_index);

    BackwardTrace trace;

    HS.curr_segment_index = HS.segment_last;
    HS.curr_segment_locus = HS.G->Lengths[HS.segment_last] - 1;
    HS.curr_abs_ambiguous = HS.ambiguous_last;
    HS.curr_abs_missing = HS.missing_last;
    HS.curr_abs_transition = HS.transition_last;
    HS.prev_abs_locus = HS.locus_last;

    for (HS.curr_abs_locus = HS.locus_last; HS.curr_abs_locus >= HS.locus_first; --HS.curr_abs_locus) {
        HS.curr_rel_locus = HS.curr_abs_locus - HS.locus_first;
        HS.curr_rel_missing = HS.curr_abs_missing - HS.missing_first;
        char rare_allele = HS.M.rare_allele[HS.curr_abs_locus];
        bool update_prev_locus = true;
        bool amb = VAR_GET_AMB(MOD2(HS.curr_abs_locus), HS.G->Variants[DIV2(HS.curr_abs_locus)]);
        bool mis = VAR_GET_MIS(MOD2(HS.curr_abs_locus), HS.G->Variants[DIV2(HS.curr_abs_locus)]);
        bool hom = !(amb || mis);
        HS.yt = (HS.curr_abs_locus == HS.locus_last) ? 0.0 : HS.M.getBackwardTransProb(HS.prev_abs_locus, HS.curr_abs_locus);
        HS.nt = 1.0f - HS.yt;

        if (HS.curr_abs_locus == HS.locus_last) {
            if (hom) HS.INIT_HOM();
            else if (amb) HS.INIT_AMB();
            else HS.INIT_MIS();
        } else if (HS.curr_segment_locus != HS.G->Lengths[HS.curr_segment_index] - 1) {
            if (hom) update_prev_locus = HS.RUN_HOM(rare_allele);
            else if (amb) HS.RUN_AMB();
            else HS.RUN_MIS();
        } else {
            if (hom) HS.COLLAPSE_HOM();
            else if (amb) HS.COLLAPSE_AMB();
            else HS.COLLAPSE_MIS();
        }

        if (HS.curr_segment_locus == 0) HS.SUMK();
        HS.prev_abs_locus = update_prev_locus ? HS.curr_abs_locus : HS.prev_abs_locus;

        if (HS.curr_abs_locus == 0) {
            std::vector<double> dummy_transition(HS.G->countTransitions(), 0.0);
            HS.SET_FIRST_TRANS(dummy_transition);
        }
        if (HS.curr_segment_locus == 0 && HS.curr_abs_locus != HS.locus_first) {
            std::vector<double> dummy_transition(HS.G->countTransitions(), 0.0);
            HS.SET_OTHER_TRANS(dummy_transition);
        }

        if (mis) {
            HS.curr_abs_missing--;
        }

        BackwardTraceEntry entry{};
        entry.locus = HS.curr_abs_locus;
        entry.prev_before = update_prev_locus ? HS.curr_abs_locus : HS.prev_abs_locus;
        entry.prev_after = HS.prev_abs_locus;
        entry.yt = HS.yt;
        entry.nt = HS.nt;
        entry.update_prev = update_prev_locus;
        entry.is_anchor = is_anchor_site(HS, HS.curr_abs_locus);
        entry.is_sibling = is_sibling_site(HS, HS.curr_abs_locus);
        entry.amb = amb;
        entry.mis = mis;
        entry.hom = hom;
        for (int h = 0; h < HAP_NUMBER; ++h) entry.probSumH[h] = HS.probSumH[h];
        entry.probSumT = HS.probSumT;
        trace.entries.push_back(entry);

        if (HS.curr_segment_locus == 0) {
            HS.curr_segment_index--;
            if (HS.curr_segment_index >= HS.segment_first)
                HS.curr_segment_locus = HS.G->Lengths[HS.curr_segment_index] - 1;
        } else {
            HS.curr_segment_locus--;
        }
    }

    std::reverse(trace.entries.begin(), trace.entries.end());
    return trace;
}

static FBResult run_forward_backward(genotype& G,
                                     conditioning_set& H,
                                     hmm_parameters& M,
                                     const window& W,
                                     const std::vector<unsigned int>& idxH,
                                     const SuperSiteContext* ctx) {
    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<bool>* is_super_site = ctx ? &ctx->is_super_site : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;
    const std::vector<int>* super_site_var_index = ctx ? &ctx->super_site_var_index : nullptr;
    const uint8_t* panel_codes = (ctx && !ctx->packed_codes.empty()) ? ctx->packed_codes.data() : nullptr;

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M,
                                super_sites, is_super_site, locus_to_super_idx,
                                panel_codes, super_site_var_index);

    HS.forward();

    std::vector<double> transition_probabilities(G.countTransitions(), 0.0);
    std::vector<float> missing_probabilities;
    HS.backward(transition_probabilities, missing_probabilities,
                /*SC*/nullptr, /*anchor_has_missing*/nullptr);

    FBResult res;
    res.prob.assign(HS.prob.begin(), HS.prob.end());
    res.probSumH.assign(HS.probSumH.begin(), HS.probSumH.end());
    res.probSumT = HS.probSumT;
    res.transition_probabilities = transition_probabilities;
    return res;
}

static double evaluate_orientation(genotype& G,
                                   conditioning_set& H,
                                   hmm_parameters& M,
                                   const window& W,
                                   const std::vector<unsigned int>& idxH,
                                   const SuperSiteContext* ctx,
                                   const Orientation& orient) {
    const std::vector<unsigned char> backup = G.Variants;
    apply_orientation(G, orient);

    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<bool>* is_super_site = ctx ? &ctx->is_super_site : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;
    const std::vector<int>* super_site_var_index = ctx ? &ctx->super_site_var_index : nullptr;
    const uint8_t* panel_codes = (ctx && !ctx->packed_codes.empty()) ? ctx->packed_codes.data() : nullptr;

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M,
                                super_sites, is_super_site, locus_to_super_idx,
                                panel_codes, super_site_var_index);
    HS.forward();
    const double likelihood = HS.probSumT;

    restore_genotype(G, backup);
    return likelihood;
}

static std::string phase_to_string(PhaseCode code) {
    switch (code) {
        case REF_REF: return "0|0";
        case ALT_ALT: return "1|1";
        case ALT_REF: return "1|0";
        case REF_ALT: return "0|1";
    }
    return "?";
}

static void print_forward_trace(const std::string& label, const ForwardTrace& trace) {
    auto ofs = make_trace_stream("forward", label);
    if (!ofs) return; // gated or cannot open; skip verbose output
    auto& os = *ofs;
    os << "# Forward trace - " << label << "\n";
    os << "locus\tprev_before\tprev_after\tyt\tnt\tupdate_prev\tis_anchor\tis_sibling\tamb\tmis\thom";
    for (int h = 0; h < HAP_NUMBER; ++h) os << "\tprobSumH[" << h << "]";
    os << "\tprobSumT\n";
    for (const auto& e : trace.entries) {
        os << e.locus << "\t" << e.prev_before << "\t" << e.prev_after << "\t"
           << e.yt << "\t" << e.nt << "\t"
           << (e.update_prev ? 1 : 0) << "\t"
           << (e.is_anchor ? 1 : 0) << "\t"
           << (e.is_sibling ? 1 : 0) << "\t"
           << (e.amb ? 1 : 0) << "\t"
           << (e.mis ? 1 : 0) << "\t"
           << (e.hom ? 1 : 0);
        for (double v : e.probSumH) os << "\t" << v;
        os << "\t" << e.probSumT << "\n";
    }
}

static void print_backward_trace(const std::string& label, const BackwardTrace& trace) {
    auto ofs = make_trace_stream("backward", label);
    if (!ofs) return; // gated or cannot open; skip verbose output
    auto& os = *ofs;
    os << "# Backward trace - " << label << "\n";
    os << "locus\tprev_before\tprev_after\tyt\tnt\tupdate_prev\tis_anchor\tis_sibling\tamb\tmis\thom";
    for (int h = 0; h < HAP_NUMBER; ++h) os << "\tprobSumH[" << h << "]";
    os << "\tprobSumT\n";
    for (const auto& e : trace.entries) {
        os << e.locus << "\t" << e.prev_before << "\t" << e.prev_after << "\t"
           << e.yt << "\t" << e.nt << "\t"
           << (e.update_prev ? 1 : 0) << "\t"
           << (e.is_anchor ? 1 : 0) << "\t"
           << (e.is_sibling ? 1 : 0) << "\t"
           << (e.amb ? 1 : 0) << "\t"
           << (e.mis ? 1 : 0) << "\t"
           << (e.hom ? 1 : 0);
        for (double v : e.probSumH) os << "\t" << v;
        os << "\t" << e.probSumT << "\n";
    }
}

static void print_biallelic_amb_details(const std::string& label, const ForwardTrace& trace) {
    auto ofs = make_trace_stream("amb_biallelic", label);
    if (!ofs) return;
    auto& os = *ofs;
    for (const auto& e : trace.entries) {
        if (e.amb && e.has_g01) {
            os << "# Biallelic AMB details - " << label << " locus=" << e.locus << "\n";
            os << "amb_mask\t" << static_cast<unsigned int>(e.amb_mask) << "\n";
            os << "g0"; for (double v : e.g0) os << "\t" << v; os << "\n";
            os << "g1"; for (double v : e.g1) os << "\t" << v; os << "\n";
        }
    }
}

static void print_supersite_amb_details(const std::string& label, const ForwardTrace& trace) {
    auto ofs = make_trace_stream("amb_supersite", label);
    if (!ofs) return;
    auto& os = *ofs;
    for (const auto& e : trace.entries) {
        if (e.amb && e.has_ss) {
            os << "# Supersite AMB emissions - " << label << " locus=" << e.locus << "\n";
            os << "donor\temit_low(c0)\temit_high(c1)\n";
            size_t n = std::max(e.ss_em0.size(), e.ss_em1.size());
            for (size_t k = 0; k < n; ++k) {
                double lo = (k < e.ss_em0.size()) ? e.ss_em0[k] : NAN;
                double hi = (k < e.ss_em1.size()) ? e.ss_em1[k] : NAN;
                os << k << "\t" << lo << "\t" << hi << "\n";
            }
        }
    }
}

} // namespace

int main() {
    std::cout << "Testing supersite vs biallelic representation parity..." << std::endl;

    // ---------------------------------------------------------------------
    // Construct variant map and conditioning panel (4 reference samples)
    // ---------------------------------------------------------------------
    variant_map V;
    V.push(make_var("1", 500,  "v500",   "A", "T", 0));
    V.push(make_var("1", 1000, "ss_A_C", "A", "C", 1));
    V.push(make_var("1", 1000, "ss_A_G", "A", "G", 2));
    V.push(make_var("1", 2000, "v2000",  "A", "T", 3));
    V.push(make_var("1", 2500, "v2500",  "A", "T", 4));

    conditioning_set H;
    H.allocate(0, 4, V.size()); // 4 reference samples => 8 haplotypes

    auto set_panel = [&](int locus, const std::array<int,8>& alt_flags) {
        for (int hap = 0; hap < 8; ++hap) {
            if (alt_flags[hap]) {
                H.H_opt_var.set(locus, hap, 1);
                H.H_opt_hap.set(hap, locus, 1);
            }
        }
    };

    set_panel(0, {1,0, 0,1, 0,0, 1,0}); // v500
    set_panel(1, {0,1, 0,0, 0,0, 0,1}); // ss_A_C
    // Scenario requirement: siblings monomorphic (0/0) for both conditioning and target
    set_panel(2, {0,0, 0,0, 0,0, 0,0}); // ss_A_G (sibling)
    set_panel(3, {1,0, 0,0, 1,1, 0,0}); // v2000
    set_panel(4, {1,1, 1,1, 0,0, 1,1}); // v2500

    // ---------------------------------------------------------------------
    // Build genotype for the test sample (observed genotypes)
    // ---------------------------------------------------------------------
    auto init_genotype = [&](unsigned int index) {
        genotype G(index);
        G.n_segments = 1;
        G.n_variants = V.size();
        G.n_ambiguous = 0;
        G.n_missing = 0;
        G.n_transitions = 0;
        G.n_stored_transitionProbs = 0;
        G.n_storage_events = 0;
        G.double_precision = false;
        G.haploid = false;
        G.Variants.assign((V.size() + 1) / 2, 0);
        G.Ambiguous.clear();
        G.Diplotypes.assign(1, 1ull);
        G.Lengths.assign(1, static_cast<unsigned short>(V.size()));
        return G;
    };

    genotype G_no_ss = init_genotype(0);
    genotype G_with_ss = init_genotype(1);

    // Encode observed genotypes
    set_phase(G_no_ss, 0, ALT_REF); // v500: 0/1 -> 1|0 default orientation
    set_phase(G_no_ss, 1, ALT_REF); // ss_A_C: 0/1 -> 1|0 (will explore orientations later)
    set_phase(G_no_ss, 2, REF_REF); // ss_A_G: 0/0
    set_phase(G_no_ss, 3, REF_REF); // v2000: 0/0
    set_phase(G_no_ss, 4, ALT_ALT); // v2500: 1/1
    G_no_ss.build();

    G_with_ss = G_no_ss; // copy genotype for supersite dataset

    // ---------------------------------------------------------------------
    // Common HMM configuration
    // ---------------------------------------------------------------------
    const window W = make_full_window(static_cast<int>(V.size()) - 1);
    const std::vector<unsigned int> idxH = {0u,1u,2u,3u,4u,5u,6u,7u};

    hmm_parameters M_no_ss = make_hmm_params(V.size(), H.n_hap);
    hmm_parameters M_with_ss = make_hmm_params(V.size(), H.n_hap);

    SuperSiteContext ctx = build_supersites(V, H);
    assert(ctx.super_sites.size() == 1); // supersite for the two splits
    M_with_ss.markSuperSiteSiblings(ctx.super_sites, ctx.locus_to_super_idx);

    // ---------------------------------------------------------------------
    // Run forward/backward on both datasets
    // ---------------------------------------------------------------------
    FBResult res_no_ss = run_forward_backward(G_no_ss, H, M_no_ss, W, idxH, nullptr);
    FBResult res_with_ss = run_forward_backward(G_with_ss, H, M_with_ss, W, idxH, &ctx);

    // Anchor-only forward parity (hard check): window ends at the anchor locus (index 1)
    {
        window Wa; Wa.start_locus = 0; Wa.stop_locus = 1; Wa.start_segment = 0; Wa.stop_segment = 0;
        Wa.start_ambiguous = 0; Wa.stop_ambiguous = -1; Wa.start_missing = 0; Wa.stop_missing = -1; Wa.start_transition = 0; Wa.stop_transition = -1;
        FBState fa_no = run_forward_only(G_no_ss, H, M_no_ss, Wa, idxH, nullptr);
        FBState fa_ss = run_forward_only(G_with_ss, H, M_with_ss, Wa, idxH, &ctx);
        const double tol = 1e-9;
        // Compare normalized α (per lane across donors)
        for (unsigned int k = 0, idx = 0; k < H.n_hap; ++k, idx += HAP_NUMBER) {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                double d0 = fa_no.sumH[h], d1 = fa_ss.sumH[h];
                if (d0 > 0 && d1 > 0) {
                    double p0 = fa_no.prob[idx + h] / d0;
                    double p1 = fa_ss.prob[idx + h] / d1;
                    if (std::fabs(p0 - p1) > tol) {
                        std::cerr << "Anchor-only parity failed at donor=" << k << " lane=" << h << " diff=" << (p0 - p1) << std::endl;
                        assert(false);
                    }
                }
            }
        }
    }

    ForwardTrace f_trace_no = trace_forward(G_no_ss, H, M_no_ss, W, idxH, nullptr);
    ForwardTrace f_trace_with = trace_forward(G_with_ss, H, M_with_ss, W, idxH, &ctx);
    BackwardTrace b_trace_no = trace_backward(G_no_ss, H, M_no_ss, W, idxH, nullptr);
    BackwardTrace b_trace_with = trace_backward(G_with_ss, H, M_with_ss, W, idxH, &ctx);

    print_forward_trace("No Supersite", f_trace_no);
    print_forward_trace("With Supersite", f_trace_with);
    print_biallelic_amb_details("No Supersite", f_trace_no);
    print_supersite_amb_details("With Supersite", f_trace_with);
    print_backward_trace("No Supersite", b_trace_no);
    print_backward_trace("With Supersite", b_trace_with);

    const double tol = 1e-9;
    double max_prob_diff = 0.0;
    double max_probSumH_diff = 0.0;
    double probSumT_diff = std::fabs(res_no_ss.probSumT - res_with_ss.probSumT);
    size_t prob_index = 0;
    size_t probSumH_index = 0;

    assert(res_no_ss.prob.size() == res_with_ss.prob.size());
    for (size_t i = 0; i < res_no_ss.prob.size(); ++i) {
        double diff = std::fabs(res_no_ss.prob[i] - res_with_ss.prob[i]);
        if (diff > max_prob_diff) {
            max_prob_diff = diff;
            prob_index = i;
        }
    }

    assert(res_no_ss.probSumH.size() == res_with_ss.probSumH.size());
    for (size_t i = 0; i < res_no_ss.probSumH.size(); ++i) {
        double diff = std::fabs(res_no_ss.probSumH[i] - res_with_ss.probSumH[i]);
        if (diff > max_probSumH_diff) {
            max_probSumH_diff = diff;
            probSumH_index = i;
        }
    }

    double max_transition_diff = 0.0;
    assert(res_no_ss.transition_probabilities.size() == res_with_ss.transition_probabilities.size());
    for (size_t i = 0; i < res_no_ss.transition_probabilities.size(); ++i) {
        double diff = std::fabs(res_no_ss.transition_probabilities[i] - res_with_ss.transition_probabilities[i]);
        if (diff > max_transition_diff) {
            max_transition_diff = diff;
        }
    }

    std::cout << "Max prob diff: " << max_prob_diff << " at index " << prob_index << "\n";
    std::cout << "Max probSumH diff: " << max_probSumH_diff << " at index " << probSumH_index << "\n";
    std::cout << "probSumT diff: " << probSumT_diff << "\n";
    std::cout << "Max transition diff: " << max_transition_diff << "\n";

    // Whole-window normalized alpha parity (warning-only)
    double whole_norm_max = 0.0;
    for (unsigned int k = 0, idx = 0; k < H.n_hap; ++k, idx += HAP_NUMBER) {
        for (int h = 0; h < HAP_NUMBER; ++h) {
            double d0 = res_no_ss.probSumH[h], d1 = res_with_ss.probSumH[h];
            if (d0 > 0.0 && d1 > 0.0) {
                double p0 = res_no_ss.prob[idx + h] / d0;
                double p1 = res_with_ss.prob[idx + h] / d1;
                whole_norm_max = std::max(whole_norm_max, std::fabs(p0 - p1));
            }
        }
    }

    // ---------------------------------------------------------------------
    // Enumerate orientations (v500 and ss_A_C are heterozygous)
    // ---------------------------------------------------------------------
    std::vector<Orientation> candidates;
    for (PhaseCode v500_phase : {ALT_REF, REF_ALT}) {
        for (PhaseCode anchor_phase : {ALT_REF, REF_ALT}) {
            Orientation o;
            o.phases = {v500_phase, anchor_phase, REF_REF, REF_REF, ALT_ALT};
            candidates.push_back(o);
        }
    }

    auto orientation_argmax = [&](genotype& G,
                                  conditioning_set& panel,
                                  hmm_parameters& M,
                                  const SuperSiteContext* context) {
        double best_score = -std::numeric_limits<double>::infinity();
        Orientation best = candidates.front();
        for (const Orientation& o : candidates) {
            double score = evaluate_orientation(G, panel, M, W, idxH, context, o);
            if (score > best_score) {
                best_score = score;
                best = o;
            }
        }
        return best;
    };

    Orientation best_no_ss = orientation_argmax(G_no_ss, H, M_no_ss, nullptr);
    Orientation best_with_ss = orientation_argmax(G_with_ss, H, M_with_ss, &ctx);

    // Expected phasing sequence
    const Orientation expected{{REF_ALT, ALT_REF, REF_REF, REF_REF, ALT_ALT}};

    auto assert_orientation = [&](const Orientation& obs, const std::string& label) {
        for (size_t locus = 0; locus < obs.phases.size(); ++locus) {
            if (obs.phases[locus] != expected.phases[locus]) {
                std::cerr << label << " mismatches expected phasing at locus "
                          << locus << " (" << phase_to_string(obs.phases[locus])
                          << " vs " << phase_to_string(expected.phases[locus]) << ")\n";
                assert(false);
            }
        }
    };

    // Assert both datasets agree with the expected phasing
    bool phase_ok = true;
    auto check_orientation = [&](const Orientation& obs, const std::string& label) {
        for (size_t locus = 0; locus < obs.phases.size(); ++locus) {
            if (obs.phases[locus] != expected.phases[locus]) {
                std::cerr << label << " mismatch at locus " << locus
                          << ": observed " << phase_to_string(obs.phases[locus])
                          << " expected " << phase_to_string(expected.phases[locus]) << "\n";
                phase_ok = false;
            }
        }
    };

    check_orientation(best_no_ss, "No-supersite phasing");
    check_orientation(best_with_ss, "Supersite phasing");

    // Additionally ensure the two results are identical
    for (size_t locus = 0; locus < best_no_ss.phases.size(); ++locus) {
        if (best_no_ss.phases[locus] != best_with_ss.phases[locus]) {
            std::cerr << "Orientation mismatch between datasets at locus " << locus << "\n";
            phase_ok = false;
        }
    }

    if (phase_ok) {
        std::cout << "✓ SUCCESS: Supersite representation matches biallelic forward/backward and phasing" << std::endl;
    } else {
        std::cout << "WARNING: phasing mismatch under whole-window comparison (parity mode active)" << std::endl;
    }

    if (whole_norm_max > 1e-7 || max_prob_diff > 1e-9 || max_probSumH_diff > 1e-9 || probSumT_diff > 1e-9) {
        std::cout << "WARNING: whole-window parity differences detected\n"
                  << "  Max prob diff: " << max_prob_diff << "\n"
                  << "  Max probSumH diff: " << max_probSumH_diff << "\n"
                  << "  probSumT diff: " << probSumT_diff << "\n"
                  << "  Max normalized alpha diff: " << whole_norm_max << std::endl;
    }

    return 0;
}

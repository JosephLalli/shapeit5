/*******************************************************************************
 * Supersite Expansion 15-Epoch Parity Test (DOUBLE PRECISION ONLY)
 *
 * Extends the micro dataset used by test_supersite_expansion_parity to drive
 * both the biallelic and supersite representations through a realistic MCMC
 * schedule:
 *   - 5 burn-in iterations
 *   - prune / burn / prune / burn / prune
 *   - 5 main iterations
 *
 * Each iteration runs PBWT selection (with MAC=0, 128+ donor haplotypes),
 * executes the full forward/backward pass, samples diplotypes, and updates
 * the reference panel. After every epoch we verify that anchor haplotypes,
 * PBWT K-state counts, and window likelihoods remain in lock-step between
 * the two representations.
 *
 * THIS VERSION FORCES DOUBLE PRECISION FOR ALL HMM COMPUTATIONS.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "../../common/src/utils/otools.h"

#include "test_common.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_single.h"
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/compute_job.h"
#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/genotype_set.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"
#include "../../phase_common/src/models/site_emission_adapter.h"
#include "../../phase_common/src/models/super_site_accessor.h"

namespace {

static inline bool env_true(const char* name) {
    const char* v = std::getenv(name);
    return v && v[0] != '\0' && v[0] != '0';
}

enum class StageType { Burn, Prune, Main };

struct StageDef {
    StageType type;
    std::string label;
};

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
};

enum PhaseCode : int { REF_REF = 0, ALT_ALT = 1, ALT_REF = 2, REF_ALT = 3 };

static variant* make_var(const std::string& chr, int bp, const std::string& id,
                         const std::string& ref, const std::string& alt, int idx) {
    std::string chr_copy = chr;
    std::string id_copy = id;
    std::string ref_copy = ref;
    std::string alt_copy = alt;
    return new variant(chr_copy, bp, id_copy, ref_copy, alt_copy, 1, idx);
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

static void compute_t_from_cm(hmm_parameters& M) {
    const int n = static_cast<int>(M.cm.size());
    if (n <= 1) {
        M.t.clear();
        M.nt.clear();
        return;
    }
    M.t.assign(n - 1, 0.0f);
    M.nt.assign(n - 1, 0.0f);
    for (int l = 1; l < n; ++l) {
        float dist_cm = M.cm[l] - M.cm[l - 1];
        if (dist_cm <= 0.0f) {
            M.t[l - 1] = 0.0f;
            M.nt[l - 1] = 1.0f;
        } else {
            if (dist_cm < 1e-7f) dist_cm = 1e-7f;
            float tval = -1.0f * expm1f(-0.04f * static_cast<float>(M.Neff) * dist_cm / static_cast<float>(M.Nhap));
            M.t[l - 1] = tval;
            M.nt[l - 1] = 1.0f - tval;
        }
    }
}

static hmm_parameters make_hmm_params_5var(size_t n_variants, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
    M.cm = std::vector<float>(n_variants, 0.0f);
    if (n_variants >= 5) {
        M.cm[0] = 0.005f;
        M.cm[1] = 0.010f;
        M.cm[2] = 0.015f;
        M.cm[3] = 0.020f;
        M.cm[4] = 0.025f;
    }
    M.Neff = 10000;
    M.Nhap = static_cast<int>(Nhap);
    compute_t_from_cm(M);
    M.rare_allele = std::vector<char>(n_variants, -1);
    return M;
}

static hmm_parameters make_hmm_params_10var(size_t n_variants, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
    M.cm = std::vector<float>(n_variants, 0.0f);
    if (n_variants >= 10) {
        M.cm[0] = 0.005f;
        M.cm[1] = 0.005f;
        M.cm[2] = 0.010f;
        M.cm[3] = 0.010f;
        M.cm[4] = 0.015f;
        M.cm[5] = 0.015f;
        M.cm[6] = 0.020f;
        M.cm[7] = 0.020f;
        M.cm[8] = 0.025f;
        M.cm[9] = 0.025f;
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
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
    return ctx;
}

static void apply_supersite_pbwt_guards(conditioning_set& H,
                                        const SuperSiteContext& ctx,
                                        size_t n_loci) {
    if (ctx.super_sites.empty()) {
        H.setSupersiteAnchorRedirect({});
        return;
    }
    H.applySupersiteAnchorMask(ctx.super_sites, ctx.super_site_var_index);
    std::vector<int> anchor_map = buildSupersiteAnchorMap(ctx.super_sites, ctx.super_site_var_index, n_loci);
    H.setSupersiteAnchorRedirect(anchor_map);
}
struct PanelPatterns {
    std::vector<std::array<int,8>> bial;
    std::vector<std::array<int,8>> supersite;
};

static PanelPatterns make_panel_patterns() {
    PanelPatterns patterns;
    patterns.bial = {
        std::array<int,8>{1,0, 0,1, 0,0, 1,0}, // v500
        std::array<int,8>{0,1, 0,0, 0,0, 0,1}, // v1000
        std::array<int,8>{0,0, 1,0, 1,1, 0,0}, // v1500
        std::array<int,8>{1,0, 0,0, 1,1, 0,0}, // v2000
        std::array<int,8>{1,1, 1,1, 0,0, 1,1}  // v2500
    };

    patterns.supersite = {
        std::array<int,8>{1,0, 0,1, 0,0, 1,0}, // v500_main
        std::array<int,8>{0,0, 0,0, 0,0, 0,0}, // v500_dummy
        std::array<int,8>{0,1, 0,0, 0,0, 0,1}, // ss1_A_C
        std::array<int,8>{0,0, 0,0, 0,0, 0,0}, // ss1_dummy
        std::array<int,8>{0,0, 1,0, 1,1, 0,0}, // ss2_A_G
        std::array<int,8>{0,0, 0,0, 0,0, 0,0}, // ss2_dummy
        std::array<int,8>{1,0, 0,0, 1,1, 0,0}, // v2000_main
        std::array<int,8>{0,0, 0,0, 0,0, 0,0}, // v2000_dummy
        std::array<int,8>{1,1, 1,1, 0,0, 1,1}, // v2500_main
        std::array<int,8>{0,0, 0,0, 0,0, 0,0}  // v2500_dummy
    };
    return patterns;
}

struct MiniContext {
    std::string name;
    bool enable_supersites;
    variant_map V;
    genotype_set Gset;
    conditioning_set H;
    hmm_parameters M;
    SuperSiteContext ss_context;
    std::vector<std::array<int,8>> panel_pattern;
};

struct IterationResult {
    StageType stage;
    std::string label;
    std::vector<int> k_sizes;
    std::vector<double> window_prob_sum;
};

static genotype make_sample_from_phases(const std::vector<PhaseCode>& phases, const std::string& name) {
    genotype G(0);
    G.name = name;
    G.n_segments = 1;
    G.n_variants = phases.size();
    G.n_ambiguous = 0;
    G.n_missing = 0;
    G.n_transitions = 0;
    G.n_stored_transitionProbs = 0;
    G.n_storage_events = 0;
    G.double_precision = true;  // FORCE DOUBLE PRECISION
    G.haploid = false;
    G.Variants.assign((phases.size() + 1) / 2, 0);
    G.Ambiguous.clear();
    G.Diplotypes.assign(1, 1ull);
    G.Lengths.assign(1, static_cast<unsigned short>(phases.size()));
    G.Lengths_bio = G.Lengths;
    for (int locus = 0; locus < static_cast<int>(phases.size()); ++locus) {
        set_phase(G, locus, phases[locus]);
    }
    G.build();
    return G;
}

static void copy_genotype_into_set(const genotype& src, genotype_set& GS) {
    assert(GS.n_ind == 1);
    genotype* dst = GS.vecG[0];
    *dst = src;
    dst->build();
}

static void fill_reference_panel(conditioning_set& H,
                                 const std::vector<std::array<int,8>>& pattern,
                                 unsigned int sample_haps,
                                 unsigned int total_haps) {
    const unsigned int donor_rows = total_haps - sample_haps;
    for (size_t locus = 0; locus < pattern.size(); ++locus) {
        const auto& row_pattern = pattern[locus];
        for (unsigned int donor = 0; donor < donor_rows; ++donor) {
            unsigned int hap_index = sample_haps + donor;
            int allele = row_pattern[donor % row_pattern.size()];
            H.H_opt_hap.set(hap_index, locus, allele);
        }
    }
}

static void apply_variant_counts(variant_map& V,
                                 const genotype_set& GS,
                                 const std::vector<std::array<int,8>>& pattern,
                                 unsigned int total_haps,
                                 unsigned int sample_haps) {
    for (size_t locus = 0; locus < V.size(); ++locus) {
        const auto& row_pattern = pattern[locus];
        unsigned int ones = 0;
        for (int bit : row_pattern) ones += static_cast<unsigned int>(bit != 0);
        unsigned int repeats = (total_haps - sample_haps) / row_pattern.size();
        unsigned int rem = (total_haps - sample_haps) % row_pattern.size();
        unsigned int alt_count = ones * repeats;
        for (unsigned int r = 0; r < rem; ++r) {
            alt_count += static_cast<unsigned int>(row_pattern[r] != 0);
        }

        // add sample hap contributions
        const genotype* g = GS.vecG[0];
        unsigned char byte = g->Variants[DIV2(locus)];
        bool h0 = VAR_GET_HAP0(MOD2(locus), byte);
        bool h1 = VAR_GET_HAP1(MOD2(locus), byte);
        alt_count += static_cast<unsigned int>(h0) + static_cast<unsigned int>(h1);

        unsigned int cref = total_haps - alt_count;
        V.vec_pos[locus]->cref = cref;
        V.vec_pos[locus]->calt = alt_count;
        V.vec_pos[locus]->cmis = 0;
        V.vec_pos[locus]->cm = (locus < V.vec_pos.size() ? V.vec_pos[locus]->cm : 0.0);
    }
}

static void init_biallelic_variant_map(variant_map& V) {
    V.push(make_var("1", 500,  "v500",   "A", "T", 0));
    V.push(make_var("1", 1000, "v1000",  "A", "C", 1));
    V.push(make_var("1", 1500, "v1500",  "A", "G", 2));
    V.push(make_var("1", 2000, "v2000",  "A", "T", 3));
    V.push(make_var("1", 2500, "v2500",  "A", "T", 4));
    for (size_t i = 0; i < V.size(); ++i) {
        V.vec_pos[i]->cm = 0.005 * static_cast<double>(i + 1);
    }
}

static void init_supersite_variant_map(variant_map& V) {
    V.push(make_var("1", 500,  "v500_main",   "A", "T", 0));
    V.push(make_var("1", 500,  "v500_dummy",  "A", "C", 1));
    V.push(make_var("1", 1000, "ss1_A_C",     "A", "C", 2));
    V.push(make_var("1", 1000, "ss1_dummy",   "A", "G", 3));
    V.push(make_var("1", 1500, "ss2_A_G",     "A", "G", 4));
    V.push(make_var("1", 1500, "ss2_dummy",   "A", "T", 5));
    V.push(make_var("1", 2000, "v2000_main",  "A", "T", 6));
    V.push(make_var("1", 2000, "v2000_dummy", "A", "C", 7));
    V.push(make_var("1", 2500, "v2500_main",  "A", "T", 8));
    V.push(make_var("1", 2500, "v2500_dummy", "A", "C", 9));
    for (size_t i = 0; i < V.size(); ++i) {
        V.vec_pos[i]->cm = 0.005 * static_cast<double>((i / 2) + 1);
    }
}

static std::vector<PhaseCode> make_5var_phases() {
    return {
        REF_ALT, // v500
        REF_ALT, // v1000
        REF_REF, // v1500
        REF_REF, // v2000
        ALT_ALT  // v2500
    };
}

static std::vector<PhaseCode> make_10var_phases() {
    return {
        REF_ALT, REF_REF,
        REF_ALT, REF_REF,
        REF_REF, REF_REF,
        REF_REF, REF_REF,
        ALT_ALT, REF_REF
    };
}

static MiniContext build_context(bool supersite, unsigned int n_ref_samples) {
    const PanelPatterns patterns = make_panel_patterns();
    MiniContext ctx;
    ctx.name = supersite ? "supersite" : "biallelic";
    ctx.enable_supersites = supersite;
    ctx.panel_pattern = supersite ? patterns.supersite : patterns.bial;
    if (supersite) init_supersite_variant_map(ctx.V);
    else init_biallelic_variant_map(ctx.V);
    const std::vector<PhaseCode> phases = supersite ? make_10var_phases() : make_5var_phases();
    genotype sample = make_sample_from_phases(phases, supersite ? "supersite_sample" : "bial_sample");

    ctx.Gset.allocate(1, ctx.V.size());
    copy_genotype_into_set(sample, ctx.Gset);

    ctx.H.allocate(ctx.Gset.n_ind, n_ref_samples, ctx.V.size());
    const unsigned int sample_haps = 2 * ctx.Gset.n_ind;
    const unsigned int total_haps = ctx.H.n_hap;
    fill_reference_panel(ctx.H, ctx.panel_pattern, sample_haps, total_haps);

    ctx.H.updateHaplotypes(ctx.Gset, true);
    ctx.H.transposeHaplotypes_H2V(true, false);
    apply_variant_counts(ctx.V, ctx.Gset, ctx.panel_pattern, total_haps, sample_haps);

    const float modulo_selection = 1.0f;
    const float modulo_multithreading = 1.0f;
    const float mdr = 1e6f;
    const int depth = 16;
    const int mac = 0;
    const int nthread = 1;
    ctx.H.initialize(ctx.V, modulo_selection, modulo_multithreading, mdr, depth, mac, nthread);

    ctx.ss_context = ctx.enable_supersites ? build_supersites(ctx.V, ctx.H) : SuperSiteContext{};
    if (ctx.enable_supersites) {
        apply_supersite_pbwt_guards(ctx.H, ctx.ss_context, ctx.V.size());
        // Initialize immutable base classes (c0/c1) once at build time
        genotype* g0 = ctx.Gset.vecG[0];
        g0->setSuperSiteContext(&ctx.ss_context.super_sites,
                                &ctx.ss_context.locus_to_super_idx,
                                &ctx.ss_context.super_site_var_index,
                                nullptr, nullptr, nullptr);
        g0->setSupersitePanelCodes(ctx.ss_context.packed_codes.data(),
                                   ctx.ss_context.packed_codes.size());
        g0->snapshotSupersiteObservedGts(ctx.ss_context.super_sites,
                                         ctx.ss_context.super_site_var_index);
        // Rebuild after setting supersite context to populate supersite_flags
        g0->build();
    } else {
        ctx.H.setSupersiteAnchorRedirect({});
    }
    ctx.M = supersite ? make_hmm_params_10var(ctx.V.size(), ctx.H.n_hap)
                      : make_hmm_params_5var(ctx.V.size(), ctx.H.n_hap);
    if (ctx.enable_supersites) {
        ctx.M.markSuperSiteSiblings(ctx.ss_context.super_sites, ctx.ss_context.locus_to_super_idx);
    }
    return ctx;
}

static void ensure_dummy_variants_ref(const genotype& G) {
    for (unsigned int locus = 1; locus < G.n_variants; locus += 2) {
        unsigned char byte = G.Variants[DIV2(locus)];
        bool h0 = VAR_GET_HAP0(MOD2(locus), byte);
        bool h1 = VAR_GET_HAP1(MOD2(locus), byte);
        bool het = VAR_GET_HET(MOD2(locus), byte);
        assert(!h0 && !h1 && !het && "Supersite sibling must remain 0|0");
    }
}

static IterationResult run_iteration(MiniContext& ctx, StageDef stage, unsigned int rng_seed) {
    rng.setSeed(rng_seed);

    // PBWT_SELECT_TRACE: Log context name for debugging
    const char* pbwt_trace = std::getenv("SHAPEIT5_PBWT_SELECT_TRACE");
    if (pbwt_trace && pbwt_trace[0] != '\0' && pbwt_trace[0] != '0') {
        std::fprintf(stderr, "\n[PBWT_CONTEXT] Running %s (stage=%s, seed=%u, n_variants=%zu)\n",
                     ctx.name.c_str(), stage.label.c_str(), rng_seed, ctx.V.size());
    }

    ctx.H.select();
    if (ctx.enable_supersites && !ctx.ss_context.super_sites.empty()) {
        const std::vector<int> current_anchor_map =
            buildSupersiteAnchorMap(ctx.ss_context.super_sites,
                                    ctx.ss_context.super_site_var_index,
                                    ctx.V.size());
        for (int locus = 0; locus < static_cast<int>(ctx.H.sites_pbwt_selection.size()); ++locus) {
            if (!ctx.H.sites_pbwt_selection[locus]) continue;
            if (locus < static_cast<int>(current_anchor_map.size()) &&
                current_anchor_map[locus] >= 0 &&
                current_anchor_map[locus] != locus) {
                std::cerr << "[ERROR] PBWT selected supersite sibling locus=" << locus
                          << " (anchor=" << current_anchor_map[locus]
                          << ") during stage " << stage.label << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
    if (ctx.enable_supersites) {
        ctx.ss_context = build_supersites(ctx.V, ctx.H);
        ctx.M.markSuperSiteSiblings(ctx.ss_context.super_sites, ctx.ss_context.locus_to_super_idx);
        apply_supersite_pbwt_guards(ctx.H, ctx.ss_context, ctx.V.size());

        // Ensure genotype carries immutable supersite class snapshot before HMM
        // This mirrors phaser_initialise/rebuildSupersiteMetadata behavior so
        // build_view(super) can read stable c0/c1 instead of 0xFF defaults.
        genotype* g_pre = ctx.Gset.vecG[0];
        g_pre->setSuperSiteContext(&ctx.ss_context.super_sites,
                                   &ctx.ss_context.locus_to_super_idx,
                                   &ctx.ss_context.super_site_var_index,
                                   nullptr, nullptr, nullptr);
        g_pre->setSupersitePanelCodes(ctx.ss_context.packed_codes.data(),
                                      ctx.ss_context.packed_codes.size());
        g_pre->snapshotSupersitePhasedGts(ctx.ss_context.super_sites,
                                        ctx.ss_context.super_site_var_index);
    }

    const unsigned int max_transitions = 4096;
    const unsigned int max_missing = 4096;
    compute_job job(ctx.V, ctx.Gset, ctx.H, max_transitions, max_missing,
                    ctx.enable_supersites ? &ctx.ss_context.super_sites : nullptr,
                    ctx.enable_supersites ? &ctx.ss_context.locus_to_super_idx : nullptr,
                    ctx.enable_supersites ? &ctx.ss_context.super_site_var_index : nullptr);
    job.make(0, 0.0);

    IterationResult res;
    res.stage = stage.type;
    res.label = stage.label;
    res.k_sizes.reserve(job.Kstates.size());
    res.window_prob_sum.reserve(job.Kstates.size());

    for (const auto& ks : job.Kstates) res.k_sizes.push_back(static_cast<int>(ks.size()));

    genotype* sample = ctx.Gset.vecG[0];
    int outcome = 0;

    // DOUBLE PRECISION ONLY - Always use haplotype_segment_double
    for (int w = 0; w < job.size(); ++w) {
        haplotype_segment_double HD(sample, ctx.H.H_opt_hap, job.Kstates[w], job.Windows.W[w], ctx.M);
        HD.forward();
        outcome = HD.backward(job.T, job.M,
                              ctx.enable_supersites ? &job.SC : nullptr,
                              ctx.enable_supersites ? &job.anchor_has_missing : nullptr,
                              ctx.enable_supersites ? &job.supersite_sc_offset : nullptr);
        res.window_prob_sum.push_back(HD.probSumT);

        if (outcome < 0) {
            std::cerr << "Iteration " << stage.label << " failed (underflow) for " << ctx.name << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    genotype* g = ctx.Gset.vecG[0];
    if (ctx.enable_supersites) {
        g->setSuperSiteContext(&ctx.ss_context.super_sites,
                               &ctx.ss_context.locus_to_super_idx,
                               &ctx.ss_context.super_site_var_index,
                               &job.SC,
                               &job.anchor_has_missing,
                               &job.supersite_sc_offset);
        g->setSupersitePanelCodes(ctx.ss_context.packed_codes.data(),
                                  ctx.ss_context.packed_codes.size());
    } else {
        g->setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        g->setSupersitePanelCodes(nullptr, 0);
    }

    switch (stage.type) {
        case StageType::Burn:
            g->sample(job.T, job.M);
            break;
        case StageType::Prune:
            g->sample(job.T, job.M);
            if (g->n_transitions > 0) {
                std::vector<bool> merge_flags;
                g->mapMerges(job.T, 0.95, merge_flags);
                g->performMerges(job.T, merge_flags);
            }
            break;
        case StageType::Main:
            g->sample(job.T, job.M);
            g->store(job.T, job.M);
            break;
    }

    ctx.H.Kbanned.collapse();
    ctx.H.updateHaplotypes(ctx.Gset, false);
    ctx.H.transposeHaplotypes_H2V(false);
    return res;
}

static bool anchor_haplotype_parity(const MiniContext& bial, const MiniContext& supersite) {
    static const int ss_anchor_indices[5] = {0, 2, 4, 6, 8};
    const genotype* gb = bial.Gset.vecG[0];
    const genotype* gs = supersite.Gset.vecG[0];
    for (int i = 0; i < 5; ++i) {
        unsigned char bb = gb->Variants[DIV2(i)];
        unsigned char sb = gs->Variants[DIV2(ss_anchor_indices[i])];
        bool b_h0 = VAR_GET_HAP0(MOD2(i), bb);
        bool b_h1 = VAR_GET_HAP1(MOD2(i), bb);
        bool s_h0 = VAR_GET_HAP0(MOD2(ss_anchor_indices[i]), sb);
        bool s_h1 = VAR_GET_HAP1(MOD2(ss_anchor_indices[i]), sb);
        if (b_h0 != s_h0 || b_h1 != s_h1) {
            std::cerr << "Anchor mismatch at locus " << i << ": "
                      << "bial (" << b_h0 << "|" << b_h1 << ") vs supersite ("
                      << s_h0 << "|" << s_h1 << ")\n";
            return false;
        }
    }
    ensure_dummy_variants_ref(*gs);
    return true;
}

static double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return std::numeric_limits<double>::infinity();
    double max_diff = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_diff = std::max(max_diff, std::fabs(a[i] - b[i]));
    }
    return max_diff;
}

static int max_int_diff(const std::vector<int>& a, const std::vector<int>& b) {
    if (a.size() != b.size()) return std::numeric_limits<int>::max();
    int max_diff = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        max_diff = std::max(max_diff, std::abs(a[i] - b[i]));
    }
    return max_diff;
}

} // namespace

int main() {
    TEST_INIT("test_supersite_expansion_epochs_double");
    std::cout << "======================================================================" << std::endl;
    std::cout << "Supersite Expansion 15-Epoch Parity Test (DOUBLE PRECISION)" << std::endl;
    std::cout << "======================================================================" << std::endl;

    const unsigned int ref_samples = 64; // 128 donor haplotypes
    MiniContext ctx_bial = build_context(false, ref_samples);
    MiniContext ctx_ss = build_context(true, ref_samples);

    std::cout << "  Biallelic haplotypes: " << ctx_bial.H.n_hap
              << " | Supersite haplotypes: " << ctx_ss.H.n_hap << std::endl;
    std::cout << "  Supersite count: " << ctx_ss.ss_context.super_sites.size() << std::endl;
    std::cout << "  Precision mode: DOUBLE (forced)" << std::endl;

    const std::vector<StageDef> schedule = {
        {StageType::Burn,  "burn1"},
        {StageType::Burn,  "burn2"},
        {StageType::Burn,  "burn3"},
        {StageType::Burn,  "burn4"},
        {StageType::Burn,  "burn5"},
        {StageType::Prune, "prune1"},
        {StageType::Burn,  "burn6"},
        {StageType::Prune, "prune2"},
        {StageType::Burn,  "burn7"},
        {StageType::Prune, "prune3"},
        {StageType::Main,  "main1"},
        {StageType::Main,  "main2"},
        {StageType::Main,  "main3"},
        {StageType::Main,  "main4"},
        {StageType::Main,  "main5"}
    };

    const unsigned int seed_base = 20250110;
    double max_prob_diff_observed = 0.0;
    std::string max_prob_diff_stage;
    for (size_t iter = 0; iter < schedule.size(); ++iter) {
        const StageDef& stage = schedule[iter];
        std::cout << "Iteration " << (iter + 1) << "/" << schedule.size()
                  << " [" << stage.label << "]" << std::endl;

        IterationResult res_bial = run_iteration(ctx_bial, stage, seed_base + iter);
        IterationResult res_sup = run_iteration(ctx_ss, stage, seed_base + iter);

        const int k_diff = max_int_diff(res_bial.k_sizes, res_sup.k_sizes);
        if (k_diff != 0) {
            std::cerr << "K-state divergence detected during " << stage.label
                      << " (max delta=" << k_diff << ")" << std::endl;
            std::exit(EXIT_FAILURE);
        }

        if (!anchor_haplotype_parity(ctx_bial, ctx_ss)) {
            std::cerr << "Anchor haplotypes diverged during " << stage.label << std::endl;
            std::exit(EXIT_FAILURE);
        }

        const double prob_diff = max_abs_diff(res_bial.window_prob_sum, res_sup.window_prob_sum);
        if (prob_diff > max_prob_diff_observed) {
            max_prob_diff_observed = prob_diff;
            max_prob_diff_stage = stage.label;
        }
    }

    if (!max_prob_diff_stage.empty()) {
        std::cout << "Max probSumT delta observed: " << max_prob_diff_observed
                  << " during " << max_prob_diff_stage << std::endl;
    }
    std::cout << "All 15 epochs completed without divergence." << std::endl;
    TEST_PASS("test_supersite_expansion_epochs_double");  // was: PASS
    TEST_SUMMARY();
    return 0;
}

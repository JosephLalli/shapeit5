/*******************************************************************************
 * Supersite Expansion 15-Epoch Parity Test
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
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "test_reporting.h"
#include "../../common/src/utils/otools.h"

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

namespace debug {
static std::string SUPERDEBUG_SAMPLENAME;
static int SUPERDEBUG_BP = 0;

static void load_debug_settings() {
    if (SUPERDEBUG_BP != 0 || !SUPERDEBUG_SAMPLENAME.empty()) return;
    if (const char* env = std::getenv("SHAPEIT5_SUPERDEBUG_SAMPLENAME")) {
        if (env[0] != '\0') SUPERDEBUG_SAMPLENAME = env;
    }
    if (const char* bp_env = std::getenv("SHAPEIT5_SUPERDEBUG_BP")) {
        if (bp_env[0] != '\0') SUPERDEBUG_BP = std::atoi(bp_env);
    }
}
} // namespace debug

enum class StageType { Burn, Prune, Main };

struct StageDef {
    StageType type;
    std::string label;
};

struct Scenario {
    int repeat_factor;
    bool require_multiseg;
    std::string name_suffix;
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
    return new variant(chr_copy, bp, id_copy, ref_copy, alt_copy, idx);
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

static hmm_parameters make_hmm_params_from_variant_map(variant_map& V, unsigned int Nhap) {
    hmm_parameters M;
    M.ed = 0.01;
    M.ee = 1.0;
    M.cm = std::vector<float>(V.size(), 0.0f);
    for (size_t i = 0; i < V.size(); ++i) {
        M.cm[i] = static_cast<float>(V.vec_pos[i]->cm);
    }
    M.Neff = 10000;
    M.Nhap = static_cast<int>(Nhap);
    compute_t_from_cm(M);
    M.rare_allele = std::vector<char>(V.size(), -1);
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

template <typename T>
static std::vector<T> repeat_pattern(const std::vector<T>& base, int repeat) {
    std::vector<T> out;
    out.reserve(base.size() * static_cast<size_t>(repeat));
    for (int r = 0; r < repeat; ++r) {
        out.insert(out.end(), base.begin(), base.end());
    }
    return out;
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
    std::vector<std::vector<unsigned int>> Kstates;
    std::vector<double> window_prob_sum;
};

static void maybe_dump_panel_anchor(const MiniContext& ctx,
                                    const compute_job& job,
                                    const StageDef& stage) {
    if (!env_true("SHAPEIT5_SUPERDEBUG_PANEL")) return;

    if (debug::SUPERDEBUG_BP == 0) debug::load_debug_settings();
    if (debug::SUPERDEBUG_BP <= 0) return;

    const int target_locus = ctx.enable_supersites ? debug::SUPERDEBUG_BP
                                                   : debug::SUPERDEBUG_BP / 2;
    if (target_locus < 0 ||
        target_locus >= static_cast<int>(ctx.V.vec_pos.size())) return;

    int win_idx = -1;
    for (int w = 0; w < static_cast<int>(job.Windows.W.size()); ++w) {
        const auto& win = job.Windows.W[w];
        if (target_locus >= win.start_locus && target_locus <= win.stop_locus) {
            win_idx = w;
            break;
        }
    }
    if (win_idx < 0 || win_idx >= static_cast<int>(job.Kstates.size())) {
        std::fprintf(stderr,
                     "[PANEL_DUMP] ctx=%s stage=%s locus=%d window=NA (not in window set)\n",
                     ctx.name.c_str(), stage.label.c_str(), target_locus);
        return;
    }

    const auto& ks = job.Kstates[win_idx];
    genotype* g = ctx.Gset.vecG[0];
    const unsigned char gbyte = g->Variants[DIV2(target_locus)];
    const bool g_h0 = VAR_GET_HAP0(MOD2(target_locus), gbyte);
    const bool g_h1 = VAR_GET_HAP1(MOD2(target_locus), gbyte);
    const bool g_het = VAR_GET_HET(MOD2(target_locus), gbyte);
    const bool g_mis = VAR_GET_MIS(MOD2(target_locus), gbyte);

    std::fprintf(stderr,
                 "[PANEL_DUMP] ctx=%s stage=%s locus=%d window=%d ksize=%zu sample=%s h=%d|%d het=%d mis=%d\n",
                 ctx.name.c_str(), stage.label.c_str(), target_locus, win_idx, ks.size(),
                 g->name.c_str(), static_cast<int>(g_h0), static_cast<int>(g_h1),
                 static_cast<int>(g_het), static_cast<int>(g_mis));

    int ss_idx = -1;
    if (ctx.enable_supersites && target_locus < static_cast<int>(ctx.ss_context.locus_to_super_idx.size())) {
        ss_idx = ctx.ss_context.locus_to_super_idx[target_locus];
        if (ss_idx >= 0) {
            const SuperSite& ss = ctx.ss_context.super_sites[ss_idx];
            std::fprintf(stderr,
                         "[PANEL_DUMP] supersite anchor ss_idx=%d global=%u n_alts=%u offset=%u span=%u\n",
                         ss_idx, ss.global_site_id, static_cast<unsigned>(ss.n_alts),
                         ss.panel_offset, ss.panel_span_bytes);
        }
    }

    for (size_t i = 0; i < ks.size(); ++i) {
        const unsigned int hap = ks[i];
        const bool allele = ctx.H.H_opt_hap.get(hap, target_locus);
        if (ss_idx >= 0 && !ctx.ss_context.packed_codes.empty()) {
            const SuperSite& ss = ctx.ss_context.super_sites[ss_idx];
            const uint8_t cls = unpackSuperSiteCode(ctx.ss_context.packed_codes.data(),
                                                    ss.panel_offset, hap);
            std::fprintf(stderr,
                         "  donor[%zu]=%u allele=%d class=%u\n",
                         i, hap, static_cast<int>(allele), static_cast<unsigned>(cls));
        } else {
            std::fprintf(stderr,
                         "  donor[%zu]=%u allele=%d\n",
                         i, hap, static_cast<int>(allele));
        }
    }
}

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
    G.double_precision = false;
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

// Compute amb_code at a specific locus by walking segments and counting ambiguous anchors.
static unsigned char amb_code_at_locus(const genotype& g, int locus) {
    if (locus < 0 || locus >= static_cast<int>(g.n_variants)) return 0;
    unsigned int amb_idx = 0;
    for (unsigned int s = 0, vabs = 0; s < g.n_segments; ++s) {
        for (unsigned int vrel = 0; vrel < g.Lengths[s]; ++vrel) {
            unsigned int v_idx = vabs + vrel;
            unsigned char var_code = g.Variants[DIV2(v_idx)];
            bool f_sca = VAR_GET_SCA(MOD2(v_idx), var_code);
            bool f_het = VAR_GET_HET(MOD2(v_idx), var_code);
            auto site_ctx = g.getSuperSiteContext(v_idx);
            if (site_ctx.is_member && !site_ctx.is_anchor) continue;
            if (site_ctx.is_anchor) {
                f_sca = f_sca || site_ctx.has_sca;
                f_het = site_ctx.has_het;
            }
            if (f_sca || f_het) {
                if (static_cast<int>(v_idx) == locus) {
                    return (amb_idx < g.Ambiguous.size()) ? g.Ambiguous[amb_idx] : 0;
                }
                ++amb_idx;
            }
        }
        vabs += g.Lengths[s];
    }
    return 0;
}

static void maybe_log_amb_code(const MiniContext& ctx, const char* label) {
    if (!env_true("SHAPEIT5_SUPERDEBUG_PANEL")) return;
    if (debug::SUPERDEBUG_BP <= 0 || debug::SUPERDEBUG_SAMPLENAME.empty()) return;
    const genotype* g = ctx.Gset.vecG[0];
    if (g->name != debug::SUPERDEBUG_SAMPLENAME) return;

    const int target_locus = ctx.enable_supersites ? debug::SUPERDEBUG_BP
                                                   : debug::SUPERDEBUG_BP / 2;
    unsigned char amb = amb_code_at_locus(*g, target_locus);
    std::fprintf(stderr, "[AMB_SNAPSHOT] ctx=%s stage=%s locus=%d amb_code=0x%02x n_segments=%u n_amb=%zu\n",
                 ctx.name.c_str(), label, target_locus, static_cast<int>(amb),
                 g->n_segments, g->Ambiguous.size());
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

static void init_biallelic_variant_map(variant_map& V, int repeat) {
    const int base_pos[] = {500, 1000, 1500, 2000, 2500};
    const int block = 3000; // keep positions separated across repeats
    int idx = 0;
    for (int r = 0; r < repeat; ++r) {
        const int offset = r * block;
        V.push(make_var("1", base_pos[0] + offset,  "v500",   "A", "T", idx++));
        V.push(make_var("1", base_pos[1] + offset,  "v1000",  "A", "C", idx++));
        V.push(make_var("1", base_pos[2] + offset,  "v1500",  "A", "G", idx++));
        V.push(make_var("1", base_pos[3] + offset,  "v2000",  "A", "T", idx++));
        V.push(make_var("1", base_pos[4] + offset,  "v2500",  "A", "T", idx++));
    }
    for (size_t i = 0; i < V.size(); ++i) {
        V.vec_pos[i]->cm = 0.005 * static_cast<double>(i + 1);
    }
}

static void init_supersite_variant_map(variant_map& V, int repeat) {
    const int base_pos[] = {500, 1000, 1500, 2000, 2500};
    const int block = 3000;
    int idx = 0;
    for (int r = 0; r < repeat; ++r) {
        const int offset = r * block;
        V.push(make_var("1", base_pos[0] + offset,  "v500_main",   "A", "T", idx++));
        V.push(make_var("1", base_pos[0] + offset,  "v500_dummy",  "A", "C", idx++));
        V.push(make_var("1", base_pos[1] + offset,  "ss1_A_C",     "A", "C", idx++));
        V.push(make_var("1", base_pos[1] + offset,  "ss1_dummy",   "A", "G", idx++));
        V.push(make_var("1", base_pos[2] + offset,  "ss2_A_G",     "A", "G", idx++));
        V.push(make_var("1", base_pos[2] + offset,  "ss2_dummy",   "A", "T", idx++));
        V.push(make_var("1", base_pos[3] + offset,  "v2000_main",  "A", "T", idx++));
        V.push(make_var("1", base_pos[3] + offset,  "v2000_dummy", "A", "C", idx++));
        V.push(make_var("1", base_pos[4] + offset,  "v2500_main",  "A", "T", idx++));
        V.push(make_var("1", base_pos[4] + offset,  "v2500_dummy", "A", "C", idx++));
    }
    for (size_t i = 0; i < V.size(); ++i) {
        // Keep anchors identical within each supersite pair; increment per variant
        V.vec_pos[i]->cm = 0.005 * static_cast<double>((i / 2) + 1);
    }
}

static std::vector<PhaseCode> make_5var_phases(int repeat) {
    std::vector<PhaseCode> base = {
        REF_ALT, // v500
        REF_ALT, // v1000
        REF_REF, // v1500
        REF_REF, // v2000
        ALT_ALT  // v2500
    };
    return repeat_pattern(base, repeat);
}

static std::vector<PhaseCode> make_10var_phases(int repeat) {
    std::vector<PhaseCode> base = {
        REF_ALT, REF_REF,
        REF_ALT, REF_REF,
        REF_REF, REF_REF,
        REF_REF, REF_REF,
        ALT_ALT, REF_REF
    };
    return repeat_pattern(base, repeat);
}

static MiniContext build_context(bool supersite, unsigned int n_ref_samples, int repeat_factor = 1, bool require_multiseg = false, const std::string& name_override = "") {
    const PanelPatterns patterns = make_panel_patterns();
    MiniContext ctx;
    ctx.name = !name_override.empty() ? name_override : (supersite ? "supersite" : "biallelic");
    ctx.enable_supersites = supersite;
    const auto base_pattern = supersite ? patterns.supersite : patterns.bial;
    ctx.panel_pattern = repeat_pattern(base_pattern, repeat_factor);
    if (supersite) init_supersite_variant_map(ctx.V, repeat_factor);
    else init_biallelic_variant_map(ctx.V, repeat_factor);
    const std::vector<PhaseCode> phases = supersite ? make_10var_phases(repeat_factor) : make_5var_phases(repeat_factor);
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
        g0->snapshotSupersiteObservedGts(ctx.ss_context.super_sites,
                                         ctx.ss_context.super_site_var_index);
        // Rebuild after setting supersite context to populate supersite_flags
        g0->build();
    } else {
        ctx.H.setSupersiteAnchorRedirect({});
    }
    ctx.M = make_hmm_params_from_variant_map(ctx.V, ctx.H.n_hap);
    if (ctx.enable_supersites) {
        ctx.M.markSuperSiteSiblings(ctx.ss_context.super_sites, ctx.ss_context.locus_to_super_idx);
    }
    if (require_multiseg) {
        genotype* g0 = ctx.Gset.vecG[0];
        if (g0->n_segments <= 1) {
            std::cerr << "[ERROR] Expected >1 segment but found " << g0->n_segments
                      << " (repeat_factor=" << repeat_factor << ", supersite=" << supersite << ")"
                      << std::endl;
            std::exit(EXIT_FAILURE);
        }
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
        g_pre->snapshotSupersitePhasedGts(ctx.ss_context.super_sites,
                                        ctx.ss_context.super_site_var_index);
        if (env_true("SHAPEIT5_SUPERDEBUG_PANEL") &&
            !debug::SUPERDEBUG_SAMPLENAME.empty() &&
            g_pre->name == debug::SUPERDEBUG_SAMPLENAME &&
            debug::SUPERDEBUG_BP > 0) {
            int locus = debug::SUPERDEBUG_BP;
            unsigned char amb = (!g_pre->Ambiguous.empty()) ? g_pre->Ambiguous[0] : 0;
            std::fprintf(stderr,
                         "[PRESWEEP_BUILD] ctx=%s locus=%d amb_code=0x%02x\n",
                         ctx.name.c_str(), locus, static_cast<int>(amb));
        }
        maybe_log_amb_code(ctx, "post_select");
    }

    const unsigned int max_transitions = 4096;
    const unsigned int max_missing = 4096;
    compute_job job(ctx.V, ctx.Gset, ctx.H, max_transitions, max_missing,
                    ctx.enable_supersites ? &ctx.ss_context.super_sites : nullptr,
                    ctx.enable_supersites ? &ctx.ss_context.locus_to_super_idx : nullptr,
                    ctx.enable_supersites ? &ctx.ss_context.super_site_var_index : nullptr);
    job.make(0, 0.0);
    if (env_true("SHAPEIT5_SUPERDEBUG_PANEL") &&
        (stage.label == "main4" || stage.label == "main5")) {
        maybe_dump_panel_anchor(ctx, job, stage);
    }

    IterationResult res;
    res.stage = stage.type;
    res.label = stage.label;
    res.k_sizes.reserve(job.Kstates.size());
    res.window_prob_sum.reserve(job.Kstates.size());

    for (const auto& ks : job.Kstates) res.k_sizes.push_back(static_cast<int>(ks.size()));
    
    // Store full Kstates for exact parity check
    res.Kstates = job.Kstates;

    genotype* sample = ctx.Gset.vecG[0];
    int outcome = 0;
    for (int w = 0; w < job.size(); ++w) {
        haplotype_segment_single HS(sample, ctx.H.H_opt_hap, job.Kstates[w], job.Windows.W[w], ctx.M,
                                    ctx.enable_supersites ? &ctx.ss_context.super_sites : nullptr,
                                    ctx.enable_supersites ? &ctx.ss_context.is_super_site : nullptr,
                                    ctx.enable_supersites ? &ctx.ss_context.locus_to_super_idx : nullptr,
                                    ctx.enable_supersites ? ctx.ss_context.packed_codes.data() : nullptr,
                                    ctx.enable_supersites ? ctx.ss_context.packed_codes.size() : 0,
                                    ctx.enable_supersites ? &ctx.ss_context.super_site_var_index : nullptr);
        HS.forward();
        outcome = HS.backward(job.T, job.M,
                              ctx.enable_supersites ? &job.SC : nullptr,
                              ctx.enable_supersites ? &job.anchor_has_missing : nullptr,
                              ctx.enable_supersites ? &job.supersite_sc_offset : nullptr);
        if (outcome != 0) {
            haplotype_segment_double HD(sample, ctx.H.H_opt_hap, job.Kstates[w], job.Windows.W[w], ctx.M,
                                        ctx.enable_supersites ? &ctx.ss_context.super_sites : nullptr,
                                        ctx.enable_supersites ? &ctx.ss_context.is_super_site : nullptr,
                                        ctx.enable_supersites ? &ctx.ss_context.locus_to_super_idx : nullptr,
                                        ctx.enable_supersites ? ctx.ss_context.packed_codes.data() : nullptr,
                                        ctx.enable_supersites ? ctx.ss_context.packed_codes.size() : 0,
                                        ctx.enable_supersites ? &ctx.ss_context.super_site_var_index : nullptr);
            HD.forward();
            outcome = HD.backward(job.T, job.M,
                                  ctx.enable_supersites ? &job.SC : nullptr,
                                  ctx.enable_supersites ? &job.anchor_has_missing : nullptr,
                                  ctx.enable_supersites ? &job.supersite_sc_offset : nullptr);
            res.window_prob_sum.push_back(HD.probSumT);
        } else {
            res.window_prob_sum.push_back(HS.probSumT);
        }
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
    } else {
        g->setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
    }

    switch (stage.type) {
        case StageType::Burn:
            g->sample(job.T, job.M);
            break;
        case StageType::Prune:
            g->sample(job.T, job.M);
            if (std::getenv("SHAPEIT5_DETAILED_ITERATION_TRACE")) {
                unsigned char amb0 = g->Ambiguous.empty() ? 0 : g->Ambiguous[0];
                std::fprintf(stderr, "[ITER_TRACE] %s prune pre-merge: n_segments=%d n_transitions=%d amb0=0x%02x\n",
                             ctx.name.c_str(), g->n_segments, g->n_transitions, static_cast<int>(amb0));
            }
            if (g->n_transitions > 0) {
                std::vector<bool> merge_flags;
                g->mapMerges(job.T, 0.95, merge_flags);
                g->performMerges(job.T, merge_flags);
            }
            if (std::getenv("SHAPEIT5_DETAILED_ITERATION_TRACE")) {
                unsigned char amb0 = g->Ambiguous.empty() ? 0 : g->Ambiguous[0];
                std::fprintf(stderr, "[ITER_TRACE] %s prune post-merge: n_segments=%d n_transitions=%d amb0=0x%02x\n",
                             ctx.name.c_str(), g->n_segments, g->n_transitions, static_cast<int>(amb0));
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
    maybe_log_amb_code(ctx, "post_transpose");
    return res;
}

static bool anchor_haplotype_parity(const MiniContext& bial, const MiniContext& supersite, int repeat_factor) {
    const genotype* gb = bial.Gset.vecG[0];
    const genotype* gs = supersite.Gset.vecG[0];
    const int anchors_per_block = 5;
    const int anchors_total = anchors_per_block * repeat_factor;
    for (int i = 0; i < anchors_total; ++i) {
        const int b_idx = i;
        const int ss_idx = i * 2; // anchors at even positions in supersite layout
        unsigned char bb = gb->Variants[DIV2(b_idx)];
        unsigned char sb = gs->Variants[DIV2(ss_idx)];
        bool b_h0 = VAR_GET_HAP0(MOD2(b_idx), bb);
        bool b_h1 = VAR_GET_HAP1(MOD2(b_idx), bb);
        bool s_h0 = VAR_GET_HAP0(MOD2(ss_idx), sb);
        bool s_h1 = VAR_GET_HAP1(MOD2(ss_idx), sb);
        if (b_h0 != s_h0 || b_h1 != s_h1) {
            std::cerr << "Anchor mismatch at locus " << b_idx << ": "
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
    TEST_INIT("test_supersite_expansion_epochs");

    std::cout << "======================================================================" << std::endl;
    std::cout << "Supersite Expansion 15-Epoch Parity Test" << std::endl;
    std::cout << "======================================================================" << std::endl;

    const unsigned int ref_samples = 64; // 128 donor haplotypes
    const std::vector<Scenario> scenarios = {
        {1, false, ""},   // original micro test
        {2, false, "x2"}, // expanded block to exercise repeated anchors
        {4, false, "x4"}, // moderate expansion
        {8, true,  "x8"}, // repeated pattern to force multiple segments
        {16, true, "x16"}, // heavy expansion for deeper epoch coverage
        {32, true, "x32"}, // larger expansion
        {64, true, "x64"}, // larger expansion
        {128, true, "x128"} // largest expansion
    };

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

    for (const auto& scenario : scenarios) {
        const int repeat_factor = scenario.repeat_factor;
        const bool require_multiseg = scenario.require_multiseg;
        const auto make_name = [&](bool supersite) {
            std::string base = supersite ? "supersite" : "biallelic";
            if (!scenario.name_suffix.empty()) return base + "_" + scenario.name_suffix;
            if (require_multiseg) return base + "_x" + std::to_string(repeat_factor);
            return base;
        };

        std::string scenario_name = "repeat_factor_" + std::to_string(repeat_factor);
        std::string scenario_desc = std::to_string(repeat_factor) + "x expansion" +
                                    (require_multiseg ? " (multi-segment)" : "");
        SCENARIO_START(scenario_name, scenario_desc);
        TEST_CONTEXT(scenario_name + " initialization");

        MiniContext ctx_bial = build_context(false, ref_samples, repeat_factor, require_multiseg,
                                             make_name(false));
        MiniContext ctx_ss = build_context(true, ref_samples, repeat_factor, require_multiseg,
                                           make_name(true));

        std::cout << "  Biallelic haplotypes: " << ctx_bial.H.n_hap
                  << " | Supersite haplotypes: " << ctx_ss.H.n_hap << std::endl;
        std::cout << "  Supersite count: " << ctx_ss.ss_context.super_sites.size() << std::endl;

        const unsigned int seed_base = 20250110;
        double max_prob_diff_observed = 0.0;
        std::string max_prob_diff_stage;
        for (size_t iter = 0; iter < schedule.size(); ++iter) {
            const StageDef& stage = schedule[iter];
            ITERATION(iter + 1, schedule.size(), stage.label);
            TEST_CONTEXT(scenario_name + ", iteration " + std::to_string(iter + 1) +
                        "/" + std::to_string(schedule.size()) + " (" + stage.label + ")");

            IterationResult res_bial = run_iteration(ctx_bial, stage, seed_base + iter);
            IterationResult res_sup = run_iteration(ctx_ss, stage, seed_base + iter);

            // Iteration-level haplotype tracing for divergence detection
            if (std::getenv("SHAPEIT5_DETAILED_ITERATION_TRACE")) {
                std::fprintf(stderr, "\n[ITER_TRACE] Iteration %zu (%s):\n",
                             iter+1, stage.label.c_str());
                const genotype* gb = ctx_bial.Gset.vecG[0];
                const genotype* gs = ctx_ss.Gset.vecG[0];

                // Log first 5 anchor haplotypes
                const bool log_all_anchors = env_true("SHAPEIT5_LOG_ALL_ANCHORS");
                const int anchors_to_log = log_all_anchors ? (repeat_factor * 5)
                                                           : std::min(5, repeat_factor * 5);
                for (int i = 0; i < anchors_to_log; ++i) {
                    int b_idx = i;
                    int ss_idx = i * 2;
                    unsigned char bb = gb->Variants[DIV2(b_idx)];
                    unsigned char sb = gs->Variants[DIV2(ss_idx)];
                    bool b_h0 = VAR_GET_HAP0(MOD2(b_idx), bb);
                    bool b_h1 = VAR_GET_HAP1(MOD2(b_idx), bb);
                    bool s_h0 = VAR_GET_HAP0(MOD2(ss_idx), sb);
                    bool s_h1 = VAR_GET_HAP1(MOD2(ss_idx), sb);

                    std::fprintf(stderr, "  anchor[%d]: bial=%d|%d ss=%d|%d %s\n",
                                 i, (int)b_h0, (int)b_h1, (int)s_h0, (int)s_h1,
                                 (b_h0 == s_h0 && b_h1 == s_h1) ? "✓" : "MISMATCH");
                }
            }

            const int k_diff = max_int_diff(res_bial.k_sizes, res_sup.k_sizes);
            if (k_diff != 0) {
                std::string reason = "K-state size divergence during " + stage.label +
                                    " (max delta=" + std::to_string(k_diff) + ")";
                std::cerr << reason << std::endl;
                SCENARIO_FAIL(scenario_name, reason);
                std::exit(EXIT_FAILURE);
            }

            // Check exact K-state parity
            if (res_bial.Kstates.size() != res_sup.Kstates.size()) {
                 std::string reason = "K-state window count mismatch during " + stage.label +
                                     ": " + std::to_string(res_bial.Kstates.size()) +
                                     " vs " + std::to_string(res_sup.Kstates.size());
                 std::cerr << reason << std::endl;
                 SCENARIO_FAIL(scenario_name, reason);
                 std::exit(EXIT_FAILURE);
            }
            for(size_t w=0; w<res_bial.Kstates.size(); ++w) {
                const auto& kb = res_bial.Kstates[w];
                const auto& ks = res_sup.Kstates[w];
                if (kb.size() != ks.size() || !std::equal(kb.begin(), kb.end(), ks.begin())) {
                    std::string reason = "K-state content divergence at window " + std::to_string(w) +
                                        " during " + stage.label;
                    std::cerr << reason << std::endl;
                    std::cerr << "  Biallelic: ";
                    for(auto x : kb) std::cerr << x << " ";
                    std::cerr << "\n  Supersite: ";
                    for(auto x : ks) std::cerr << x << " ";
                    std::cerr << "\n";
                    SCENARIO_FAIL(scenario_name, reason);
                    std::exit(EXIT_FAILURE);
                }
            }

            if (!anchor_haplotype_parity(ctx_bial, ctx_ss, repeat_factor)) {
                std::string reason = "Anchor haplotypes diverged during " + stage.label;
                std::cerr << reason << std::endl;
                SCENARIO_FAIL(scenario_name, reason);
                std::exit(EXIT_FAILURE);
            }

            const double prob_diff = max_abs_diff(res_bial.window_prob_sum, res_sup.window_prob_sum);
            if (prob_diff > max_prob_diff_observed) {
                max_prob_diff_observed = prob_diff;
                max_prob_diff_stage = stage.label;
            }

            // Trace ambiguous mask + segment count after each iteration when detailed tracing is enabled
            unsigned char amb0_b = ctx_bial.Gset.vecG[0]->Ambiguous.empty() ? 0 : ctx_bial.Gset.vecG[0]->Ambiguous[0];
            unsigned char amb0_s = ctx_ss.Gset.vecG[0]->Ambiguous.empty() ? 0 : ctx_ss.Gset.vecG[0]->Ambiguous[0];
            std::fprintf(stderr,
                         "[ITER_TRACE] amb_state post-%s: bial(seg=%d,amb0=0x%02x) supersite(seg=%d,amb0=0x%02x)\n",
                         stage.label.c_str(),
                         ctx_bial.Gset.vecG[0]->n_segments, static_cast<int>(amb0_b),
                         ctx_ss.Gset.vecG[0]->n_segments, static_cast<int>(amb0_s));
        }

        if (!max_prob_diff_stage.empty()) {
            std::cout << "Max probSumT delta observed: " << max_prob_diff_observed
                      << " during " << max_prob_diff_stage << std::endl;
        }

        std::string details = std::to_string(schedule.size()) + " iterations";
        if (!max_prob_diff_stage.empty()) {
            details += ", max_delta=" + std::to_string(max_prob_diff_observed);
        }
        SCENARIO_PASS(scenario_name, details);
    }

    std::cout << "\nAll scenarios completed without divergence." << std::endl;
    TEST_SUMMARY();
    return 0;
}

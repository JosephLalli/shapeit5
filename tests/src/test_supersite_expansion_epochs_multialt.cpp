/*******************************************************************************
 * Supersite Expansion Epochs Test with True Multiallelic Anchors (Gap D)
 *
 * This test extends test_supersite_expansion_epochs.cpp by using GENUINE
 * multiallelic anchors where the sample is heterozygous between DISTINCT
 * ALT classes (e.g., ALT1|ALT2, ALT2|ALT3) instead of dummy 0/0 siblings.
 *
 * Key differences from test_supersite_expansion_epochs.cpp:
 * - Sample genotypes include true multiallelic HET states (ALT1|ALT2)
 * - Panel patterns ensure each ALT class has distinct donor haplotypes
 * - Validates mutual exclusivity at supersite anchors across epochs
 * - Tracks per-anchor (c0,c1) class codes via getSampleSuperSiteAlleleCode()
 * - Runs shorter schedule (5 burn + 2 prune + 5 main = 12 iterations)
 *
 * Purpose: Stress-test Phase 3 SC-based imputation with segment pruning in
 * realistic multiallelic settings. This may expose bugs not visible when
 * using dummy siblings with 0|0 genotypes.
 ******************************************************************************/

#include <algorithm>
#include <array>
#include <cassert>
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

// Panel patterns for multiallelic test
// Each supersite has 3 splits (ALT1, ALT2, ALT3)
// Donor haplotypes segregate into distinct blocks for each ALT class
struct PanelPatterns {
    std::vector<std::array<int,8>> pattern;
};

static PanelPatterns make_multiallelic_panel_patterns() {
    PanelPatterns patterns;
    // ss1 at bp=1000: ALT1(C), ALT2(G), ALT3(T)
    // Donors 0,1: prefer ALT1 (C)
    // Donors 2,3: prefer ALT2 (G)
    // Donors 4,5: prefer ALT3 (T)
    // Donors 6,7: REF
    patterns.pattern = {
        std::array<int,8>{1,1, 0,0, 0,0, 0,0}, // ss1_A_C (ALT1): donors 0,1 carry
        std::array<int,8>{0,0, 1,1, 0,0, 0,0}, // ss1_A_G (ALT2): donors 2,3 carry
        std::array<int,8>{0,0, 0,0, 1,1, 0,0}, // ss1_A_T (ALT3): donors 4,5 carry

        // ss2 at bp=2000: similar pattern but rotated
        std::array<int,8>{0,0, 1,1, 0,0, 0,0}, // ss2_A_C: donors 2,3 carry
        std::array<int,8>{0,0, 0,0, 1,1, 0,0}, // ss2_A_G: donors 4,5 carry
        std::array<int,8>{0,0, 0,0, 0,0, 1,1}, // ss2_A_T: donors 6,7 carry

        // v3000: simple biallelic
        std::array<int,8>{1,0, 0,1, 1,0, 0,1}  // v3000
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

static void init_multiallelic_variant_map(variant_map& V, int repeat) {
    int idx = 0;
    const int block = 5000;
    for (int r = 0; r < repeat; ++r) {
        const int offset = r * block;
        // Supersite 1: bp=1000, 3 ALTs (C, G, T)
        V.push(make_var("1", 1000 + offset, "ss1_A_C", "A", "C", idx++));
        V.push(make_var("1", 1000 + offset, "ss1_A_G", "A", "G", idx++));
        V.push(make_var("1", 1000 + offset, "ss1_A_T", "A", "T", idx++));

        // Supersite 2: bp=2000, 3 ALTs (C, G, T)
        V.push(make_var("1", 2000 + offset, "ss2_A_C", "A", "C", idx++));
        V.push(make_var("1", 2000 + offset, "ss2_A_G", "A", "G", idx++));
        V.push(make_var("1", 2000 + offset, "ss2_A_T", "A", "T", idx++));

        // Simple biallelic variant
        V.push(make_var("1", 3000 + offset, "v3000", "A", "T", idx++));
    }
    for (size_t i = 0; i < V.size(); ++i) {
        V.vec_pos[i]->cm = 0.01 * static_cast<double>(i + 1);
    }
}

// Sample phases: ALT1|ALT2 at ss1, ALT2|ALT3 at ss2, REF|ALT at v3000
static std::vector<PhaseCode> make_multiallelic_phases(int repeat) {
    std::vector<PhaseCode> base = {
        // ss1: ALT1|ALT2 means split0=ALT|REF, split1=REF|ALT, split2=REF|REF
        ALT_REF,  // ss1_A_C: hap0=ALT1
        REF_ALT,  // ss1_A_G: hap1=ALT2
        REF_REF,  // ss1_A_T: neither

        // ss2: ALT2|ALT3 means split0=REF|REF, split1=ALT|REF, split2=REF|ALT
        REF_REF,  // ss2_A_C: neither
        ALT_REF,  // ss2_A_G: hap0=ALT2
        REF_ALT,  // ss2_A_T: hap1=ALT3

        // v3000: simple het
        REF_ALT   // v3000: REF|ALT
    };
    return repeat_pattern(base, repeat);
}

static MiniContext build_multiallelic_context(unsigned int n_ref_samples, int repeat_factor = 1) {
    const PanelPatterns patterns = make_multiallelic_panel_patterns();
    MiniContext ctx;
    ctx.name = "multiallelic";
    ctx.panel_pattern = repeat_pattern(patterns.pattern, repeat_factor);

    init_multiallelic_variant_map(ctx.V, repeat_factor);
    const std::vector<PhaseCode> phases = make_multiallelic_phases(repeat_factor);
    genotype sample = make_sample_from_phases(phases, "multiallelic_sample");

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

    ctx.ss_context = build_supersites(ctx.V, ctx.H);
    apply_supersite_pbwt_guards(ctx.H, ctx.ss_context, ctx.V.size());

    genotype* g0 = ctx.Gset.vecG[0];
    g0->setSuperSiteContext(&ctx.ss_context.super_sites,
                            &ctx.ss_context.locus_to_super_idx,
                            &ctx.ss_context.super_site_var_index,
                            nullptr, nullptr, nullptr);
    g0->snapshotSupersiteBaseClasses(ctx.ss_context.super_sites,
                                     ctx.ss_context.super_site_var_index);
    g0->build();

    ctx.M = make_hmm_params_from_variant_map(ctx.V, ctx.H.n_hap);
    ctx.M.markSuperSiteSiblings(ctx.ss_context.super_sites, ctx.ss_context.locus_to_super_idx);

    return ctx;
}

// Check mutual exclusivity: exactly one ALT per haplotype across all splits of a supersite
static bool check_mutual_exclusivity(const genotype& G, const SuperSite& ss,
                                     const std::vector<int>& super_site_var_index,
                                     int ss_idx) {
    for (int hap = 0; hap < 2; ++hap) {
        int alt_count = 0;
        for (uint32_t i = 0; i < ss.var_count; ++i) {
            int v_idx = super_site_var_index[ss.var_start + i];
            unsigned char byte = G.Variants[DIV2(v_idx)];
            bool carries = (hap == 0) ? VAR_GET_HAP0(MOD2(v_idx), byte)
                                      : VAR_GET_HAP1(MOD2(v_idx), byte);
            if (carries) alt_count++;
        }
        if (alt_count > 1) {
            std::cerr << "[MUTUAL_EXCLUSIVITY_VIOLATION] ss_idx=" << ss_idx
                      << " hap=" << hap << " alt_count=" << alt_count << std::endl;
            return false;
        }
    }
    return true;
}

// Check class code consistency with projected haplotype bits
static bool check_class_consistency(const genotype& G, const SuperSite& ss,
                                    const std::vector<int>& super_site_var_index,
                                    int ss_idx) {
    uint8_t c0 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 0);
    uint8_t c1 = getSampleSuperSiteAlleleCode(&G, ss, super_site_var_index, 1);

    // Verify each haplotype code matches the split pattern
    for (int hap = 0; hap < 2; ++hap) {
        uint8_t code = (hap == 0) ? c0 : c1;

        if (code == SUPERSITE_CODE_MISSING) continue; // Missing is allowed
        if (code == SUPERSITE_CODE_REF) {
            // All splits should be REF for this hap
            for (uint32_t i = 0; i < ss.var_count; ++i) {
                int v_idx = super_site_var_index[ss.var_start + i];
                unsigned char byte = G.Variants[DIV2(v_idx)];
                bool carries = (hap == 0) ? VAR_GET_HAP0(MOD2(v_idx), byte)
                                          : VAR_GET_HAP1(MOD2(v_idx), byte);
                if (carries) {
                    std::cerr << "[CLASS_CONSISTENCY_VIOLATION] ss_idx=" << ss_idx
                              << " hap=" << hap << " code=REF but split " << i
                              << " carries ALT" << std::endl;
                    return false;
                }
            }
        } else {
            // code-1 should be the index of the split carrying ALT
            uint32_t expected_split = code - 1;
            for (uint32_t i = 0; i < ss.var_count; ++i) {
                int v_idx = super_site_var_index[ss.var_start + i];
                unsigned char byte = G.Variants[DIV2(v_idx)];
                bool carries = (hap == 0) ? VAR_GET_HAP0(MOD2(v_idx), byte)
                                          : VAR_GET_HAP1(MOD2(v_idx), byte);
                if (i == expected_split && !carries) {
                    std::cerr << "[CLASS_CONSISTENCY_VIOLATION] ss_idx=" << ss_idx
                              << " hap=" << hap << " code=" << (int)code
                              << " but split " << i << " does not carry ALT" << std::endl;
                    return false;
                }
                if (i != expected_split && carries) {
                    std::cerr << "[CLASS_CONSISTENCY_VIOLATION] ss_idx=" << ss_idx
                              << " hap=" << hap << " code=" << (int)code
                              << " but split " << i << " unexpectedly carries ALT" << std::endl;
                    return false;
                }
            }
        }
    }
    return true;
}

static IterationResult run_iteration(MiniContext& ctx, StageDef stage, unsigned int rng_seed) {
    rng.setSeed(rng_seed);

    ctx.H.select();

    // Check PBWT didn't select siblings
    if (!ctx.ss_context.super_sites.empty()) {
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

    ctx.ss_context = build_supersites(ctx.V, ctx.H);
    ctx.M.markSuperSiteSiblings(ctx.ss_context.super_sites, ctx.ss_context.locus_to_super_idx);
    apply_supersite_pbwt_guards(ctx.H, ctx.ss_context, ctx.V.size());

    genotype* g_pre = ctx.Gset.vecG[0];
    g_pre->setSuperSiteContext(&ctx.ss_context.super_sites,
                               &ctx.ss_context.locus_to_super_idx,
                               &ctx.ss_context.super_site_var_index,
                               nullptr, nullptr, nullptr);
    g_pre->snapshotSupersiteClasses(ctx.ss_context.super_sites,
                                    ctx.ss_context.super_site_var_index);

    const unsigned int max_transitions = 4096;
    const unsigned int max_missing = 4096;
    compute_job job(ctx.V, ctx.Gset, ctx.H, max_transitions, max_missing,
                    &ctx.ss_context.super_sites,
                    &ctx.ss_context.locus_to_super_idx,
                    &ctx.ss_context.super_site_var_index);
    job.make(0, 0.0);

    IterationResult res;
    res.stage = stage.type;
    res.label = stage.label;
    res.k_sizes.reserve(job.Kstates.size());
    res.window_prob_sum.reserve(job.Kstates.size());

    for (const auto& ks : job.Kstates) res.k_sizes.push_back(static_cast<int>(ks.size()));
    res.Kstates = job.Kstates;

    genotype* sample = ctx.Gset.vecG[0];
    int outcome = 0;
    for (int w = 0; w < job.size(); ++w) {
        haplotype_segment_single HS(sample, ctx.H.H_opt_hap, job.Kstates[w], job.Windows.W[w], ctx.M,
                                    &ctx.ss_context.super_sites,
                                    &ctx.ss_context.is_super_site,
                                    &ctx.ss_context.locus_to_super_idx,
                                    ctx.ss_context.packed_codes.data(),
                                    ctx.ss_context.packed_codes.size(),
                                    &ctx.ss_context.super_site_var_index);
        HS.forward();
        outcome = HS.backward(job.T, job.M, &job.SC, &job.anchor_has_missing, &job.supersite_sc_offset);
        if (outcome != 0) {
            haplotype_segment_double HD(sample, ctx.H.H_opt_hap, job.Kstates[w], job.Windows.W[w], ctx.M,
                                        &ctx.ss_context.super_sites,
                                        &ctx.ss_context.is_super_site,
                                        &ctx.ss_context.locus_to_super_idx,
                                        ctx.ss_context.packed_codes.data(),
                                        ctx.ss_context.packed_codes.size(),
                                        &ctx.ss_context.super_site_var_index);
            HD.forward();
            outcome = HD.backward(job.T, job.M, &job.SC, &job.anchor_has_missing, &job.supersite_sc_offset);
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
    g->setSuperSiteContext(&ctx.ss_context.super_sites,
                           &ctx.ss_context.locus_to_super_idx,
                           &ctx.ss_context.super_site_var_index,
                           &job.SC,
                           &job.anchor_has_missing,
                           &job.supersite_sc_offset);

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

} // namespace

int main() {
    TEST_INIT("test_supersite_expansion_epochs_multialt");

    std::cout << "======================================================================" << std::endl;
    std::cout << "Supersite Expansion Epochs Test with Multiallelic Anchors (Gap D)" << std::endl;
    std::cout << "======================================================================" << std::endl;

    const unsigned int ref_samples = 64; // 128 donor haplotypes (matching existing epoch tests)
    const int repeat_factor = 1;

    const std::vector<StageDef> schedule = {
        {StageType::Burn,  "burn1"},
        {StageType::Burn,  "burn2"},
        {StageType::Burn,  "burn3"},
        {StageType::Burn,  "burn4"},
        {StageType::Burn,  "burn5"},
        {StageType::Prune, "prune1"},
        {StageType::Burn,  "burn6"},
        {StageType::Prune, "prune2"},
        {StageType::Main,  "main1"},
        {StageType::Main,  "main2"},
        {StageType::Main,  "main3"},
        {StageType::Main,  "main4"},
        {StageType::Main,  "main5"}
    };

    SCENARIO_START("multiallelic_epochs", "Multiallelic anchors across 12 epochs");
    TEST_CONTEXT("multiallelic_epochs initialization");

    MiniContext ctx = build_multiallelic_context(ref_samples, repeat_factor);

    std::cout << "  Total haplotypes: " << ctx.H.n_hap << std::endl;
    std::cout << "  Supersite count: " << ctx.ss_context.super_sites.size() << std::endl;
    std::cout << "  Total variants: " << ctx.V.size() << std::endl;

    const unsigned int seed_base = 20251201;

    for (size_t iter = 0; iter < schedule.size(); ++iter) {
        const StageDef& stage = schedule[iter];
        ITERATION(iter + 1, schedule.size(), stage.label);
        TEST_CONTEXT("multiallelic_epochs, iteration " + std::to_string(iter + 1) +
                    "/" + std::to_string(schedule.size()) + " (" + stage.label + ")");

        IterationResult res = run_iteration(ctx, stage, seed_base + iter);

        // Validate mutual exclusivity and class consistency at all supersites
        const genotype* g = ctx.Gset.vecG[0];
        for (size_t ss_idx = 0; ss_idx < ctx.ss_context.super_sites.size(); ++ss_idx) {
            const SuperSite& ss = ctx.ss_context.super_sites[ss_idx];

            if (!check_mutual_exclusivity(*g, ss, ctx.ss_context.super_site_var_index, ss_idx)) {
                std::string reason = "Mutual exclusivity violated at ss_idx=" +
                                    std::to_string(ss_idx) + " during " + stage.label;
                std::cerr << reason << std::endl;
                SCENARIO_FAIL("multiallelic_epochs", reason);
                std::exit(EXIT_FAILURE);
            }

            if (!check_class_consistency(*g, ss, ctx.ss_context.super_site_var_index, ss_idx)) {
                std::string reason = "Class consistency violated at ss_idx=" +
                                    std::to_string(ss_idx) + " during " + stage.label;
                std::cerr << reason << std::endl;
                SCENARIO_FAIL("multiallelic_epochs", reason);
                std::exit(EXIT_FAILURE);
            }
        }

        // Optional: Trace class codes at each epoch
        if (env_true("SHAPEIT5_TEST_TRACE")) {
            std::fprintf(stderr, "\n[CLASS_TRACE] Iteration %zu (%s):\n",
                        iter+1, stage.label.c_str());
            for (size_t ss_idx = 0; ss_idx < ctx.ss_context.super_sites.size(); ++ss_idx) {
                const SuperSite& ss = ctx.ss_context.super_sites[ss_idx];
                uint8_t c0 = getSampleSuperSiteAlleleCode(g, ss,
                                ctx.ss_context.super_site_var_index, 0);
                uint8_t c1 = getSampleSuperSiteAlleleCode(g, ss,
                                ctx.ss_context.super_site_var_index, 1);
                std::fprintf(stderr, "  ss[%zu]: c0=%d c1=%d\n",
                            ss_idx, (int)c0, (int)c1);
            }
        }
    }

    SCENARIO_PASS("multiallelic_epochs", std::to_string(schedule.size()) +
                 " iterations, mutual exclusivity maintained");

    std::cout << "\nAll epochs completed with mutual exclusivity and class consistency." << std::endl;
    TEST_SUMMARY();
    return 0;
}

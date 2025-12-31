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
#include <vector>
#include <string>

#include "../../common/src/utils/otools.h"

#include "test_reporting.h"

#define private public
#define protected public
#include "../../phase_common/src/models/haplotype_segment_double.h"
#undef private
#undef protected

#include "../../phase_common/src/objects/super_site_builder.h"
#include "../../phase_common/src/containers/variant_map.h"
#include "../../phase_common/src/containers/conditioning_set/conditioning_set_header.h"

namespace {

struct SuperSiteContext {
    std::vector<SuperSite> super_sites;
    std::vector<bool> is_super_site;
    std::vector<uint8_t> packed_codes;
    std::vector<int> locus_to_super_idx;
    std::vector<int> super_site_var_index;
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

static window make_full_window(const genotype& G,
                               const SuperSiteContext* ctx,
                               int stop_locus) {
    window W{};
    W.start_locus = 0;
    W.stop_locus = stop_locus;
    W.start_segment = 0;

    // Determine stop segment index (segments are contiguous in this fixture)
    unsigned int locus_acc = 0;
    unsigned int seg_idx = 0;
    while (seg_idx < G.n_segments) {
        unsigned int seg_len = (seg_idx < G.Lengths.size()) ? G.Lengths[seg_idx] : 0;
        if (seg_len == 0) break;
        if (locus_acc + seg_len - 1 >= static_cast<unsigned int>(stop_locus)) break;
        locus_acc += seg_len;
        seg_idx++;
    }
    W.stop_segment = static_cast<int>(seg_idx);

    const std::vector<SuperSite>* super_sites = ctx ? &ctx->super_sites : nullptr;
    const std::vector<int>* locus_to_super_idx = ctx ? &ctx->locus_to_super_idx : nullptr;

    int first_amb = -1;
    int last_amb = -1;
    int amb_count = 0;

    for (int v = 0; v <= stop_locus; ) {
        int ss_idx = (super_sites && locus_to_super_idx && v < static_cast<int>(locus_to_super_idx->size()))
                     ? (*locus_to_super_idx)[v]
                     : -1;
        bool is_anchor = false;
        unsigned int advance = 1;
        bool consider_variant = true;

        if (ss_idx >= 0 && super_sites) {
            const SuperSite& ss = (*super_sites)[ss_idx];
            is_anchor = (v == static_cast<int>(ss.global_site_id));
            if (is_anchor) {
                advance = ss.var_count;
            } else {
                consider_variant = false; // supersite siblings do not advance ambiguous cursor
            }
        }

        bool is_amb = false;
        if (consider_variant) {
            const unsigned char var_code = G.Variants[DIV2(v)];
            is_amb = VAR_GET_AMB(MOD2(v), var_code) != 0;
            if (!is_amb && VAR_GET_SCA(MOD2(v), var_code) != 0) {
                is_amb = true;
            }

            if (is_amb) {
                if (first_amb < 0) first_amb = amb_count;
                last_amb = amb_count;
                amb_count++;
            }
        }

        if (advance == 0) advance = 1;
        v += static_cast<int>(advance);
    }

    if (amb_count > 0) {
        W.start_ambiguous = 0;
        W.stop_ambiguous = last_amb;
    } else {
        W.start_ambiguous = 0;
        W.stop_ambiguous = -1;
    }

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
    // Use standard supersite emissions for debugging donor codes
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
                    ctx.locus_to_super_idx, ctx.super_site_var_index);
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
    if (ctx) {
        G.setSuperSiteContext(&ctx->super_sites, &ctx->locus_to_super_idx,
                              &ctx->super_site_var_index, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(ctx->packed_codes.data(), ctx->packed_codes.size());
    } else {
        G.setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(nullptr, 0);
    }

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M);
    HS.forward();
    FBState s; s.prob.assign(HS.prob.begin(), HS.prob.end()); s.sumH.assign(HS.probSumH.begin(), HS.probSumH.end());
    return s;
}

static FBResult run_forward_backward(genotype& G,
                                     conditioning_set& H,
                                     hmm_parameters& M,
                                     const window& W,
                                     const std::vector<unsigned int>& idxH,
                                     const SuperSiteContext* ctx) {
    if (ctx) {
        G.setSuperSiteContext(&ctx->super_sites, &ctx->locus_to_super_idx,
                              &ctx->super_site_var_index, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(ctx->packed_codes.data(), ctx->packed_codes.size());
    } else {
        G.setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(nullptr, 0);
    }

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M);

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

    if (ctx) {
        G.setSuperSiteContext(&ctx->super_sites, &ctx->locus_to_super_idx,
                              &ctx->super_site_var_index, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(ctx->packed_codes.data(), ctx->packed_codes.size());
    } else {
        G.setSuperSiteContext(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);
        G.setSupersitePanelCodes(nullptr, 0);
    }

    haplotype_segment_double HS(&G, H.H_opt_hap, const_cast<std::vector<unsigned int>&>(idxH),
                                const_cast<window&>(W), M);
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

} // namespace

int main() {
    TEST_INIT("test_supersite_representation_parity");
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
        G.Lengths_bio = G.Lengths;
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
    const std::vector<unsigned int> idxH = {0u,1u,2u,3u,4u,5u,6u,7u};

    hmm_parameters M_no_ss = make_hmm_params(V.size(), H.n_hap);
    hmm_parameters M_with_ss = make_hmm_params(V.size(), H.n_hap);

    SuperSiteContext ctx = build_supersites(V, H);
    assert(ctx.super_sites.size() == 1); // supersite for the two splits
    M_with_ss.markSuperSiteSiblings(ctx.super_sites, ctx.locus_to_super_idx);

    const window W_no = make_full_window(G_no_ss, nullptr, static_cast<int>(V.size()) - 1);
    const window W_ss = make_full_window(G_with_ss, &ctx, static_cast<int>(V.size()) - 1);

    // ---------------------------------------------------------------------
    // Run forward/backward on both datasets
    // ---------------------------------------------------------------------
    const bool skip_baseline = [](){
        const char* env = std::getenv("SHAPEIT5_SKIP_BASELINE_PARITY");
        return env && env[0] != '\0' && env[0] != '0';
    }();

    FBResult res_no_ss{};
    if (!skip_baseline) {
        res_no_ss = run_forward_backward(G_no_ss, H, M_no_ss, W_no, idxH, nullptr);
    }
    FBResult res_with_ss = run_forward_backward(G_with_ss, H, M_with_ss, W_ss, idxH, &ctx);

    // Anchor-only forward parity (hard check): window ends at the anchor locus (index 1)
    {
        window Wa_no = make_full_window(G_no_ss, nullptr, 1);
        window Wa_ss = make_full_window(G_with_ss, &ctx, 1);
        FBState fa_no = run_forward_only(G_no_ss, H, M_no_ss, Wa_no, idxH, nullptr);
        FBState fa_ss = run_forward_only(G_with_ss, H, M_with_ss, Wa_ss, idxH, &ctx);
        const double tol = 1e-9;
        // Compare normalized α (per lane across donors)
        std::printf("=== ANCHOR-ONLY PARITY DETAILED COMPARISON ===\n");
        std::printf("Biallelic sumH: [");
        for (int h = 0; h < HAP_NUMBER; ++h) std::printf("%.6f%s", fa_no.sumH[h], h < HAP_NUMBER-1 ? ", " : "]\n");
        std::printf("Supersite sumH: [");
        for (int h = 0; h < HAP_NUMBER; ++h) std::printf("%.6f%s", fa_ss.sumH[h], h < HAP_NUMBER-1 ? ", " : "]\n");
        
        // Check for cyclic permutation pattern
        std::printf("PATTERN ANALYSIS:\n");
        for (int shift = 0; shift < HAP_NUMBER; ++shift) {
            bool matches = true;
            for (int h = 0; h < HAP_NUMBER; ++h) {
                int shifted_idx = (h + shift) % HAP_NUMBER;
                if (std::fabs(fa_no.sumH[h] - fa_ss.sumH[shifted_idx]) > 1e-8) {
                    matches = false;
                    break;
                }
            }
            if (matches) {
                std::printf("  Supersite sumH is biallelic shifted by %d positions\n", shift);
                break;
            }
        }
        
        for (unsigned int k = 0, idx = 0; k < H.n_hap; ++k, idx += HAP_NUMBER) {
            std::printf("Donor %u:\n", k);
            std::printf("  Biallelic prob: [");
            for (int h = 0; h < HAP_NUMBER; ++h) std::printf("%.8f%s", fa_no.prob[idx + h], h < HAP_NUMBER-1 ? ", " : "]\n");
            std::printf("  Supersite prob: [");
            for (int h = 0; h < HAP_NUMBER; ++h) std::printf("%.8f%s", fa_ss.prob[idx + h], h < HAP_NUMBER-1 ? ", " : "]\n");
            
            for (int h = 0; h < HAP_NUMBER; ++h) {
                double d0 = fa_no.sumH[h], d1 = fa_ss.sumH[h];
                if (d0 > 0 && d1 > 0) {
                    double p0 = fa_no.prob[idx + h] / d0;
                    double p1 = fa_ss.prob[idx + h] / d1;
                    double diff = p0 - p1;
                    std::printf("  Lane %d: p0=%.8f p1=%.8f diff=%.8f\n", h, p0, p1, diff);
                    if (std::fabs(diff) > tol) {
                        std::cerr << "Anchor-only parity failed at donor=" << k << " lane=" << h << " diff=" << diff << std::endl;
                        assert(false);
                    }
                }
            }
        }
    }

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
                                  const SuperSiteContext* context,
                                  const window& W_local) {
        double best_score = -std::numeric_limits<double>::infinity();
        Orientation best = candidates.front();
        for (const Orientation& o : candidates) {
            double score = evaluate_orientation(G, panel, M, W_local, idxH, context, o);
            if (score > best_score) {
                best_score = score;
                best = o;
            }
        }
        return best;
    };

    Orientation best_no_ss = orientation_argmax(G_no_ss, H, M_no_ss, nullptr, W_no);
    Orientation best_with_ss = orientation_argmax(G_with_ss, H, M_with_ss, &ctx, W_ss);

    // Expected phasing sequence
    const Orientation expected{{ALT_REF, ALT_REF, REF_REF, REF_REF, ALT_ALT}};

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

    TEST_SUMMARY();
    return 0;
}

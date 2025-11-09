#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <containers/bitmatrix.h>
#include <models/site_emission_types.h>
#include <models/super_site_accessor.h>
#include <objects/genotype/genotype_header.h>
#include <utils/otools.h>

class genotype;

class BiallelicEmissionAdapter {
public:
    explicit BiallelicEmissionAdapter(const genotype* G, bitmatrix* Hvar)
        : G_(G), Hvar_(Hvar) {}

    void build_view(int abs_locus, int curr_abs_ambiguous, SiteView& view) const;

    void build_match_mask(const SiteView& view,
                          unsigned int n_cond_haps,
                          int matrix_row,
                          MatchMask& mask) const;

private:
    const genotype* G_{nullptr};
    bitmatrix* Hvar_{nullptr};
};

class SupersiteEmissionAdapter {
public:
    SupersiteEmissionAdapter(const genotype* G,
                              const std::vector<SuperSite>* super_sites,
                              const std::vector<int>* locus_to_super_idx,
                              const std::vector<int>* super_site_var_index,
                              const uint8_t* panel_codes,
                                                            const std::vector<unsigned int>* cond_idx,
                                                            size_t panel_codes_size = 0)
        : G_(G),
          super_sites_(super_sites),
          locus_to_super_idx_(locus_to_super_idx),
          super_site_var_index_(super_site_var_index),
                    panel_codes_(panel_codes),
                    cond_idx_(cond_idx),
                    panel_codes_size_(panel_codes_size) {}

    bool build_view(int abs_locus,
                    int curr_abs_ambiguous,
                    SiteView& view) const;

    void build_match_mask(const SiteView& view,
                          unsigned int n_cond_haps,
                          bool use_anchor_split_semantics,
                          MatchMask& mask) const;

private:
    const genotype* G_{nullptr};
    const std::vector<SuperSite>* super_sites_{nullptr};
    const std::vector<int>* locus_to_super_idx_{nullptr};
    const std::vector<int>* super_site_var_index_{nullptr};
    const uint8_t* panel_codes_{nullptr};
    const std::vector<unsigned int>* cond_idx_{nullptr};
    size_t panel_codes_size_{0};
};

// =============================
// Inline implementations
// =============================

inline void BiallelicEmissionAdapter::build_view(int abs_locus,
                                                 int curr_abs_ambiguous,
                                                 SiteView& view) const {
    view.kind = SiteKind::Biallelic;
    view.supersite = nullptr;
    view.supersite_index = -1;
    view.locus = abs_locus;
    view.anchor_class = 0;
    view.sample_class0 = 0;
    view.sample_class1 = 0;

    if (!G_) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        view.sample_class0 = SUPERSITE_CODE_MISSING;
        view.sample_class1 = SUPERSITE_CODE_MISSING;
        view.amb_mask = 0u;
        return;
    }

    const unsigned char variant = G_->Variants[DIV2(abs_locus)];
    const bool is_missing = VAR_GET_MIS(MOD2(abs_locus), variant);
    const bool is_ambiguous = VAR_GET_AMB(MOD2(abs_locus), variant);

    if (is_missing) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        view.sample_class0 = SUPERSITE_CODE_MISSING;
        view.sample_class1 = SUPERSITE_CODE_MISSING;
        view.amb_mask = 0u;
        return;
    }

    if (is_ambiguous) {
        view.emit_kind = EmitKind::Amb;
        uint8_t amb_code = 0u;
        if (curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) {
            amb_code = G_->Ambiguous[curr_abs_ambiguous];
        }
        view.amb_mask = amb_code;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool wants_alt = ((amb_code >> h) & 1U) != 0;
            view.lane_class[h] = wants_alt ? 1u : 0u;
        }
        view.sample_class0 = 0u;
        view.sample_class1 = 1u;
        return;
    }

    view.emit_kind = EmitKind::Hom;
    const bool hap0_is_alt = VAR_GET_HAP0(MOD2(abs_locus), variant);
    const uint8_t allele = hap0_is_alt ? 1u : 0u;
    std::fill(std::begin(view.lane_class), std::end(view.lane_class), allele);
    view.sample_class0 = allele;
    view.sample_class1 = allele;
    view.amb_mask = 0u;
}

inline void BiallelicEmissionAdapter::build_match_mask(const SiteView& view,
                                                       unsigned int n_cond_haps,
                                                       int matrix_row,
                                                       MatchMask& mask) const {
    const std::size_t total_entries = static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER;
    mask.resize(total_entries);

    if (view.emit_kind == EmitKind::Mis) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), MatchMask::kMatch);
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), true);
        return;
    }

    if (!Hvar_) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(0));
        return;
    }

    std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), MatchMask::kMismatch);
    std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);

    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const bool donor_alt = Hvar_->get(matrix_row, k);
        const uint8_t donor_code = donor_alt ? 1u : 0u;
        const std::size_t base = static_cast<std::size_t>(k) * HAP_NUMBER;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool match = (donor_code == view.lane_class[h]);
            mask.by_donor_lane[base + h] = match ? MatchMask::kMatch : MatchMask::kMismatch;
            mask.any_match_lane[h] = mask.any_match_lane[h] || match;
        }
    }
}

inline bool SupersiteEmissionAdapter::build_view(int abs_locus,
                                                 int curr_abs_ambiguous,
                                                 SiteView& view) const {
    view.locus = abs_locus;
    view.anchor_class = 0;
    view.supersite_index = -1;
    std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
    view.sample_class0 = 0u;
    view.sample_class1 = 0u;

    if (!super_sites_ || !locus_to_super_idx_ || abs_locus < 0 || abs_locus >= static_cast<int>(locus_to_super_idx_->size())) {
        view.kind = SiteKind::Biallelic;
        view.supersite = nullptr;
        return false;
    }

    const int ss_idx = (*locus_to_super_idx_)[abs_locus];
    if (ss_idx < 0 || ss_idx >= static_cast<int>(super_sites_->size())) {
        view.kind = SiteKind::Biallelic;
        view.supersite = nullptr;
        return false;
    }

    const SuperSite& ss = (*super_sites_)[ss_idx];
    view.supersite = &ss;
    view.supersite_index = ss_idx;
    const bool is_anchor = (abs_locus == static_cast<int>(ss.global_site_id));
    view.kind = is_anchor ? SiteKind::SuperAnchor : SiteKind::SuperSibling;

    // Optional trace
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    if (tr && tr[0] != '\0' && tr[0] != '0') {
        std::fprintf(stdout, "build_view locus=%d ss_idx=%d is_anchor=%d var_count=%u panel_off=%u\n",
                     abs_locus, ss_idx, (int)is_anchor, (unsigned)ss.var_count, (unsigned)ss.panel_offset);
    }

    if (!is_anchor) {
        view.emit_kind = EmitKind::Mis;
        view.sample_class0 = SUPERSITE_CODE_MISSING;
        view.sample_class1 = SUPERSITE_CODE_MISSING;
        view.amb_mask = 0u;
        return true;
    }

    uint8_t c0 = SUPERSITE_CODE_MISSING;
    uint8_t c1 = SUPERSITE_CODE_MISSING;
    SSClass cls = SSClass::MIS;
    if (G_ && super_site_var_index_) {
        cls = classify_supersite(G_, ss, *super_site_var_index_, c0, c1);
    }

    switch (cls) {
        case SSClass::MIS:
            view.emit_kind = EmitKind::Mis;
            view.sample_class0 = SUPERSITE_CODE_MISSING;
            view.sample_class1 = SUPERSITE_CODE_MISSING;
            view.amb_mask = 0u;
            break;
        case SSClass::HOM:
            view.emit_kind = EmitKind::Hom;
            std::fill(std::begin(view.lane_class), std::end(view.lane_class), c0);
            view.sample_class0 = c0;
            view.sample_class1 = c0;
            view.amb_mask = 0u;
            break;
        case SSClass::AMB: {
            view.emit_kind = EmitKind::Amb;
            uint8_t amb_mask = 0u;
            if (G_ && curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) {
                amb_mask = G_->Ambiguous[curr_abs_ambiguous];
            }
            view.amb_mask = amb_mask;
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const bool choose_c1 = ((amb_mask >> h) & 1U) != 0;
                view.lane_class[h] = choose_c1 ? c1 : c0;
            }
            view.sample_class0 = c0;
            view.sample_class1 = c1;
            break;
        }
    }

    // Record anchor ALT class (1..n_alts) for downstream split semantics
    if (super_site_var_index_) {
        for (uint16_t ai = 0; ai < ss.var_count; ++ai) {
            const int member_idx = (*super_site_var_index_)[ss.var_start + ai];
            if (member_idx == abs_locus) {
                view.anchor_class = static_cast<uint8_t>(ai + 1);
                break;
            }
        }
    }

    return true;
}

inline void SupersiteEmissionAdapter::build_match_mask(const SiteView& view,
                                                        unsigned int n_cond_haps,
                                                        bool use_anchor_split_semantics,
                                                        MatchMask& mask) const {
    const std::size_t total_entries = static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER;
    mask.resize(total_entries);

    if (view.emit_kind == EmitKind::Mis || !view.supersite) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), MatchMask::kMatch);
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), true);
        return;
    }

    if (!panel_codes_ || !cond_idx_) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), MatchMask::kMismatch);
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);
        return;
    }

    const SuperSite& ss = *view.supersite;
    // Parity mode (anchor-split) requires valid pointers
    if (use_anchor_split_semantics) {
        assert(panel_codes_ != nullptr && "panel_codes must be set when using anchor-split semantics");
        assert(cond_idx_ != nullptr && "cond_idx must be set when using anchor-split semantics");
    }
    if (!supersite_debug::validate_panel_span(ss, panel_codes_size_, view.supersite_index, "SupersiteEmissionAdapter::build_match_mask")) {
        std::abort();
    }
    std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), MatchMask::kMismatch);
    std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);

    const bool use_split = use_anchor_split_semantics && view.anchor_class != 0;
    // Optional trace
    const char* tr = std::getenv("SHAPEIT5_TEST_TRACE");
    if (tr && tr[0] != '\0' && tr[0] != '0') {
        std::fprintf(stdout,
                     "build_match_mask locus=%d ss_idx=%d use_split=%d anchor_class=%u emit=%d amb_mask=%u n_cond_haps=%u\n",
                     view.locus, view.supersite_index, (int)use_split, (unsigned)view.anchor_class,
                     (int)view.emit_kind, (unsigned)view.amb_mask, n_cond_haps);
    }

    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const unsigned int hap_idx = (*cond_idx_)[k];
        if (!supersite_debug::validate_panel_byte(ss, hap_idx, view.supersite_index, "SupersiteEmissionAdapter::build_match_mask")) {
            std::abort();
        }
        const uint8_t donor_code = unpackSuperSiteCode(panel_codes_, ss.panel_offset, hap_idx);
        const std::size_t base = static_cast<std::size_t>(k) * HAP_NUMBER;
        if (tr && tr[0] != '\0' && tr[0] != '0') {
            std::fprintf(stdout, "  donor k=%u hap_idx=%u code=%u\n", k, hap_idx, (unsigned)donor_code);
        }

        if (use_split) {
            const uint8_t donor_alt_flag = (donor_code == view.anchor_class) ? 1u : 0u;
            // Determine expected flag per lane: for AMB, use amb_mask; for HOM, use constant from sample_class0
            const bool is_amb = (view.emit_kind == EmitKind::Amb);
            const uint8_t hom_expected = (view.sample_class0 == view.anchor_class) ? 1u : 0u;
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const uint8_t expected_flag = is_amb ? (((view.amb_mask >> h) & 1U) ? 1u : 0u)
                                                    : hom_expected;
                const bool match = (donor_alt_flag == expected_flag);
                mask.by_donor_lane[base + h] = match ? MatchMask::kMatch : MatchMask::kMismatch;
                mask.any_match_lane[h] = mask.any_match_lane[h] || match;
            }
        } else {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const bool match = (donor_code == view.lane_class[h]);
                mask.by_donor_lane[base + h] = match ? MatchMask::kMatch : MatchMask::kMismatch;
                mask.any_match_lane[h] = mask.any_match_lane[h] || match;
            }
        }
    }
}

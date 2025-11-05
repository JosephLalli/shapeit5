#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <vector>

#include <containers/bitmatrix.h>
#include <models/site_emission_types.h>
#include <models/super_site_accessor.h>
#include <objects/genotype/genotype_header.h>
#include <utils/otools.h>

class genotype;

class BiallelicEmissionAdapter {
public:
    explicit BiallelicEmissionAdapter(const genotype* G, const bitmatrix* Hvar)
        : G_(G), Hvar_(Hvar) {}

    void build_view(int abs_locus, int curr_abs_ambiguous, SiteView& view) const;

    void build_match_mask(const SiteView& view,
                          unsigned int n_cond_haps,
                          int matrix_row,
                          MatchMask& mask) const;

private:
    const genotype* G_{nullptr};
    const bitmatrix* Hvar_{nullptr};
};

class SupersiteEmissionAdapter {
public:
    SupersiteEmissionAdapter(const genotype* G,
                              const std::vector<SuperSite>* super_sites,
                              const std::vector<int>* locus_to_super_idx,
                              const std::vector<int>* super_site_var_index,
                              const uint8_t* panel_codes,
                              const std::vector<unsigned int>* cond_idx)
        : G_(G),
          super_sites_(super_sites),
          locus_to_super_idx_(locus_to_super_idx),
          super_site_var_index_(super_site_var_index),
          panel_codes_(panel_codes),
          cond_idx_(cond_idx) {}

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
};

// =============================
// Inline implementations
// =============================

inline void BiallelicEmissionAdapter::build_view(int abs_locus,
                                                 int curr_abs_ambiguous,
                                                 SiteView& view) const {
    view.kind = SiteKind::Biallelic;
    view.supersite = nullptr;
    view.locus = abs_locus;
    view.anchor_class = 0;

    if (!G_) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        return;
    }

    const unsigned char variant = G_->Variants[DIV2(abs_locus)];
    const bool is_missing = VAR_GET_MIS(MOD2(abs_locus), variant);
    const bool is_ambiguous = VAR_GET_AMB(MOD2(abs_locus), variant);

    if (is_missing) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        return;
    }

    if (is_ambiguous) {
        view.emit_kind = EmitKind::Amb;
        uint8_t amb_code = 0u;
        if (curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) {
            amb_code = G_->Ambiguous[curr_abs_ambiguous];
        }
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool wants_alt = ((amb_code >> h) & 1U) != 0;
            view.lane_class[h] = wants_alt ? 1u : 0u;
        }
        return;
    }

    view.emit_kind = EmitKind::Hom;
    const bool hap0_is_alt = VAR_GET_HAP0(MOD2(abs_locus), variant);
    const uint8_t allele = hap0_is_alt ? 1u : 0u;
    std::fill(std::begin(view.lane_class), std::end(view.lane_class), allele);
}

inline void BiallelicEmissionAdapter::build_match_mask(const SiteView& view,
                                                       unsigned int n_cond_haps,
                                                       int matrix_row,
                                                       MatchMask& mask) const {
    const std::size_t total_entries = static_cast<std::size_t>(n_cond_haps) * HAP_NUMBER;
    mask.resize(total_entries);

    if (view.emit_kind == EmitKind::Mis) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(1));
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), true);
        return;
    }

    if (!Hvar_) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(0));
        return;
    }

    std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(0));
    std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);

    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const bool donor_alt = Hvar_->get(matrix_row, k);
        const uint8_t donor_code = donor_alt ? 1u : 0u;
        const std::size_t base = static_cast<std::size_t>(k) * HAP_NUMBER;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool match = (donor_code == view.lane_class[h]);
            mask.by_donor_lane[base + h] = static_cast<uint8_t>(match);
            mask.any_match_lane[h] = mask.any_match_lane[h] || match;
        }
    }
}

inline bool SupersiteEmissionAdapter::build_view(int abs_locus,
                                                 int curr_abs_ambiguous,
                                                 SiteView& view) const {
    view.locus = abs_locus;
    view.anchor_class = 0;
    std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);

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
    const bool is_anchor = (abs_locus == static_cast<int>(ss.global_site_id));
    view.kind = is_anchor ? SiteKind::SuperAnchor : SiteKind::SuperSibling;

    if (!is_anchor) {
        view.emit_kind = EmitKind::Mis;
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
            break;
        case SSClass::HOM:
            view.emit_kind = EmitKind::Hom;
            std::fill(std::begin(view.lane_class), std::end(view.lane_class), c0);
            break;
        case SSClass::AMB: {
            view.emit_kind = EmitKind::Amb;
            uint8_t amb_mask = 0u;
            if (G_ && curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) {
                amb_mask = G_->Ambiguous[curr_abs_ambiguous];
            }
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const bool choose_c1 = ((amb_mask >> h) & 1U) != 0;
                view.lane_class[h] = choose_c1 ? c1 : c0;
            }
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
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(1));
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), true);
        return;
    }

    if (!panel_codes_ || !cond_idx_) {
        std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(0));
        std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);
        return;
    }

    const SuperSite& ss = *view.supersite;
    std::fill(mask.by_donor_lane.begin(), mask.by_donor_lane.end(), static_cast<uint8_t>(0));
    std::fill(std::begin(mask.any_match_lane), std::end(mask.any_match_lane), false);

    const bool use_split = use_anchor_split_semantics && view.anchor_class != 0;

    for (unsigned int k = 0; k < n_cond_haps; ++k) {
        const unsigned int hap_idx = (*cond_idx_)[k];
        const uint8_t donor_code = unpackSuperSiteCode(panel_codes_, ss.panel_offset, hap_idx);
        const std::size_t base = static_cast<std::size_t>(k) * HAP_NUMBER;

        if (use_split) {
            const uint8_t donor_alt_flag = (donor_code == view.anchor_class) ? 1u : 0u;
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const uint8_t expected_flag = (view.lane_class[h] == view.anchor_class) ? 1u : 0u;
                const bool match = (donor_alt_flag == expected_flag);
                mask.by_donor_lane[base + h] = static_cast<uint8_t>(match);
                mask.any_match_lane[h] = mask.any_match_lane[h] || match;
            }
        } else {
            for (int h = 0; h < HAP_NUMBER; ++h) {
                const bool match = (donor_code == view.lane_class[h]);
                mask.by_donor_lane[base + h] = static_cast<uint8_t>(match);
                mask.any_match_lane[h] = mask.any_match_lane[h] || match;
            }
        }
    }
}

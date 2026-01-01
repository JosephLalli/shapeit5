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

private:
    const genotype* G_{nullptr};
    bitmatrix* Hvar_{nullptr};
};

class SupersiteEmissionAdapter {
public:
    SupersiteEmissionAdapter(const genotype* G,
                              const std::vector<SuperSite>* super_sites,
                              const std::vector<int>* locus_to_super_idx,
                              const std::vector<int>* super_site_var_index)
        : G_(G),
          super_sites_(super_sites),
          locus_to_super_idx_(locus_to_super_idx),
          super_site_var_index_(super_site_var_index) {}

    bool build_view(int abs_locus,
                    int curr_abs_ambiguous,
                    SiteView& view) const;

private:
    const genotype* G_{nullptr};
    const std::vector<SuperSite>* super_sites_{nullptr};
    const std::vector<int>* locus_to_super_idx_{nullptr};
    const std::vector<int>* super_site_var_index_{nullptr};
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
    view.sample_class0 = 0;
    view.sample_class1 = 0;

    if (!G_) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        view.sample_class0 = 0u;
        view.sample_class1 = 0u;
        view.amb_mask = 0u;
        return;
    }

    const unsigned char variant = G_->Variants[DIV2(abs_locus)];
    const bool is_missing = VAR_GET_MIS(MOD2(abs_locus), variant);
    const bool is_ambiguous = VAR_GET_AMB(MOD2(abs_locus), variant);

    if (is_missing) {
        view.emit_kind = EmitKind::Mis;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0);
        view.sample_class0 = 0u;
        view.sample_class1 = 0u;
        view.amb_mask = 0u;
        return;
    }

    if (is_ambiguous) {
        view.emit_kind = EmitKind::Amb;
        assert((curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) && "Ambiguous index out of range - site_emission_adapter.h, BiallelicEmissionAdapter::build_view");
        uint8_t amb_code = G_->Ambiguous[curr_abs_ambiguous];
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

inline bool SupersiteEmissionAdapter::build_view(int abs_locus,
                                                 int curr_abs_ambiguous,
                                                 SiteView& view) const {
    view.locus = abs_locus;
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

    if (!is_anchor) {
        view.emit_kind = EmitKind::Mis;
        view.sample_class0 = 0u;
        view.sample_class1 = 0u;
        view.amb_mask = 0u;
        return true;
    }

    uint8_t c0 = 0u;
    uint8_t c1 = 0u;
    if (G_) {
        // Use immutable observed genotype (c0,c1) for emissions
        G_->getSupersiteObservedGt(ss_idx, c0, c1);
    }

    const bool anchor_missing = G_ ? G_->supersiteIsMissing(ss_idx) : true;
    const bool anchor_amb = (!anchor_missing && (c0 != c1));

    if (anchor_missing) {
        // Missing anchor: emissions uninformative, lanes uniform (mirror bial MIS)
        view.emit_kind = EmitKind::Mis;
        view.sample_class0 = 0u;
        view.sample_class1 = 0u;
        view.amb_mask = 0u;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), 0u);
        return true;
    }

    if (anchor_amb) {
        view.emit_kind = EmitKind::Amb;
        uint8_t amb_mask = 0u;
        if (G_) {
            assert((curr_abs_ambiguous >= 0 && curr_abs_ambiguous < static_cast<int>(G_->Ambiguous.size())) && "Ambiguous index out of range - site_emission_adapter.h, SupersiteEmissionAdapter::build_view");
            amb_mask = G_->Ambiguous[curr_abs_ambiguous];
        }
        view.amb_mask = amb_mask;
        for (int h = 0; h < HAP_NUMBER; ++h) {
            const bool choose_c1 = ((amb_mask >> h) & 1U) != 0;
            view.lane_class[h] = choose_c1 ? c1 : c0;
        }
        view.sample_class0 = c0;
        view.sample_class1 = c1;
    } else {
        view.emit_kind = EmitKind::Hom;
        view.sample_class0 = c0;
        view.sample_class1 = c0;
        view.amb_mask = 0u;
        std::fill(std::begin(view.lane_class), std::end(view.lane_class), c0);
    }

    return true;
}

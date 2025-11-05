#pragma once

#include <cstdint>

#include <utils/otools.h> // aligned_vector32
#include <models/super_site_accessor.h> // SuperSite

enum class SiteKind : uint8_t {
    Biallelic = 0,
    SuperAnchor = 1,
    SuperSibling = 2,
};

enum class EmitKind : uint8_t {
    Hom = 0,
    Amb = 1,
    Mis = 2,
};

struct SiteView {
    SiteKind kind{SiteKind::Biallelic};
    int locus{-1};
    const SuperSite* supersite{nullptr};
    EmitKind emit_kind{EmitKind::Hom};
    uint8_t lane_class[HAP_NUMBER] = {0};
    uint8_t anchor_class{0}; // Meaningful for supersite anchors (ALT code 0..15)
};

struct MatchMask {
    aligned_vector32<uint8_t> by_donor_lane;
    bool any_match_lane[HAP_NUMBER] = {false};

    MatchMask() = default;
    explicit MatchMask(std::size_t total_entries)
        : by_donor_lane(total_entries, 0) {}

    void resize(std::size_t total_entries) {
        by_donor_lane.assign(total_entries, 0);
        std::fill(std::begin(any_match_lane), std::end(any_match_lane), false);
    }
};

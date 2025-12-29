#pragma once

#include <cstdint>

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
    int supersite_index{-1};
    const SuperSite* supersite{nullptr};
    EmitKind emit_kind{EmitKind::Hom};
    uint8_t lane_class[HAP_NUMBER] = {0};
    // Immutable supersite classes for emissions (c0/c1). These are not the sampled h0/h1.
    uint8_t sample_class0{0};
    uint8_t sample_class1{0};
    // New: biallelic ambiguous orientation bits for this site (Ambiguous mask)
    // Used to unify supersite anchor emissions with biallelic split semantics
    uint8_t amb_mask{0};
};

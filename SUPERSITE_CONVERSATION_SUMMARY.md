# Supersite Implementation - Conversation Summary

## Overview
This document summarizes a comprehensive technical discussion about the `super-sites-support` branch of the SHAPEIT5 project, which implements multiallelic variant phasing by treating split biallelic records as single "supersites" within the forward/backward HMM algorithm.

## Key Architectural Insights

### What Supersites Are
- **Problem**: When multiallelic sites (e.g., STRs with alleles 3,4,5,6) are split into multiple biallelic records for BCF compatibility, vanilla SHAPEIT5 phases each split independently, potentially creating biologically impossible haplotypes (e.g., both allele 4 and 6 on the same chromosome).
- **Solution**: Supersites group all split records at the same genomic position (`chr:bp`) and treat them as a single HMM locus, ensuring consistent phasing across all alleles.

### How Supersites Integrate with the HMM

#### Preprocessing (`buildSuperSites()`)
Called once per iteration after PBWT conditioning haplotype selection:
1. **Groups variants** by `(chr, bp)` to identify split multiallelic sites
2. **Assigns 4-bit allele codes** (0=REF, 1..15=ALT index) to each conditioning haplotype
   - Uses "first ALT seen wins" rule to canonicalize alleles
   - Packs 2 codes per byte for cache efficiency
3. **Creates SuperSite structs** with:
   - `global_site_id`: anchor variant index (where DP runs)
   - `panel_offset`: byte offset into packed codes
   - `var_start/var_count`: indices into flat member list
4. **Builds lookup tables**:
   - `locus_to_super_idx[v]`: maps variant→supersite (-1 if not in supersite)
   - `super_site_var_index`: flattened list of member variant indices

#### Forward Pass Integration
At each locus, the forward algorithm:

1. **Checks** if `locus_to_super_idx[curr_abs_locus]` indicates a supersite
2. **Classifies** the supersite genotype for this sample:
   - **MIS**: Both haplotypes missing → use MIS emission (neutral)
   - **HOM**: Both haplotypes same allele (c0==c1) → penalize mismatching donors
   - **AMB**: Two different alleles (c0≠c1) → allow both orientations, sum emissions
3. **Unpacks** conditioning haplotype codes once per supersite
4. **Computes emissions** by comparing each donor's 4-bit code to sample's allele(s)
5. **Applies** transition math (nt/yt) and emission exactly like biallelic sites
6. **Gates** DP to run only at `ss.global_site_id` (anchor); sibling splits are skipped

#### Backward Pass Integration
- Uses **same recursion** as biallelic: sample donor haplotypes via forward probs + nt/yt
- **Decodes** sampled donor haplotype's allele at supersite via 4-bit code
- **Projects** consistent phase back to each split record (e.g., hapA=4, hapB=6 → "4 vs ref" gets 1|0, "6 vs ref" gets 0|1)

### Key Design Principles
1. **Minimal divergence**: Supersite logic reuses existing HMM state arrays, AVX2 patterns, and transition math
2. **O(1) allele lookup**: Precomputed 4-bit codes avoid per-locus scanning of split records
3. **Anchor gating**: DP runs once per biological site, not once per split, preventing double-counting
4. **Classification correctness**: Supersite genotype classification (HOM/AMB/MIS) drives emission, not per-split flags

## Critical Bugs Identified

### 1. INIT_AMB Overlapping Store (Data Corruption)
**Problem**: Two 256-bit stores at `&prob[i]` and `&prob[i+4]` write 16 floats, corrupting the next block.
```cpp
_mm256_store_ps(&prob[i], _prob0);    // 8 floats at i
_mm256_store_ps(&prob[i+4], _prob1);  // 8 floats at i+4 = 16 total
```
**Fix**: Build single `_combined` vector from 128-bit halves, store once (RUN_AMB pattern).

### 2. AMB Lane Mask Broken
**Problem**: Supersite AMB splits emissions by low/high 128-bit halves (lanes 0-3 vs 4-7), breaking the original lane semantics where `amb_code` determines which lanes want hap0 vs hap1.
```cpp
// WRONG: assumes lanes 0-3=hap0, 4-7=hap1
_prob_low  = _mm_mul_ps(_prob_low,  _emit0);
_prob_hi   = _mm_mul_ps(_prob_hi,   _emit1);
```
**Fix**: Build `want_c0_mask` and `want_c1_mask` from `amb_code` bits, blend scalar emissions per lane.

### 3. Per-Split Classification
**Problem**: Uses biallelic `amb/hom/mis` flags to route supersite logic, misclassifying cases like 5|6 (looks HOM at 4-split anchor).
**Fix**: Centralize classification via `getSampleSuperSiteAlleleCode` for both haplotypes when `ss_idx>=0`.

### 4. Double Counting at Siblings
**Problem**: DP runs at every supersite member variant, not just the anchor.
**Fix**: Gate `if (curr_abs_locus != ss.global_site_id) return;` before emission/transition in all SS_* functions.

### 5. Missing HOM Fallback
**Problem**: AMB path doesn't route to HOM when `c0==c1`.
**Fix**: Add equality check in `classify_supersite`, route to `SS_INIT/RUN/COLLAPSE_HOM`.

### 6. Repeated Unpacking (Performance)
**Problem**: Decodes panel codes every INIT/RUN/COLLAPSE call.
**Fix**: Cache `ss_cond_codes` per supersite index, invalidate on `H.select()`. Saves 10-20% in supersite-heavy windows.

## Coding Pattern Alignment Recommendations

### Current Problem
Supersite logic is embedded in biallelic functions via `if (ss_idx>=0)` branches, diverging from the clean taxonomy of the original code.

### Target Pattern
Mirror the biallelic structure exactly:

#### Function Taxonomy
Create separate functions for each genotype type:
```cpp
// Separate functions matching biallelic signatures
inline void SS_INIT_HOM();
inline void SS_INIT_AMB();
inline void SS_INIT_MIS();

inline bool SS_RUN_HOM(char rare_allele);
inline void SS_RUN_AMB();
inline void SS_RUN_MIS();

inline void SS_COLLAPSE_HOM();
inline void SS_COLLAPSE_AMB();
inline void SS_COLLAPSE_MIS();
```

#### Dispatch Pattern
At the top of each existing function:
```cpp
inline void INIT_HOM() {
    int ss_idx = (super_sites && locus_to_super_idx) 
                 ? (*locus_to_super_idx)[curr_abs_locus] : -1;
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        // Anchor gate
        if (curr_abs_locus != (int)ss.global_site_id) return;
        
        // Classify supersite
        uint8_t c0, c1;
        auto cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
        
        // Route to appropriate SS_* function
        switch (cls) {
            case SSClass::HOM: return SS_INIT_HOM();
            case SSClass::AMB: return SS_INIT_AMB();
            case SSClass::MIS: return SS_INIT_MIS();
        }
    }
    // Original biallelic logic continues...
}
```

#### Lane Mask Semantics
For AMB paths, preserve the lane mask from `amb_code`:
```cpp
// Build masks from amb_code (which lanes want c0 vs c1)
__m256 want_c0_mask = ...; // from ~amb_code bits
__m256 want_c1_mask = ...; // from amb_code bits

for each donor k:
    float e0 = (ss_cond_codes[k] == c0) ? 1.0f : mismatch;
    float e1 = (ss_cond_codes[k] == c1) ? 1.0f : mismatch;
    
    __m256 _v0 = _mm256_set1_ps(e0);
    __m256 _v1 = _mm256_set1_ps(e1);
    
    // Blend per lane based on which allele that lane wants
    __m256 _emit = _mm256_or_ps(
        _mm256_and_ps(want_c0_mask, _v0),
        _mm256_and_ps(want_c1_mask, _v1)
    );
    
    // Standard HMM update
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
    _prob = _mm256_mul_ps(_prob, _emit);
```

#### Naming Consistency
Use **identical** local variable names as biallelic code:
- AVX registers: `_sum`, `_prob`, `_emit`, `_factor`, `_tFreq`, `_nt`, `_mismatch`
- Class members: `probSumH`, `probSumT`, `probSumK` (no new accumulators)
- Supersite-specific: `ss_emissions`, `ss_emissions_h1`, `ss_cond_codes`

#### Caching Pattern
```cpp
// Class members
std::vector<bool> ss_cached;
aligned_vector32<uint8_t> ss_cond_codes;

// Helper
inline void ss_load_cond_codes(const SuperSite& ss, int ss_idx) {
    if (!ss_cached[ss_idx]) {
        for (int k=0; k<(int)n_cond_haps; k++) {
            unsigned int gh = (*cond_idx)[k];
            ss_cond_codes[k] = unpackSuperSiteCode(panel_codes, ss.panel_offset, gh);
        }
        ss_cached[ss_idx] = true;
    }
}
```

## Benefits of Alignment

### Maintainability
- **Easier review**: Diff SS_INIT_HOM vs INIT_HOM line-by-line
- **Clear separation**: Supersite logic not scattered through biallelic branches
- **Familiar structure**: New contributors recognize the INIT/RUN/COLLAPSE pattern

### Correctness
- **Forced parity**: Identical structure prevents subtle divergences (e.g., missing normalization steps)
- **Type safety**: Inline functions catch signature mismatches that macros hide
- **Testability**: Can unit-test SS_* functions independently with synthetic inputs

### Performance
- **Cache-friendly**: Same memory access patterns, same vectorization
- **No overhead when disabled**: Supersite pointers null → branches never taken
- **Explicit caching**: `ss_load_cond_codes` helper makes optimization obvious

## Memory Layout and SIMD Constraints

### 8-Lane Blocks (HAP_NUMBER=8)
- Each block is exactly 8 floats (256 bits)
- **One** `_mm256_store_ps` per block
- Overlapping stores corrupt adjacent state

### Aligned Vectors
- All probability arrays use `aligned_vector32<T>` (32-byte aligned)
- Required for `_mm256_load_ps` / `_mm256_store_ps`
- Scratch buffers must be class members, not local temporaries (avoid heap churn)

### Lane Semantics
- `amb_code` is an 8-bit mask where bit `h` indicates lane `h`'s orientation
- Biallelic AMB: `g0[h]` = emission if lane wants hap0, `g1[h]` if wants hap1
- Supersite AMB: must preserve this → blend two scalars per lane, not split halves

## Next Actions

See the newly created todo lists:
1. **Bug Fixes**: 8 items addressing overlapping store, lane masks, classification, anchor gating, caching
2. **Refactoring**: 13 items for extracting SS_* functions, unifying naming, implementing helpers, adding tests

All items are documented with specific file locations, function names, and code patterns to follow.

## References
- Original conversation: `untitled:Untitled-1`
- Updated documentation: `.github/copilot-instructions.md`
- Bug todo list: Created in workspace
- Refactoring todo list: Created in workspace

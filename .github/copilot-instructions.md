# SHAPEIT5 AI Coding Agent Instructions

## Project Overview
SHAPEIT5 is a high-performance C++17 haplotype phasing toolkit for genomic data (WGS/array), optimized with AVX2 SIMD and processing BCF/VCF via HTSlib. The codebase includes:
- **phase_common**: Phase common variants (SNP arrays)
- **phase_rare**: Phase rare variants onto common scaffolds using Li & Stephens HMM
- **ligate/switch/simulate/xcftools**: Supporting tools
- **tests/**: Unit test harness for supersite integration (experimental multi-allelic phasing)

Current development focuses on optional **supersite phasing** for multi-allelic loci (`--enable-supersites`), maintaining backward compatibility with default biallelic paths.

## Architecture & Data Flow

### HMM Core (`phase_common/src/models/`)
- **`haplotype_segment_single`** / **`haplotype_segment_double`**: Forward-backward HMM implementations (float/double precision)
  - Single precision (float) is default for performance; double precision used on underflow
  - Each uses `aligned_vector32<T>` (32-byte aligned via Boost) for SIMD-friendly Alpha/Beta arrays
  - **INIT/RUN/COLLAPSE** state transitions encoded as inlined functions (one set per genotype type: HOM/AMB/MIS)
  - Operates on 8-lane vectors (HAP_NUMBER=8) using AVX2 intrinsics
  - Supersite support gated by `if (super_sites && locus_to_super_idx)` checks – avoids overhead when disabled
- **Transition model**: Li & Stephens HMM with recombination distance-based state switching
  - `nt` = probability of no transition (stay on same donor haplotype)
  - `yt` = probability of recombination (switch to new donor, spread across K states)
  - Transition: `trans[h] = nt * prob_prev[h] + yt * AlphaSum_prev`
- **Emission model**: per-donor-haplotype match/mismatch probabilities
  - HOM: ~1.0 if donor matches, ~ε (M.ed/M.ee) if mismatches
  - AMB: uses lane mask (amb_code) to allow both 0|1 and 1|0 orientations
  - MIS: ~1.0 for all donors (uninformative, imputed via flanking sites)

### Supersite Data Model (`phase_common/src/objects/super_site_builder.cpp`)
- **Purpose**: Treat multiallelic loci (split into multiple biallelic records) as single HMM positions
- **`buildSuperSites()`** preprocessing (called once per iteration after PBWT selection):
  - Groups split multi-allelic records by `chr:bp`
  - Assigns 4-bit allele codes per conditioning haplotype: 0=REF, 1..15=which ALT (first match wins)
  - Packs codes densely: 2 codes per byte in `packed_allele_codes` (uint8_t array)
  - Chunks large multiallelic sites (>SUPERSITE_MAX_ALTS=15) into multiple SuperSite records
- **SuperSite struct fields**:
  - `global_site_id`: anchor variant index (first split record) – the locus where DP actually runs
  - `panel_offset`: byte offset into `packed_allele_codes` for this supersite's haplotype codes
  - `var_start`, `var_count`: indices into `super_site_var_index` (flat list of member variant IDs)
  - `n_alts`: number of alternate alleles in this chunk
- **Lookup structures** (all sized to total variant count):
  - `locus_to_super_idx[v]`: maps variant index v → supersite index s (or -1 if not in supersite)
  - `super_site_var_index`: flattened array of variant indices (CSR-style, indexed by SuperSite.var_start/count)
  - `is_super_site[v]`: legacy boolean flag (superseded by locus_to_super_idx)
- **Runtime allele lookup**:
  - Panel: `unpackSuperSiteCode(panel_codes, ss.panel_offset, hap_idx)` → O(1) 4-bit code
  - Sample: `getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap)` → infer from split records

### Emission Kernels
- **Biallelic**: direct AVX2 ops on conditioning haplotypes (bitmatrix Hvar access)
  - HOM: broadcast match/mismatch scalar per donor
  - AMB: precompute g0/g1 emission vectors using amb_code lane mask, blend with donor allele
  - MIS: no emission penalty, all donors get weight ~1.0
- **Supersite**: multiallelic emission via 4-bit code comparison
  - Precompute per-donor emissions using `precomputeSuperSiteEmissions_FloatScalar`
  - HOM: scalar match if `donor_code == sample_code`, else mismatch
  - AMB: **must preserve lane mask semantics** – blend emissions for c0/c1 per lane using amb_code
    - Low/high half split breaks original model; use bitwise masks like biallelic AMB
  - MIS: fallback to standard MIS emission (all ~1.0)
  - **Scratch buffers**: `ss_emissions`, `ss_emissions_h1` must be persistent `aligned_vector32` members to avoid heap churn in backward pass
  - **Classification**: per-supersite, per-sample → HOM (c0==c1), AMB (c0≠c1), or MIS (both missing)
- **Anchor gating**: only run DP at `ss.global_site_id`; sibling splits skip emission/transition to avoid double counting

### Threading & Precision Fallback
- PBWT conditioning selection → window segmentation → multi-threaded HMM
- If float HMM underflows (`outcome != 0`), retry with `haplotype_segment_double`
- Statistics tracked: `statH` (K states), `statS` (window size), underflow counts

## Build System

### Standard Build
```bash
make -j$(nproc)                    # Top-level: phase_common, phase_rare, switch, ligate, simulate, xcftools, tests
make -C phase_common               # Single module
make -C tests                      # Unit tests only
```

### Targets & Environment
- **`desktop`** (default): assumes `~/.linuxbrew` for htslib/boost
- **`static_exe`**: produces `bin/<module>_static` linked with `.a` libs
- **`debug`**: `-g` flags, different htslib path (see `makefile_common.mk`)

### Required Dependencies
- **htslib** 1.12+ (VCF/BCF I/O, static linking for exe)
- **Boost**: `program_options`, `iostreams`, `serialization` (dev + static libs)
- **Compiler**: g++ with C++17, AVX2 support (`-mavx2 -mfma`)

**Critical**: Tests require `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib` at runtime for dynamic boost libs.

### Makefile Structure
- **`common/makefile_common.mk`**: Shared flags, library discovery, git commit tracing (`__COMMIT_ID__`)
- Each module includes this via `-include ../common/makefile_common.mk`
- **`tests/makefile`**: Links against `phase_common` objects (`haplotype_segment_single.o`, `super_site_builder.o`, etc.)

## Testing Conventions

### Unit Tests (`tests/`)
```bash
make -C tests                      # Build all test binaries
make -C tests test-run             # Build + run all (sets LD_LIBRARY_PATH)
LD_LIBRARY_PATH=$HOME/.linuxbrew/lib tests/bin/test_supersite_hmm_states
```

- **Test harnesses** use `#define private public` to expose internals for state validation
- **`test_toolbox.cpp`** defines global objects (`rng`, `vrb`, etc.) via `_DECLARE_TOOLBOX_HERE` (instantiate exactly once)
- Synthetic data: no external fixtures yet (placeholders under `tests/data/`)

### Integration Tests (`test/scripts/`)
```bash
bash test/scripts/phase.wgs.unrelated.sh
```
- Compare phased BCF output against `expected/` golden files using `bcftools view -G` diffs
- Require `bcftools` in PATH; tests clean up temp dirs via `trap`

### CI (`.github/workflows/build.yml`)
- Ubuntu 20.04, installs dependencies from apt + wget htslib source
- Runs `make -j$(nproc)` then `test/scripts/phase*.sh` regression suite

## Code Patterns & Conventions

### Memory Alignment & SIMD
- **Always use** `aligned_vector32<T>` for arrays consumed by AVX2 intrinsics
- Avoid `std::vector<float>` in hot paths (causes allocator churn); hoist to class members
- Scratch buffers must be pre-allocated: `aligned_vector32<float> ss_emissions;` as class member, not local temporary
- Each block operates on HAP_NUMBER=8 lanes (256 bits = 8×float or 4×double)
- One `_mm256_store_ps` per block; overlapping stores corrupt adjacent state

### Supersite Coding Patterns (CRITICAL for Maintainability)
**Goal**: Supersite code must exactly mirror biallelic code structure to minimize divergence and facilitate review/testing.

#### Function Taxonomy
- Create separate `SS_INIT_HOM()`, `SS_INIT_AMB()`, `SS_INIT_MIS()` functions (not inline branches)
- Mirror signatures and return types: `SS_RUN_HOM(char rare_allele)` returns `bool` like `RUN_HOM`
- Dispatch pattern in existing functions:
  ```cpp
  inline void INIT_HOM() {
      int ss_idx = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[curr_abs_locus] : -1;
      if (ss_idx >= 0) {
          const SuperSite& ss = (*super_sites)[ss_idx];
          // Anchor gate: only run DP at global_site_id
          if (curr_abs_locus != (int)ss.global_site_id) return; // sibling: no-op
          
          uint8_t c0, c1;
          auto cls = classify_supersite(G, ss, *super_site_var_index, c0, c1);
          switch (cls) {
            case SSClass::HOM: return SS_INIT_HOM();
            case SSClass::AMB: return SS_INIT_AMB();
            case SSClass::MIS: return SS_INIT_MIS();
          }
      }
      // ... original biallelic logic
  }
  ```

#### Classification Logic
- **Centralize** per-supersite classification; do NOT use per-split amb/hom/mis flags
  ```cpp
  enum class SSClass { MIS, HOM, AMB };
  
  inline SSClass classify_supersite(const genotype* G, const SuperSite& ss,
                                     const std::vector<int>& idx,
                                     uint8_t& c0, uint8_t& c1) {
      c0 = getSampleSuperSiteAlleleCode(G, ss, idx, 0);
      c1 = getSampleSuperSiteAlleleCode(G, ss, idx, 1);
      if (c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING) return SSClass::MIS;
      if (c0 == c1) return SSClass::HOM;
      return SSClass::AMB;
  }
  ```

#### Lane Mask Semantics (AMB path)
- **Never split lanes by low/high halves** for different emissions
- **Preserve amb_code lane semantics**:
  - Biallelic: each lane `h` "wants" hap0 or hap1 allele based on `HAP_GET(amb_code, h)`
  - Supersite: build two masks (want_c0, want_c1) from amb_code, blend scalar emissions per lane
  ```cpp
  // Pseudo-code pattern for SS_RUN_AMB:
  for each donor k:
      float e0 = (ss_cond_codes[k] == c0) ? 1.0f : mismatch;
      float e1 = (ss_cond_codes[k] == c1) ? 1.0f : mismatch;
      __m256 _v0 = _mm256_set1_ps(e0);
      __m256 _v1 = _mm256_set1_ps(e1);
      __m256 _emit = _mm256_or_ps(_mm256_and_ps(want_c0_mask, _v0),
                                  _mm256_and_ps(want_c1_mask, _v1));
      _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq);
      _prob = _mm256_mul_ps(_prob, _emit);
  ```
- **Do not** use separate 128-bit low/high emissions (`_emit0`, `_emit1`) unless it's purely for packing, not semantics

#### Anchor Gating
- Only execute DP (transition + emission) at `ss.global_site_id` (anchor)
- Sibling splits: skip DP but maintain bookkeeping (`AlphaMissing`, `curr_abs_missing++`, etc.)
  ```cpp
  if (ss_idx >= 0 && curr_abs_locus != (int)ss.global_site_id) {
      // This is a sibling: no DP update, but still track missing if needed
      if (mis) { AlphaMissing[curr_rel_missing] = prob; ... }
      return (appropriate default);
  }
  ```

#### Caching & Performance
- Cache `ss_cond_codes` per supersite (don't recompute per INIT/RUN/COLLAPSE)
  - Use `std::vector<bool> ss_cached` and `inline void ss_load_cond_codes(ss)` helper
- Precompute emission vectors once, not per function
- Reuse same AVX broadcast pattern as biallelic: `_factor`, `_tFreq`, `_nt`, `_mismatch`

#### Naming Consistency
- Local AVX vars: `_sum`, `_prob`, `_emit`, `_factor`, `_tFreq`, `_nt`, `_mismatch` (identical to biallelic)
- Class members: `ss_emissions`, `ss_emissions_h1`, `ss_cond_codes` (all `aligned_vector32`)
- No new accumulators unless strictly necessary; reuse `probSumH`, `probSumT`, `probSumK`

### Error Handling & Underflow
- HMM `forward()`/`backward()` return `int`: 0=success, non-zero=underflow
- Caller (`phaser_algorithm.cpp`) retries with double precision on failure
- **No exceptions** in hot paths; use asserts for invariant checks

### Global Singletons (`_DECLARE_TOOLBOX_HERE`)
- **Exactly one** `.cpp` defines this macro to instantiate `rng`, `stb`, `alg`, `vrb`, `tac`
- All other files declare as `extern` (via `common/src/utils/otools.h`)
- In tests: `test_toolbox.cpp` owns the declaration

### Commit Tracing
- Makefile injects `__COMMIT_ID__` and `__COMMIT_DATE__` via `-D` flags (from `git rev-parse`)
- Displayed in binary `--version` output

## Debugging Workflows

### Valgrind / GDB
```bash
make -C phase_common debug         # -g flags, no optimization
gdb --args phase_common/bin/phase_common --input test/wgs/target.unrelated.bcf ...
valgrind --leak-check=full tests/bin/test_supersite_hmm
```

### Logging
- Use `vrb.bullet()`, `vrb.title()` for structured console output
- HMM state dumps: inspect `Alpha`, `AlphaSum` members via debugger (arrays aligned to 32 bytes)

### Common Pitfalls
1. **Heap allocations in RUN_AMB**: Temporary `std::vector<float>` inside loop causes `posix_memalign` crashes → hoist to member
2. **Missing LD_LIBRARY_PATH**: Tests segfault if boost libs not found at runtime
3. **Uninitialized supersite pointers**: Always gate with null checks (`if (super_sites)`)
4. **AVX2 alignment**: Segfaults if loading from non-32-byte-aligned address with `_mm256_load_ps`

## Current Sprint (AGENTS.md Context)

### Known Bugs (High Priority)
1. **INIT_AMB supersite overlapping store**: Two 256-bit stores at `&prob[i]` and `&prob[i+4]` write 16 floats, corrupting next block
   - Fix: Build single `_combined` vector from low/high 128-bit halves (like RUN_AMB does)
2. **AMB supersite lane mask broken**: Currently splits emissions by low/high halves, not by amb_code semantics
   - Fix: Build per-lane want_c0/want_c1 masks from amb_code, blend scalar emissions
3. **Per-split classification drives supersite**: Uses biallelic amb/hom/mis flags instead of supersite genotype
   - Fix: Always classify via `getSampleSuperSiteAlleleCode` for both haplotypes when ss_idx>=0
4. **Double counting at sibling loci**: DP runs at every supersite member, not just anchor
   - Fix: Gate `if (curr_abs_locus != ss.global_site_id) return;` before emission/transition
5. **Missing c0==c1 HOM fallback**: AMB path doesn't route to HOM when both haplotypes carry same allele
   - Fix: Add equality check in classification, route to SS_INIT/RUN/COLLAPSE_HOM
6. **Repeated `ss_cond_codes` unpacking**: Decodes panel codes every INIT/RUN/COLLAPSE (hot path)
   - Fix: Cache per supersite index, invalidate on H.select()

### Outstanding Refactoring
1. **Function taxonomy misalignment**: Supersite logic embedded in biallelic functions via if-branches
   - Target: Separate `SS_INIT_HOM`, `SS_RUN_HOM`, etc. with identical signatures/returns
2. **Macro hygiene**: `super_site_macros.h` mixes double/float, references nonexistent G members
   - Target: Delete macros or reimplement as inline functions matching real API
3. **Naming inconsistency**: Some paths use `_emit0/_emit1` split, others use `_combined`
   - Target: Uniform naming matching biallelic conventions

### Next Steps (Priority Order)
1. **Fix INIT_AMB overlapping store** (correctness blocker)
2. **Fix AMB lane mask semantics** (correctness blocker)
3. **Implement anchor gating** (prevents double counting)
4. **Centralize classification** (handles edge cases like HOM masquerading as AMB)
5. **Extract SS_* functions** (maintainability / review)
6. **Add ss_cond_codes caching** (performance, 10–20% in SS-heavy windows)
7. **Backward pass validation**: Re-enable assertions in `test_supersite_hmm_states` once heap churn resolved
8. **Smoke test fixtures**: Populate `tests/data/` with real VCF examples for integration testing

## Key Files Reference

| Path | Purpose |
|------|---------|
| `common/src/utils/otools.h` | Global typedefs (`aligned_vector32`), toolbox externs |
| `phase_common/src/models/haplotype_segment_single.{h,cpp}` | Float HMM core + supersite integration |
| `phase_common/src/models/super_site_macros.h` | Emission macro templates (INIT/RUN/COLLAPSE) |
| `phase_common/src/objects/super_site_builder.{h,cpp}` | Multi-allelic grouping + panel encoding |
| `phase_common/src/phaser/phaser_algorithm.cpp` | Window segmentation, threading, precision fallback |
| `tests/makefile` | Unit test build rules, external object linking |
| `test/scripts/lib/test_utils.sh` | BCF comparison helpers for integration tests |

## External Resources
- **HTSlib docs**: https://www.htslib.org/doc/
- **AVX2 intrinsics**: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/
- **Li & Stephens model**: Genetics 2003, doi:10.1534/genetics.103.015354

## Todo Lists

### Bug Fixes (Priority Order)

#### 1. Fix INIT_AMB Overlapping Store
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Function**: `INIT_AMB()` supersite branch  
**Issue**: Two 256-bit stores at `&prob[i]` and `&prob[i+4]` write 16 floats total, corrupting the next 8-lane block.  
**Fix Pattern**:
```cpp
// Replace current overlapping stores:
_mm256_store_ps(&prob[i], _prob0);
_mm256_store_ps(&prob[i+4], _prob1);

// With single combined store (matching RUN_AMB pattern):
__m128 _prob_low = _mm256_castps256_ps128(_prob0);
__m128 _prob_hi  = _mm256_extractf128_ps(_prob1, 0);
__m256 _combined = _mm256_castps128_ps256(_prob_low);
_combined = _mm256_insertf128_ps(_combined, _prob_hi, 1);
_mm256_store_ps(&prob[i], _combined);
```

#### 2. Fix AMB Lane Mask Semantics
**Files**: `phase_common/src/models/haplotype_segment_single.{h,cpp}`  
**Functions**: All `SS_INIT_AMB`, `SS_RUN_AMB`, `SS_COLLAPSE_AMB` paths  
**Issue**: Currently splits emissions by low/high 128-bit halves instead of using `amb_code` lane semantics. Each lane should blend c0/c1 emissions based on which haplotype orientation it represents.  
**Fix Pattern**:
```cpp
// Build lane masks from amb_code (each bit indicates if lane wants hap0 or hap1)
__m256 want_c0_mask = build_amb_mask_for_c0(amb_code);
__m256 want_c1_mask = build_amb_mask_for_c1(amb_code);

// Per donor k:
float e0 = (ss_cond_codes[k] == c0) ? 1.0f : mismatch;
float e1 = (ss_cond_codes[k] == c1) ? 1.0f : mismatch;
__m256 _v0 = _mm256_set1_ps(e0);
__m256 _v1 = _mm256_set1_ps(e1);
__m256 _emit = _mm256_or_ps(_mm256_and_ps(want_c0_mask, _v0),
                            _mm256_and_ps(want_c1_mask, _v1));
```

#### 3. Centralize Supersite Classification
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Location**: Add new helper function  
**Issue**: Currently uses per-split biallelic `amb/hom/mis` flags to drive supersite logic, causing misclassification (e.g., 5|6 at 4-split looks HOM).  
**Fix Pattern**:
```cpp
enum class SSClass { MIS, HOM, AMB };

inline SSClass classify_supersite(const genotype* G, const SuperSite& ss,
                                   const std::vector<int>& idx,
                                   uint8_t& c0, uint8_t& c1) {
    c0 = getSampleSuperSiteAlleleCode(G, ss, idx, 0);
    c1 = getSampleSuperSiteAlleleCode(G, ss, idx, 1);
    if (c0 == SUPERSITE_CODE_MISSING && c1 == SUPERSITE_CODE_MISSING) 
        return SSClass::MIS;
    if (c0 == c1) 
        return SSClass::HOM;
    return SSClass::AMB;
}
```
**Usage**: Replace all `if (hom/amb/mis)` checks in supersite branches with `classify_supersite()` dispatch.

#### 4. Implement Anchor Gating
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Functions**: All `INIT_*`, `RUN_*`, `COLLAPSE_*` supersite branches  
**Issue**: DP runs at every supersite member variant instead of just the anchor (`global_site_id`), causing double-counting of recombination events.  
**Fix Pattern**:
```cpp
inline void INIT_HOM() {
    int ss_idx = (super_sites && locus_to_super_idx) ? (*locus_to_super_idx)[curr_abs_locus] : -1;
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        
        // Anchor gate: only run DP at global_site_id
        if (curr_abs_locus != (int)ss.global_site_id) {
            // This is a sibling split: skip DP but maintain bookkeeping
            if (mis) { 
                AlphaMissing[curr_rel_missing] = prob; 
                curr_abs_missing++; 
                curr_rel_missing++; 
            }
            return;
        }
        
        // ... proceed with supersite DP
    }
    // ... original biallelic logic
}
```
**Apply to**: All 9 functions (INIT/RUN/COLLAPSE × HOM/AMB/MIS)

#### 5. Add HOM Fallback for c0==c1
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Location**: `classify_supersite()` helper (from item 3)  
**Issue**: AMB path doesn't route to HOM when both haplotypes carry the same alternate allele (e.g., 2|2).  
**Fix**: Already included in `classify_supersite()` implementation above (`if (c0 == c1) return SSClass::HOM`). Ensure dispatch logic routes HOM classification to `SS_INIT/RUN/COLLAPSE_HOM()`.

#### 6. Cache ss_cond_codes Per Supersite
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Members**: Add `std::vector<bool> ss_cached;` and update `ss_cond_codes` usage  
**Issue**: Unpacks panel codes repeatedly in every INIT/RUN/COLLAPSE call (hot path).  
**Fix Pattern**:
```cpp
// In class definition:
aligned_vector32<uint8_t> ss_cond_codes;
std::vector<bool> ss_cached;

// Helper function:
inline void ss_load_cond_codes(const SuperSite& ss, int ss_idx) {
    if (ss_cached[ss_idx]) return;
    
    ss_cond_codes.resize(H.n_haps);
    for (int k = 0; k < H.n_haps; k++) {
        ss_cond_codes[k] = unpackSuperSiteCode(
            packed_allele_codes->data(), 
            ss.panel_offset, 
            k
        );
    }
    ss_cached[ss_idx] = true;
}

// In INIT/RUN/COLLAPSE:
if (ss_idx >= 0) {
    const SuperSite& ss = (*super_sites)[ss_idx];
    ss_load_cond_codes(ss, ss_idx);
    // ... use ss_cond_codes[k] instead of unpacking inline
}
```
**Invalidation**: Reset `ss_cached` after `H.select()` (new conditioning haplotypes).

#### 7. Fix TRANS_DIP_ADD Typo
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Function**: Likely in transition probability computation  
**Issue**: Typo in variable name or macro.  
**Fix**: Search for `TRANS_DIP_ADD` and correct to intended name (likely `TRANS_DIP_ADP` or similar based on context).

#### 8. Verify rare_allele Handling in SS_RUN_HOM
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Function**: `SS_RUN_HOM(char rare_allele)` (when extracted)  
**Issue**: Need to verify that `rare_allele` parameter is correctly derived from supersite sample code, not from biallelic split record.  
**Verification**: Ensure `rare_allele` comes from supersite classification (c0 or c1 from `getSampleSuperSiteAlleleCode`), not from `G->al0/al1`.

---

### Refactoring Tasks

#### 1. Extract SS_INIT_HOM Function
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Create standalone `inline void SS_INIT_HOM()` matching biallelic `INIT_HOM()` structure.  
**Pattern**:
```cpp
inline void SS_INIT_HOM() {
    // Assumes supersite context already established (ss, c0, c1 available)
    // Mirror biallelic INIT_HOM logic but use ss_cond_codes instead of H.get()
    // ... implementation
}
```
**Call site**: From `INIT_HOM()` after classification dispatch.

#### 2. Extract SS_INIT_AMB Function
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Create standalone `inline void SS_INIT_AMB()`.  
**Requirements**: Fix lane mask semantics (bug #2), implement proper `_combined` store pattern (bug #1).

#### 3. Extract SS_INIT_MIS Function
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Create standalone `inline void SS_INIT_MIS()`.  
**Pattern**: Should be simplest – all donors get emission ~1.0.

#### 4. Extract SS_RUN_* Functions (HOM/AMB/MIS)
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Create `inline bool SS_RUN_HOM(char rare_allele)`, `inline bool SS_RUN_AMB()`, `inline bool SS_RUN_MIS()`.  
**Signature**: Match biallelic `RUN_*` returns (bool for segment continuation).

#### 5. Extract SS_COLLAPSE_* Functions (HOM/AMB/MIS)
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Create `inline void SS_COLLAPSE_HOM()`, `inline void SS_COLLAPSE_AMB()`, `inline void SS_COLLAPSE_MIS()`.  
**Pattern**: Handle segment boundary updates for supersite positions.

#### 6. Add Dispatcher Logic to Existing Functions
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Functions**: `INIT_HOM`, `INIT_AMB`, `INIT_MIS`, `RUN_HOM`, `RUN_AMB`, `RUN_MIS`, `COLLAPSE_HOM`, `COLLAPSE_AMB`, `COLLAPSE_MIS`  
**Pattern**: (See "Function Taxonomy" in Supersite Coding Patterns section above)  
**Goal**: Each existing function checks for supersite, applies anchor gating, classifies, then dispatches to `SS_*` variant or proceeds with biallelic logic.

#### 7. Unify Naming Conventions
**Files**: All supersite-related code in `haplotype_segment_single.{h,cpp}`  
**Goal**: Ensure all AVX variables follow biallelic naming: `_sum`, `_prob`, `_emit`, `_factor`, `_tFreq`, `_nt`, `_mismatch`.  
**Action**: Remove inconsistent names like `_emit0`, `_emit1`, `_prob0`, `_prob1` unless they're intermediate 128-bit extracts (document if kept).

#### 8. Implement Lane Mask Semantics Properly
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Location**: All `SS_*_AMB` functions  
**Goal**: Build `want_c0_mask` and `want_c1_mask` from `amb_code`, blend per-donor scalar emissions per lane.  
**Helper**: Consider adding `inline __m256 build_amb_mask_for_allele(unsigned char amb_code, int allele_idx)` utility.

#### 9. Add ss_load_cond_codes Helper
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Goal**: Centralize caching logic (bug #6).  
**Pattern**: See bug fix #6 above.

#### 10. Build AMB Lane Mask Helpers
**File**: `phase_common/src/models/haplotype_segment_single.h` or new utility header  
**Goal**: Create helpers to construct AVX masks from `amb_code` bits.  
**Example**:
```cpp
inline __m256 build_amb_mask_for_c0(unsigned char amb_code) {
    // For each lane h, set all bits to 1 if HAP_GET(amb_code, h) == 0
    // (i.e., lane wants haplotype 0 allele = c0)
    unsigned int mask_bits = 0;
    for (int h = 0; h < 8; h++) {
        if (!HAP_GET(amb_code, h)) {
            mask_bits |= (0xFF << (h * 4)); // Set 32 bits for this float
        }
    }
    return _mm256_castsi256_ps(_mm256_set1_epi32(mask_bits));
}
```
**Note**: Optimize if needed; this is illustrative.

#### 11. Clean Up super_site_macros.h
**File**: `phase_common/src/models/super_site_macros.h`  
**Goal**: Either delete (if superseded by extracted functions) or reimplement as inline functions matching real API.  
**Issue**: Current macros mix double/float, reference nonexistent `G` members.  
**Decision**: If macros still used, convert to templates or inline functions with proper type safety.

#### 12. Verify Pattern Alignment with Biallelic Code
**Files**: Compare `haplotype_segment_single.h` supersite paths with biallelic implementations  
**Goal**: Ensure every `SS_*` function mirrors its biallelic counterpart in structure, variable names, and AVX patterns.  
**Checklist**:
- [ ] Same local variable names (`_sum`, `_prob`, etc.)
- [ ] Same transition computation order (nt/yt calculation)
- [ ] Same emission broadcast pattern
- [ ] Same accumulator usage (`probSumH`, `probSumT`, `probSumK`)
- [ ] Same return semantics (bool for RUN, void for INIT/COLLAPSE)

#### 13. Ensure Bookkeeping Parity (AlphaMissing, counters)
**File**: `phase_common/src/models/haplotype_segment_single.h`  
**Functions**: All `SS_*` functions  
**Goal**: Ensure missing data tracking matches biallelic paths.  
**Pattern**: When `SSClass::MIS`, update `AlphaMissing[curr_rel_missing]`, increment `curr_abs_missing` and `curr_rel_missing`.  
**Anchor gate**: At sibling loci, still update missing counters even though DP is skipped.

#### 14. Add Comprehensive Unit Tests
**File**: `tests/src/test_supersite_hmm_states.cpp` (expand existing)  
**Goal**: Test all 9 `SS_*` functions with known inputs/outputs.  
**Coverage**:
- [ ] SS_INIT_HOM/AMB/MIS with various c0/c1 combinations
- [ ] SS_RUN_HOM/AMB/MIS with transition probabilities
- [ ] SS_COLLAPSE_HOM/AMB/MIS at segment boundaries
- [ ] Anchor gating (sibling loci produce no-op)
- [ ] Lane mask semantics (each of 8 lanes gets correct emission)
- [ ] Caching (repeated calls to same supersite use cached codes)
- [ ] Backward pass (re-enable assertions once heap churn fixed)

---

### Priority Execution Order

**Phase 1 - Critical Correctness** (blocking production use):
1. Bug #1: Fix INIT_AMB overlapping store
2. Bug #2: Fix AMB lane mask semantics
3. Bug #4: Implement anchor gating
4. Bug #3: Centralize classification

**Phase 2 - Refactoring for Maintainability** (enables review/testing):
5. Refactor #1-5: Extract all 9 `SS_*` functions
6. Refactor #6: Add dispatcher logic
7. Refactor #7: Unify naming conventions
8. Refactor #8: Implement lane mask semantics properly (completes bug #2)

**Phase 3 - Performance & Polish**:
9. Bug #6: Add ss_cond_codes caching
10. Refactor #9: Add caching helper
11. Bug #5: Verify HOM fallback (should be automatic from bug #3 fix)
12. Bug #7-8: Fix typo, verify rare_allele

**Phase 4 - Validation**:
13. Refactor #12-13: Verify pattern alignment and bookkeeping parity
14. Refactor #14: Add comprehensive unit tests
15. Re-enable backward pass assertions in `test_supersite_hmm_states`
16. Smoke test with real multiallelic VCF fixtures

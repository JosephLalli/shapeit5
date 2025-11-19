# Multi-Alphabet PBWT Implementation Plan

**Branch**: `multi-alphabet-pbwt`
**Goal**: Implement true multi-alphabet PBWT (Durbin-style) for supersite neighbor selection
**Status**: Planning phase
**Estimated effort**: 30-40 hours (~1 week)

---

## Executive Summary

Replace the current anchor-gated binary PBWT approach with a true multi-alphabet PBWT that uses C-way stable partitioning on the full multi-allelic state (4-bit class codes). This addresses the fundamental flaw where supersite neighbor selection currently only considers the anchor variant's REF/ALT status, ignoring sibling variants and losing allelic diversity information.

**Key insight**: SHAPEIT5 already has 4-bit class codes per haplotype per supersite (`packed_allele_codes`). The missing piece is generalizing the PBWT column update from 2-way partition to C-way partition where C ≤ 16.

---

## Background & Problem Statement

### Current Implementation Flaw

From `conditioning_set_selection.cpp:115-131`, the PBWT update partitions haplotypes using:

```cpp
if (!H_opt_var.get(l, alookup)) {  // Binary: 0 or 1
    A[u] = alookup; C[u] = p; u++; // REF partition
} else {
    B[v] = alookup; D[v] = q; v++; // ALT partition
}
```

For supersites:
- Only the **anchor** variant participates in PBWT (siblings are masked via `sites_pbwt_evaluation`)
- Anchor selection is based on MAF, not representativeness
- Haplotypes that match on sibling alleles but differ on anchor are excluded from neighbor set
- This violates the principle of selecting neighbors based on overall sequence similarity

### Proposed Solution

Use the full multi-allelic state as a single symbol:
- At supersite anchors: partition by 4-bit class code (REF=0, ALT1=1, ALT2=2, ...)
- At biallelic sites: preserve existing binary partition (zero regression)
- Complexity: O(M + K) per column where K ≤ 16, essentially same as O(M)

---

## Phase 1: Branch Setup & Prerequisites

### 1.1 Create Feature Branch

```bash
cd /mnt/d/shapeit5
git checkout where-are-the-errors
git checkout -b multi-alphabet-pbwt
git status  # Verify clean or commit untracked files
```

### 1.2 Baseline Testing

**Before any changes**, establish baseline:

1. **Compile current code**:
   ```bash
   make clean
   make -j8
   ```

2. **Run existing tests** (if available):
   ```bash
   # Record baseline test results
   ./run_tests.sh > baseline_tests.log 2>&1
   ```

3. **Create test dataset** with high supersite density:
   - Use 1000 Genomes Project data or similar
   - Filter for chromosomal region with >10% multi-allelic sites
   - Small enough for fast iteration (~1000 samples, 10k variants)

4. **Measure baseline PBWT neighbor quality**:
   - At supersite anchors, compute: % neighbors sharing full allelic state vs. anchor-only
   - Expected baseline: low concordance (design flaw)

---

## Phase 2: Core Implementation

### 2.1 Symbol Extraction Infrastructure

#### File: `conditioning_set_header.h`

**Location**: `/mnt/d/shapeit5/phase_common/src/containers/conditioning_set/conditioning_set_header.h`

**Add around line 80** (in `private:` section):

```cpp
// Multi-alphabet PBWT support
inline uint8_t getHaplotypeSymbol(int hap_idx, int locus_j) const;
```

**Add around line 50** (member variables, if not present):

```cpp
// Supersite context for MA-PBWT
const std::vector<SuperSite>* super_sites;
const std::vector<int>* locus_to_super_idx;
const std::vector<uint8_t>* packed_allele_codes;
```

#### File: `conditioning_set_selection.cpp`

**Location**: `/mnt/d/shapeit5/phase_common/src/containers/conditioning_set/conditioning_set_selection.cpp`

**Add before `select(int chunk)` method** (around line 98):

```cpp
inline uint8_t conditioning_set::getHaplotypeSymbol(int hap_idx, int locus_j) const {
    // Check if this locus is a supersite anchor
    if (supersite_anchor_redirect_enabled &&
        locus_j < static_cast<int>(locus_to_super_idx->size())) {
        int ss_idx = (*locus_to_super_idx)[locus_j];

        if (ss_idx >= 0) {
            // This is a supersite anchor - extract 4-bit class code
            const SuperSite& ss = (*super_sites)[ss_idx];
            return unpackSuperSiteCode(
                packed_allele_codes->data(),
                ss.panel_offset,
                hap_idx
            );
        }
    }

    // Standard biallelic site - extract from bitmatrix
    return H_opt_var.get(locus_j, hap_idx);  // Returns 0 or 1
}
```

**Critical details**:
- `unpackSuperSiteCode()` is defined in `super_site_accessor.h:150-178`
- Extracts 4-bit code from packed byte array (2 codes per byte)
- Returns 0 (REF), 1 (ALT1), 2 (ALT2), ..., or 255 (MISSING)

---

### 2.2 Multi-Alphabet PBWT Column Update

#### File: `conditioning_set_selection.cpp`

**Add new method before `select(int chunk)`** (around line 98):

```cpp
/**
 * Perform C-way stable partition for multi-allelic site.
 *
 * This generalizes the binary PBWT partition to handle C symbols (C ≤ 16).
 * Uses counting sort for O(M + C) complexity.
 *
 * @param locus_j Current locus index
 * @param K Number of distinct allele classes at this site (including REF)
 * @param A IN/OUT: Permutation array (haplotype indices in sorted order)
 * @param A_new Scratch buffer for new permutation
 * @param C IN/OUT: Divergence array (last position each hap differed from predecessor)
 * @param C_new Scratch buffer for new divergence
 * @param prev_rank Auxiliary array mapping hap_idx -> position in old A
 */
void conditioning_set::updateColumnMultiAlphabet(
    int locus_j,
    int K,
    std::vector<int>& A,
    std::vector<int>& A_new,
    std::vector<int>& C,
    std::vector<int>& C_new,
    std::vector<int>& prev_rank
) {
    int M = A.size();

    // Step 1: Count occurrences of each symbol (counting sort)
    int cnt[16] = {0};
    std::vector<uint8_t> symbols(M);  // Cache symbols to avoid repeated extraction

    for (int i = 0; i < M; i++) {
        int h = A[i];
        uint8_t symbol = getHaplotypeSymbol(h, locus_j);

        // Handle missing data: remap code 255 -> K (place at end)
        if (symbol == 255) symbol = K;

        symbols[i] = symbol;
        if (symbol <= K) cnt[symbol]++;
    }

    // Step 2: Compute prefix sums to get starting offset for each symbol block
    int offset[17];  // K+1 for missing data
    int running = 0;
    for (int s = 0; s <= K; s++) {
        offset[s] = running;
        running += cnt[s];
    }

    // Step 3: Stable partition - place haplotypes into A_new maintaining order within each symbol
    for (int i = 0; i < M; i++) {
        int h = A[i];
        uint8_t symbol = symbols[i];

        if (symbol <= K) {
            int pos = offset[symbol]++;
            A_new[pos] = h;
            prev_rank[h] = i;  // Remember this haplotype's position in old A
        }
    }

    // Step 4: Update divergence array
    // CRITICAL: This is the subtle part that maintains PBWT correctness
    C_new[0] = locus_j;  // First element always diverges at current locus (no predecessor)

    for (int i = 1; i < M; i++) {
        int h_curr = A_new[i];
        int h_prev = A_new[i-1];

        uint8_t s_curr = symbols_at_new_position[i];  // Need to cache or recompute
        uint8_t s_prev = symbols_at_new_position[i-1];

        if (s_curr == s_prev) {
            // Symbols match at this locus
            // Divergence should be inherited from wherever these two haplotypes
            // last differed in the previous PBWT order

            // Find where h_curr was positioned in the old A array
            int old_rank_curr = prev_rank[h_curr];

            // Carry forward the divergence from that old position
            // This represents "these haps have been identical since position C[old_rank_curr]"
            C_new[i] = C[old_rank_curr];

        } else {
            // Symbols differ at this locus - divergence resets here
            C_new[i] = locus_j;
        }
    }

    // Swap buffers (avoids copy)
    A.swap(A_new);
    C.swap(C_new);
}
```

**CRITICAL BUG FIX in above code**:

The code references `symbols_at_new_position[i]` which doesn't exist. We need to either:

**Option A**: Cache symbols indexed by haplotype ID:
```cpp
// After Step 3, create reverse lookup
std::vector<uint8_t> symbol_by_hap(M * 2);  // Max hap index
for (int i = 0; i < M; i++) {
    symbol_by_hap[A[i]] = symbols[i];
}

// In Step 4:
uint8_t s_curr = symbol_by_hap[h_curr];
uint8_t s_prev = symbol_by_hap[h_prev];
```

**Option B** (cleaner): Build symbol array for new positions:
```cpp
// After Step 3:
std::vector<uint8_t> symbols_new(M);
for (int i = 0; i < M; i++) {
    symbols_new[i] = getHaplotypeSymbol(A_new[i], locus_j);
    if (symbols_new[i] == 255) symbols_new[i] = K;
}

// In Step 4:
uint8_t s_curr = symbols_new[i];
uint8_t s_prev = symbols_new[i-1];
```

**Recommended: Option B** (more cache-friendly, simpler logic).

#### Corrected Full Implementation:

```cpp
void conditioning_set::updateColumnMultiAlphabet(
    int locus_j,
    int K,
    std::vector<int>& A,
    std::vector<int>& A_new,
    std::vector<int>& C,
    std::vector<int>& C_new,
    std::vector<int>& prev_rank
) {
    int M = A.size();

    // Step 1: Extract and cache symbols for all haplotypes
    std::vector<uint8_t> symbols(M);
    int cnt[17] = {0};  // K+1 for missing

    for (int i = 0; i < M; i++) {
        int h = A[i];
        uint8_t symbol = getHaplotypeSymbol(h, locus_j);
        if (symbol == 255) symbol = K;  // Remap missing to end
        symbols[i] = symbol;
        cnt[symbol]++;
    }

    // Step 2: Prefix sums
    int offset[17];
    int running = 0;
    for (int s = 0; s <= K; s++) {
        offset[s] = running;
        running += cnt[s];
    }

    // Step 3: Stable partition
    for (int i = 0; i < M; i++) {
        int h = A[i];
        uint8_t symbol = symbols[i];
        int pos = offset[symbol]++;
        A_new[pos] = h;
        prev_rank[h] = i;  // Track old position
    }

    // Step 4: Extract symbols at new positions
    std::vector<uint8_t> symbols_new(M);
    for (int i = 0; i < M; i++) {
        uint8_t symbol = getHaplotypeSymbol(A_new[i], locus_j);
        if (symbol == 255) symbol = K;
        symbols_new[i] = symbol;
    }

    // Step 5: Update divergence
    C_new[0] = locus_j;
    for (int i = 1; i < M; i++) {
        if (symbols_new[i] == symbols_new[i-1]) {
            // Match: inherit divergence from old position of current hap
            int h_curr = A_new[i];
            int old_rank = prev_rank[h_curr];
            C_new[i] = C[old_rank];
        } else {
            // Mismatch: diverge here
            C_new[i] = locus_j;
        }
    }

    // Swap
    A.swap(A_new);
    C.swap(C_new);
}
```

**Why divergence tracking is subtle**:

The divergence array `C[i]` means: "the position where `A[i]` last differed from `A[i-1]` in the co-lexicographic order". When symbols match at the current locus:

1. The two adjacent haplotypes in the new order (`A_new[i]` and `A_new[i-1]`) did NOT diverge at this locus
2. But they might have diverged at some earlier locus
3. We need to inherit the divergence from the old PBWT order
4. `prev_rank[A_new[i]]` tells us where `A_new[i]` was in the old `A` array
5. `C[prev_rank[...]]` gives us the old divergence value to carry forward

---

### 2.3 Modify Main PBWT Loop

#### File: `conditioning_set_selection.cpp`

**In `select(int chunk)` method, replace lines 107-136**:

**Current code**:
```cpp
for (int l = W.start_locus ; l <= W.stop_locus ; l ++) {
    int u = 0, v = 0, p = l + 1, q = l + 1;
    for (int i = 0 ; i < n_haps ; i ++) {
        int alookup = A[i];
        if (C[i] > p) p = C[i];
        if (D[i] > q) q = D[i];
        if (!H_opt_var.get(l, alookup)) {
            A[u] = alookup; C[u] = p; p = 0; u++;
        } else {
            B[v] = alookup; D[v] = q; q = 0; v++;
        }
    }
    for (int i = 0 ; i < v ; i ++) {
        A[u+i] = B[i]; C[u+i] = D[i];
    }
    if (sites_pbwt_selection[l]) store(l, A, C);
}
```

**Replace with**:

```cpp
// Allocate auxiliary buffers ONCE at start of method (BEFORE loop)
// Add after line 106 (before the for loop):
std::vector<int> A_new(n_haps);
std::vector<int> C_new(n_haps);
std::vector<int> prev_rank(n_haps * 2);  // Max haplotype index

for (int l = W.start_locus ; l <= W.stop_locus ; l ++) {
    // Determine if this is a multi-allelic supersite
    int ss_idx = -1;
    int K = 2;  // Default: biallelic (REF + ALT)

    if (supersite_anchor_redirect_enabled &&
        l < static_cast<int>(locus_to_super_idx->size())) {
        ss_idx = (*locus_to_super_idx)[l];
        if (ss_idx >= 0) {
            const SuperSite& ss = (*super_sites)[ss_idx];
            K = ss.n_classes;  // 1 + n_alts (e.g., 3 for REF + ALT1 + ALT2)
        }
    }

    if (K > 2) {
        // Multi-alphabet PBWT update
        updateColumnMultiAlphabet(l, K, A, A_new, C, C_new, prev_rank);
    } else {
        // Standard binary partition (PRESERVE EXACT CURRENT BEHAVIOR)
        int u = 0, v = 0, p = l + 1, q = l + 1;
        for (int i = 0 ; i < n_haps ; i ++) {
            int alookup = A[i];
            if (C[i] > p) p = C[i];
            if (D[i] > q) q = D[i];
            if (!H_opt_var.get(l, alookup)) {
                A[u] = alookup; C[u] = p; p = 0; u++;
            } else {
                B[v] = alookup; D[v] = q; q = 0; v++;
            }
        }
        for (int i = 0 ; i < v ; i ++) {
            A[u+i] = B[i]; C[u+i] = D[i];
        }
    }

    // Neighbor selection (unchanged)
    if (sites_pbwt_selection[l]) store(l, A, C);
}
```

**Critical details**:
- Buffer allocation OUTSIDE loop (avoid repeated malloc)
- Exact preservation of binary path when `K <= 2` (zero regression)
- `store()` function unchanged (agnostic to how `A` and `C` were built)

---

### 2.4 Connect Supersite Data to conditioning_set

#### File: `conditioning_set_header.h`

**Add public method** (around line 70):

```cpp
void setSupersiteContext(
    const std::vector<SuperSite>* ss,
    const std::vector<int>* locus_map,
    const std::vector<uint8_t>* codes
);
```

#### File: `conditioning_set_managment.cpp`

**Add implementation** (end of file, around line 155):

```cpp
void conditioning_set::setSupersiteContext(
    const std::vector<SuperSite>* ss,
    const std::vector<int>* locus_map,
    const std::vector<uint8_t>* codes
) {
    super_sites = ss;
    locus_to_super_idx = locus_map;
    packed_allele_codes = codes;
}
```

#### File: `phaser/phaser_algorithm.cpp` or wherever PBWT is initialized

**Find where `H.initialize()` is called** and add:

```cpp
// After H.initialize(...)
H.setSupersiteContext(&super_sites, &locus_to_super_idx, &packed_allele_codes);
```

**Verify** these vectors exist in `phaser` class (from exploration, they should).

---

## Phase 3: Edge Case Handling

### 3.1 Missing Data (Code 255)

**Already handled** in `updateColumnMultiAlphabet()`:

```cpp
if (symbol == 255) symbol = K;  // Remap missing to class K (end of partition)
```

**Rationale**:
- Missing data doesn't match any valid class
- Placing at end ensures missing haplotypes cluster together
- Divergence tracking still works (missing vs missing = match, missing vs valid = mismatch)

### 3.2 Degenerate Cases

**K = 1** (all haplotypes REF at supersite):
- Should never occur in practice (why would it be multi-allelic?)
- If it happens: `updateColumnMultiAlphabet()` degenerates to no-op (all go in class 0 bucket)

**K = 2 at a supersite**:
- Biallelic path handles it correctly
- Example: REF + ALT1 (ALT2 has zero frequency)

**K > 16**:
- Cannot happen: `SUPERSITE_MAX_ALTS = 15` (line 70 of `super_site_accessor.h`)
- Plus REF → max K = 16

### 3.3 IBD2 Exclusion

**In `store()` function** (lines 138-178 of `conditioning_set_selection.cpp`):

```cpp
if (!Kbanned.noIBD2(tar_hap, A[i], l)) continue;  // Skip IBD2 neighbors
```

**Analysis**:
- `Kbanned.noIBD2()` checks if two haplotypes share IBD2 at locus `l`
- Uses haplotype indices and locus index only
- **No dependence on how A was sorted**
- Works identically for binary and multi-alphabet PBWT

**Conclusion**: No changes needed.

### 3.4 Thread Safety

**PBWT state per thread**:
- `A`, `B`, `C`, `D` are local to `select(int chunk)`
- Each thread processes a different chunk (non-overlapping loci)
- Auxiliary buffers (`A_new`, `C_new`, `prev_rank`) are also local

**Global data (read-only during PBWT)**:
- `H_opt_var`, `super_sites`, `packed_allele_codes`
- All const during neighbor selection phase

**Conclusion**: Thread-safe by design.

---

## Phase 4: Testing & Validation

### 4.1 Unit Tests

**Create**: `/mnt/d/shapeit5/test/test_ma_pbwt.cpp`

#### Test 1: Biallelic Unchanged

```cpp
TEST(MA_PBWT, BiallelicUnchanged) {
    // Create small biallelic dataset
    // Run binary PBWT
    // Run MA-PBWT with K=2
    // Assert: A and C arrays identical
}
```

#### Test 2: Triallelic Sorting

```cpp
TEST(MA_PBWT, TriallelicPartition) {
    // Manual example: 10 haplotypes, symbols [0,0,1,1,1,2,2,0,1,2]
    // Expected partition: [0,0,0 | 1,1,1,1 | 2,2,2] (stable sort)
    // Verify A is correctly partitioned
}
```

#### Test 3: Divergence Propagation

```cpp
TEST(MA_PBWT, DivergenceTracking) {
    // Two columns:
    //   Col 0: [0,0,1,1] → A = [h0,h1,h2,h3], C = [0,0,0,0]
    //   Col 1: [0,1,0,1] (triallelic) → A should update, C should track correctly
    // Verify C values match expected divergence positions
}
```

#### Test 4: Missing Data

```cpp
TEST(MA_PBWT, MissingHandling) {
    // Symbols: [0, 1, 255, 2, 255]
    // Expected partition: [0 | 1 | 2 | 255, 255]
    // Verify missing placed at end
}
```

#### Test 5: Edge Case K=16

```cpp
TEST(MA_PBWT, MaxAlphabet) {
    // 16 classes (REF + 15 ALTs)
    // Verify no overflow, correct partition
}
```

### 4.2 Integration Tests

#### Dataset 1: High Supersite Density

**Source**: 1000 Genomes Project, chr6 MHC region (highly polymorphic)

**Pipeline**:
```bash
# Extract test region
bcftools view -r 6:28000000-33000000 1000G.bcf -Ob -o test_mhc.bcf

# Run SHAPEIT5 with baseline (current code)
./phase_common --input test_mhc.bcf --output baseline_phased.bcf

# Run SHAPEIT5 with MA-PBWT
git checkout multi-alphabet-pbwt
make clean && make -j8
./phase_common --input test_mhc.bcf --output ma_pbwt_phased.bcf

# Compare switch error rates
./switch_error.py baseline_phased.bcf ma_pbwt_phased.bcf truth.bcf
```

**Expected**: MA-PBWT ≥ baseline accuracy (ideally better at supersites)

#### Dataset 2: Typical Dataset

**Source**: Standard WGS, ~95% biallelic

**Metrics**:
1. **Runtime**: MA-PBWT should be within 5% of baseline
2. **Memory**: No significant increase
3. **Accuracy**: No regression

#### Neighbor Quality Metric

**At each supersite anchor**:

```python
for each sample S:
    full_state = get_multiallelic_genotype(S, supersite)

    # Get K neighbors selected by PBWT
    neighbors = get_pbwt_neighbors(S, supersite)

    # Count how many neighbors match full state vs anchor-only
    full_match = sum(1 for N in neighbors if get_multiallelic_genotype(N, supersite) == full_state)
    anchor_match = sum(1 for N in neighbors if get_anchor_allele(N, supersite) == get_anchor_allele(S, supersite))

    concordance = full_match / len(neighbors)
```

**Hypothesis**:
- Baseline (binary PBWT): ~50% full concordance (random among anchor-matching neighbors)
- MA-PBWT: >70% full concordance (actively selects full-state matches)

### 4.3 Regression Testing

**Test suite**:
1. Run all existing SHAPEIT5 tests
2. Compare output to baseline on:
   - Standard biallelic VCF
   - Mixed biallelic + multi-allelic VCF
   - Edge cases: singletons, all-missing sites, etc.

**Success criteria**:
- Biallelic-only datasets: bit-identical output
- Mixed datasets: no accuracy regression

---

## Phase 5: Performance Optimization

### 5.1 Profiling

**Tools**:
```bash
# Compile with profiling
make clean
CXXFLAGS="-O3 -g -pg" make -j8

# Run on test dataset
./phase_common --input test.bcf --output out.bcf

# Analyze
gperf gmon.out
```

**Likely hotspots**:
1. `getHaplotypeSymbol()` - called M times per locus
2. `unpackSuperSiteCode()` - bit unpacking overhead
3. Symbol array allocation in `updateColumnMultiAlphabet()`

### 5.2 Optimization: Batch Symbol Extraction

**Current** (in `updateColumnMultiAlphabet()`):
```cpp
for (int i = 0; i < M; i++) {
    symbols[i] = getHaplotypeSymbol(A[i], locus_j);
}
```

**Issue**: Each call does supersite lookup + byte unpacking

**Optimized**: Precompute symbol lookup table per locus:

```cpp
// In select(int chunk), before locus loop:
std::vector<uint8_t> symbol_lut(n_haps * 2);  // Lookup table for current locus

// In locus loop, before partition:
if (K > 2) {
    // Precompute symbols for all haplotypes at this supersite
    const SuperSite& ss = (*super_sites)[ss_idx];
    for (int h = 0; h < n_haps * 2; h++) {
        symbol_lut[h] = unpackSuperSiteCode(packed_allele_codes->data(), ss.panel_offset, h);
    }

    // Pass symbol_lut to updateColumnMultiAlphabet, use direct lookup instead of getHaplotypeSymbol()
}
```

**Benefit**: Single pass over `packed_allele_codes` per locus instead of M random accesses.

### 5.3 Optimization: Cache-Friendly Layout

**Current**: `packed_allele_codes` is locus-major (all haplotypes for locus L, then locus L+1, ...)

**Potential**: Transpose to haplotype-major for sequential access during PBWT sweep

**Analysis**:
- PBWT sweeps locus-by-locus → current layout is optimal
- Transposing would help HMM (which sweeps haplotype-by-haplotype) but hurt PBWT
- **Don't change** unless profiling shows major cache misses

### 5.4 SIMD Opportunities

**Counting sort** could use SIMD for K ≤ 8:
- Vectorize symbol extraction
- Parallel histogram computation

**Effort vs Reward**:
- High implementation complexity
- Marginal benefit (K-way partition is already O(M+K) with small K)
- **Not recommended** unless profiling shows >10% time in counting sort

---

## Phase 6: Documentation & Code Quality

### 6.1 Inline Comments

**Add detailed comments** in:

1. **`updateColumnMultiAlphabet()`**:
   - Explain counting sort algorithm
   - Clarify divergence inheritance logic (subtle!)
   - Document missing data handling

2. **`getHaplotypeSymbol()`**:
   - Why we check supersite context
   - Bit unpacking format reference

3. **Main PBWT loop**:
   - Rationale for K > 2 vs K <= 2 split
   - Performance considerations

### 6.2 High-Level Documentation

**Update** (or create): `/mnt/d/shapeit5/docs/PBWT_ALGORITHM.md`

Contents:
- Overview of PBWT in SHAPEIT5
- Binary PBWT algorithm
- Multi-alphabet extension for supersites
- Complexity analysis
- References to Durbin 2014 and Li & Stephens 2003

### 6.3 Commit Strategy

**Atomic commits** (easy to review and revert):

```bash
# Commit 1
git add conditioning_set_header.h conditioning_set_selection.cpp
git commit -m "Add symbol extraction infrastructure for multi-alphabet PBWT

- Add getHaplotypeSymbol() method to extract 4-bit supersite codes
- Add supersite context pointers to conditioning_set
- No functional changes yet"

# Commit 2
git add conditioning_set_selection.cpp
git commit -m "Implement C-way PBWT partition with divergence tracking

- Add updateColumnMultiAlphabet() for K-way stable partition
- Use counting sort for O(M+K) complexity
- Correctly propagate divergence array when symbols match"

# Commit 3
git add conditioning_set_selection.cpp conditioning_set_managment.cpp
git commit -m "Integrate multi-alphabet PBWT into neighbor selection

- Modify select(int chunk) to dispatch binary vs MA-PBWT
- Preserve exact binary behavior for K<=2 (zero regression)
- Connect supersite data via setSupersiteContext()"

# Commit 4
git add test/test_ma_pbwt.cpp
git commit -m "Add unit tests for multi-alphabet PBWT

- Test biallelic path unchanged
- Test triallelic partition correctness
- Test divergence propagation
- Test missing data handling"

# Commit 5
git add test/integration/
git commit -m "Add integration tests and benchmarks

- High supersite density dataset
- Neighbor quality metrics
- Performance regression tests"
```

### 6.4 Code Review Checklist

Before requesting review:

- [ ] All tests pass
- [ ] No compiler warnings
- [ ] Valgrind clean (no memory leaks)
- [ ] Code formatted consistently
- [ ] Comments explain "why" not just "what"
- [ ] Performance benchmarks documented
- [ ] Git history clean (no "fix typo" commits)

---

## Risk Assessment & Mitigation

### Risk 1: Divergence Tracking Bug 🔴 HIGH RISK

**Symptom**: Incorrect `C` values → wrong neighbors selected → poor phasing accuracy

**Root cause**: Divergence inheritance logic is subtle (see Phase 2.2)

**Mitigation**:
1. **Unit test** with manual calculation of expected C values
2. **Visualization**: Print A and C arrays for small test case, verify by hand
3. **Assertion**: Add debug assert that C[i] <= current_locus (divergence can't be in future)

**Test case**:
```cpp
// Column 0: symbols [0, 0, 1, 1]
//   A = [h0, h1, h2, h3]
//   C = [0, 0, 0, 0]
// Column 1: symbols [0, 1, 0, 1] (h0=0, h1=1, h2=0, h3=1)
//   A_new should be: [h0, h2, h1, h3]  (0s first, 1s second, stable)
//   C_new should be: [1, C[old_rank[h2]], 1, C[old_rank[h3]]]
//           = [1, C[2]=0, 1, C[3]=0]
//           = [1, 0, 1, 0]
```

### Risk 2: Performance Regression ⚠️ MEDIUM RISK

**Symptom**: PBWT neighbor selection takes >5% longer

**Root cause**: Symbol extraction overhead, especially `unpackSuperSiteCode()`

**Mitigation**:
1. **Profile early**: Baseline vs MA-PBWT runtime on typical dataset
2. **Optimize hot path**: See Phase 5.2 (batch symbol extraction)
3. **Binary fast path**: K <= 2 uses original code (zero overhead)

**Acceptable overhead**: <5% on typical datasets (95% biallelic)

### Risk 3: Neighbor Quality Doesn't Improve ⚠️ MEDIUM RISK

**Symptom**: MA-PBWT selects neighbors with same full-state concordance as baseline

**Possible causes**:
1. Implementation bug (e.g., not actually using multi-allelic codes)
2. Supersite structure doesn't provide useful signal
3. K too small (most supersites are biallelic)

**Mitigation**:
1. **Sanity check**: Print symbols extracted at supersites, verify they're not all 0/1
2. **Direct test**: Manually construct high-K supersite, verify partition
3. **Biological validation**: Use MHC region (known high diversity)

**Fallback**: If no improvement, flag this as "infrastructure for future work" (still valuable to have MA-PBWT capability)

### Risk 4: Memory Bloat ⚠️ LOW RISK

**Symptom**: Memory usage increases significantly

**Analysis**:
- Persistent: Still just `A` and `C` (size M)
- Temporary per chunk: `A_new`, `C_new`, `prev_rank`, `symbols` (total ~4M ints + M bytes)
- For M=10k haplotypes: ~160 KB per thread

**Mitigation**: Negligible vs. main memory consumers (haplotype storage, HMM state)

### Risk 5: Breaking Biallelic Behavior ✅ ZERO RISK

**Why zero**: Binary path is **unchanged** when `K <= 2`

**Verification**:
```cpp
if (K > 2) {
    // New code
} else {
    // Exact copy of original code
}
```

As long as conditional is correct, impossible to break biallelic sites.

---

## Success Criteria

### Must Have ✅

1. **Zero regression**: All existing tests pass
2. **Correctness**: MA-PBWT unit tests pass
3. **No performance degradation**: <5% overhead on biallelic-only datasets
4. **Memory bounded**: <2% increase in peak memory usage

### Should Have 🎯

1. **Neighbor quality improvement**: ≥20% increase in full-state concordance at supersites
2. **Phasing accuracy**: No decrease in switch error rate on mixed datasets
3. **Code quality**: Clean git history, comprehensive comments
4. **Documentation**: Updated algorithm description

### Nice to Have 🌟

1. **Phasing accuracy improvement**: Measurable decrease in switch error rate at supersites
2. **Performance optimization**: <1% overhead on typical datasets
3. **Scalability**: Tested on large dataset (100k samples, 1M variants)

---

## Timeline Estimate

| Phase | Task | Estimated Time |
|-------|------|---------------|
| 1 | Branch setup & baseline testing | 1 hour |
| 2.1 | Symbol extraction infrastructure | 2 hours |
| 2.2 | Multi-alphabet partition implementation | 4 hours |
| 2.3 | Integrate into main PBWT loop | 2 hours |
| 2.4 | Connect supersite data | 1 hour |
| 3 | Edge case handling | 3 hours |
| 4.1 | Unit tests | 4 hours |
| 4.2 | Integration tests | 6 hours |
| 4.3 | Regression testing | 2 hours |
| 5 | Profiling & optimization | 4-8 hours |
| 6 | Documentation & cleanup | 2 hours |

**Total**: 31-35 hours (~1 week of focused work)

**Breakdown by risk**:
- Core implementation (low risk, well-defined): 10 hours
- Testing & validation (medium risk, may reveal bugs): 12 hours
- Optimization (variable, depends on profiling): 4-8 hours
- Documentation: 2 hours

---

## References

1. **Durbin, R. (2014)**: "Efficient haplotype matching and storage using the positional Burrows-Wheeler transform (PBWT)". *Bioinformatics* 30(9):1266-1272.
   - Original PBWT algorithm for binary sequences
   - Multi-alphabet extension mentioned but not detailed

2. **Li, N. & Stephens, M. (2003)**: "Modeling linkage disequilibrium and identifying recombination hotspots using single-nucleotide polymorphism data". *Genetics* 165(4):2213-2233.
   - Li & Stephens HMM used in SHAPEIT5

3. **SHAPEIT5 codebase**: `/mnt/d/shapeit5/SHAPEIT5_algorithms_description.md`
   - Section 2.1.2: "Anchor-gated PBWT" (current flaw)
   - Section 2.2: Supersite infrastructure

4. **Counting sort**: CLRS "Introduction to Algorithms" Chapter 8.2
   - O(n + k) stable sort for small integer keys
   - Used for C-way partition

---

## Appendix A: Key File Locations

| File | Lines | Purpose |
|------|-------|---------|
| `conditioning_set_header.h` | 88 | Class definition, add symbol extraction method |
| `conditioning_set_selection.cpp` | 230 | **Main PBWT loop** - modify lines 107-136 |
| `conditioning_set_managment.cpp` | 155 | Initialization, add setSupersiteContext() |
| `super_site_accessor.h` | 268 | SuperSite struct, unpackSuperSiteCode() |
| `super_site_builder.cpp` | 193 | Builds packed_allele_codes (no changes) |
| `bitmatrix.h` | 76 | Binary storage (keep unchanged) |
| `haplotype_segment_single.h` | 1695 | HMM (uses PBWT neighbors, no changes) |

---

## Appendix B: Debugging Tips

### Print PBWT State

Add to `select(int chunk)`:

```cpp
if (chunk == 0 && l == 100) {  // Debug first chunk, locus 100
    std::cerr << "Locus " << l << " (K=" << K << ")\n";
    std::cerr << "A: ";
    for (int i = 0; i < 10; i++) std::cerr << A[i] << " ";
    std::cerr << "\nC: ";
    for (int i = 0; i < 10; i++) std::cerr << C[i] << " ";
    std::cerr << "\n";
}
```

### Verify Symbol Extraction

```cpp
if (ss_idx >= 0) {
    const SuperSite& ss = (*super_sites)[ss_idx];
    std::cerr << "Supersite " << ss_idx << ": " << (int)ss.n_classes << " classes\n";
    for (int h = 0; h < 10; h++) {
        uint8_t symbol = getHaplotypeSymbol(h, l);
        std::cerr << "  Hap " << h << ": symbol " << (int)symbol << "\n";
    }
}
```

### Check Divergence Invariant

```cpp
// In updateColumnMultiAlphabet, after C_new computation:
for (int i = 0; i < M; i++) {
    assert(C_new[i] <= locus_j);  // Divergence can't be in future
}
```

---

## Appendix C: Testing Datasets

### Dataset 1: Synthetic High-K Supersite

Create VCF with 1000 samples, 100 variants, including:
- Variant 50: 10-allelic (REF + 9 ALTs), frequencies [50%, 10%, 5%, 5%, 5%, 5%, 5%, 5%, 5%, 5%]
- Check PBWT neighbors for samples with rare alleles (ALT9)
- Expect: MA-PBWT selects other ALT9 carriers as neighbors

### Dataset 2: 1000 Genomes MHC

```bash
# Download
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz

# Extract MHC region
bcftools view -r 6:28000000-33000000 -Ob -o mhc.bcf ALL.chr6...vcf.gz

# Count supersites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' mhc.bcf | awk -F'\t' '{n=split($4,a,","); if(n>1) count++} END {print count}'
```

Expected: >1000 multi-allelic sites in MHC region

### Dataset 3: Biallelic-Only Control

```bash
# Create biallelic-only VCF
bcftools view -m2 -M2 original.vcf.gz -Ob -o biallelic_only.bcf
```

Expected: Bit-identical output from baseline vs MA-PBWT

---

## Status Tracking

**Current status**: Planning complete, ready for implementation

**Next steps**:
1. Create branch `multi-alphabet-pbwt`
2. Implement Phase 2.1 (symbol extraction)
3. Write unit test for symbol extraction
4. Implement Phase 2.2 (partition function)
5. Write unit test for partition
6. ... (continue iteratively)

**Estimated completion**: TBD

---

*This plan document is version-controlled in the SHAPEIT5 repository for future reference.*

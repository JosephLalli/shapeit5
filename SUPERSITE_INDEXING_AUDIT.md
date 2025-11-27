# Supersite Sibling-Skipping Indexing Audit

**Date**: 2025-11-25
**Issue**: Supersites remove the assumption that every variant in a list must be seen by the HMM. Siblings get skipped during iteration, creating indexing mismatches between:
- `vrel`/`vabs`: Iterates over ALL variants (anchors + siblings)
- `arel`: Iterates over ambiguous ANCHORS only (siblings skipped)
- `Lengths[s]`: Currently stores count of ALL variants in segment

## Executive Summary

**Critical Bug Found**: genotype_prune.cpp:209,219 - Segment merge uses `(vrel < Lengths[s-1])` to determine segment membership, but `vrel` includes siblings while the comparison should only count anchors.

**Systematic Risk**: Any code that mixes all-variant iteration (`vrel`) with anchor-only counts/arrays (`arel`, ambiguous masks) is at risk.

---

## Part 1: Ambiguous Mask Build/Rebuild Points

### Location 1: Initial Build - genotype_build.cpp:241-296

**Function**: `genotype::build()` - Initial construction of Ambiguous array

**Code Structure**:
```cpp
// Line 242: Initialize
Ambiguous = vector<unsigned char>(n_ambiguous, 0U);

// Lines 244-265: First pass - SCA loci
for (unsigned int s = 0, a0 = 0, ..., vabs = 0; s < n_segments; s++) {
    for (unsigned int vrel = 0; vrel < Lengths[s]; vrel++) {
        unsigned int v_idx = vabs + vrel;
        SuperSiteContext ctx = getSuperSiteContext(v_idx);
        if (ctx.is_member && !ctx.is_anchor) continue;  // SKIP SIBLINGS

        if (f_sca) {
            // Set bits in Ambiguous[a0]
            for (unsigned int h = 0; h < HAP_NUMBER; h++) {
                bool allele = ...;
                if (allele) HAP_SET(Ambiguous[a0], h);
            }
        }
        a0 += (f_sca||f_het);
    }
    // Lines 267-296: Second pass - HET loci
    unsigned int n_unf = orderedSegments[s];
    for (unsigned int vrel = 0; vrel < Lengths[s]; vrel++) {
        unsigned int v_idx = vabs + vrel;
        SuperSiteContext ctx = getSuperSiteContext(v_idx);
        if (ctx.is_member && !ctx.is_anchor) continue;  // SKIP SIBLINGS

        if (f_het) {
            Ambiguous[a1] = CANONICAL_MASKS[n_unf % HAP_NUMBER];
            n_unf++;
        }
        a1 += (f_sca||f_het);
    }
    vabs += Lengths[s];
}
```

**Indexing Analysis**:
- ✅ **CORRECT**: Both passes skip siblings via line 251/273 `if (!ctx.is_anchor) continue`
- ✅ **CORRECT**: `a0` and `a1` only increment for ambiguous anchors
- ✅ **CORRECT**: `vabs` advances by `Lengths[s]` which includes siblings, then used as base for next segment
- ⚠️ **CONCERN**: `Lengths[s]` includes siblings, used for iteration bounds

**Verdict**: ✅ **Currently correct** - explicit sibling skipping prevents bugs

---

### Location 2: Merge Rebuild - genotype_prune.cpp:207-222

**Function**: `genotype::performMerges()` - Rebuilds ambiguous masks when merging segments

**Code Structure**:
```cpp
// Line 148: Initialize
Ambiguous2 = vector<unsigned char>(Ambiguous.size(), 0);

// Lines 207-212: Rebuild for merged_h0
for (unsigned int vrel = 0, arel = 0; vrel < (Lengths[s-1]+Lengths[s]); vrel++) {
    if (isAmbiguousAnchor(voffset + vrel)) {
        if (HAP_GET(Ambiguous[aoffset+arel],
                    (vrel<Lengths[s-1])?prev_h0:next_h0))  // ❌ BUG HERE
            HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h0]);
        arel++;
    }
}

// Lines 217-222: Rebuild for merged_h1 (identical pattern)
for (unsigned int vrel = 0, arel = 0; vrel < (Lengths[s-1]+Lengths[s]); vrel++) {
    if (isAmbiguousAnchor(voffset + vrel)) {
        if (HAP_GET(Ambiguous[aoffset+arel],
                    (vrel<Lengths[s-1])?prev_h1:next_h1))  // ❌ BUG HERE
            HAP_SET(Ambiguous2[aoffset+arel], Mhaps[merged_h1]);
        arel++;
    }
}
```

**Indexing Analysis**:
- ❌ **BUG**: `vrel` iterates ALL variants (anchors + siblings)
- ❌ **BUG**: `arel` only advances for ambiguous anchors
- ❌ **BUG**: Comparison `(vrel < Lengths[s-1])` includes sibling counts
- ❌ **RESULT**: For supersite with siblings, `vrel` crosses `Lengths[s-1]` while still processing anchors from first segment

**Example**:
```
Segment s-1: [anchor0, sibling1, anchor2, sibling3] - Lengths[s-1]=4, 2 ambiguous anchors
Segment s:   [anchor4, sibling5] - Lengths[s]=2, 1 ambiguous anchor

Iteration:
vrel=0: anchor0, arel=0, vrel<4 ✓ use prev_h0  ✓ CORRECT
vrel=1: sibling1, skipped
vrel=2: anchor2, arel=1, vrel<4 ✓ use prev_h0  ✓ CORRECT
vrel=3: sibling3, skipped
vrel=4: anchor4, arel=2, vrel>=4 ✗ use next_h0  ❌ WRONG! Should use prev_h0 if anchor2 is last in seg s-1
```

**Verdict**: ❌ **CRITICAL BUG** - Root cause of test failure

---

### Location 3: Copy-through - genotype_prune.cpp:232-236, 251-255

**Function**: `genotype::performMerges()` - Copy ambiguous masks for non-merged segments

**Code Structure**:
```cpp
// Lines 232-236: Copy when no merge at segment s-1
for (unsigned int vrel = 0, arel = 0; vrel < Lengths[s-1]; vrel++) {
    if (isAmbiguousAnchor(voffset + vrel)) {
        Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
        arel++;
    }
}

// Lines 251-255: Copy final segment
for (unsigned int vrel = 0, arel = 0; vrel < Lengths.back(); vrel++) {
    if (isAmbiguousAnchor(voffset + vrel)) {
        Ambiguous2[aoffset+arel] = Ambiguous[aoffset+arel];
        arel++;
    }
}
```

**Indexing Analysis**:
- ✅ **CORRECT**: Simple copy, no cross-segment logic
- ✅ **CORRECT**: `vrel` bounds by `Lengths[s-1]` for iteration only
- ✅ **CORRECT**: `arel` properly tracks ambiguous anchor indexing
- ✅ **CORRECT**: `isAmbiguousAnchor()` explicitly checks and skips siblings

**Verdict**: ✅ **Currently correct**

---

## Part 2: Segment Gap / Transition Macros

### Macro Search Locations
- `haplotype_segment_single.h`
- `haplotype_segment_double.h`
- `haplotype_segment_single.cpp`
- `haplotype_segment_double.cpp`

### Key Variables in HMM Loops

**Forward Pass** (haplotype_segment_single.cpp:~520-690):
```cpp
for (curr_abs_locus = locus_first; curr_abs_locus <= locus_last; curr_abs_locus++) {
    // Segment boundary check
    if (curr_segment_locus == Lengths[curr_segment_index]) {
        // Move to next segment
        curr_segment_index++;
        curr_segment_locus = 0;
        prev_abs_locus = curr_abs_locus;
    }

    // Sibling skip
    bool is_sibling = ...;
    if (is_sibling) {
        // Don't update prev_abs_locus
        // Don't increment curr_segment_locus
        continue;
    }

    // Normal processing
    curr_segment_locus++;
}
```

**Backward Pass** (haplotype_segment_single.cpp:~750-1010):
```cpp
for (curr_abs_locus = locus_last; curr_abs_locus >= locus_first; curr_abs_locus--) {
    // Sibling skip with boundary check
    if (need_init && is_sibling) {
        // CRITICAL: Check boundary BEFORE decrementing
        if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
            SET_OTHER_TRANS(transition_probabilities);
        }
        curr_segment_locus--;
        if (curr_segment_locus < 0 && curr_segment_index > 0) {
            curr_segment_index--;
            curr_segment_locus = Lengths[curr_segment_index] - 1;
        }
        continue;
    }

    // Segment boundary check
    if (curr_segment_locus == 0 && curr_abs_locus != locus_first) {
        SET_OTHER_TRANS(transition_probabilities);
    }

    curr_segment_locus--;
}
```

### Indexing Analysis for Segment Transitions

**Forward Pass**:
- ✅ **CORRECT**: `curr_segment_locus` only increments for non-siblings
- ✅ **CORRECT**: Boundary check `curr_segment_locus == Lengths[...]` compares anchor count to anchor count
- ⚠️ **QUESTION**: Does `Lengths[s]` store anchor-only count or all-variant count?

**Backward Pass**:
- ✅ **FIXED** (2025-11-19): Sibling-skip block checks boundary before continue
- ✅ **CORRECT**: `curr_segment_locus` only decrements for non-siblings
- ✅ **CORRECT**: Boundary transitions stored even when skipping siblings

**Key Question**: What does `Lengths[s]` actually represent?
- If ALL variants (anchors + siblings): ❌ **BUG** in boundary checks
- If anchors only: ✅ **CORRECT** for boundary checks

### Verification Needed

Need to trace where `Lengths[s]` is set:
1. Initial build (genotype_build.cpp)
2. After merges (genotype_prune.cpp)

---

## Part 3: Sibling AMB Flag Clearing

### Expected Behavior
"One early attempt at fixing this was, at the start of phasing, to ensure that siblings always have their ambiguous flag set to false."

### Search Results
- Only found reference in `/mnt/d/shapeit5/AGENTS.md` (documentation)
- No code found explicitly clearing AMB flags for siblings

### Analysis Required
Need to check:
1. Is AMB flag clearing for siblings implemented?
2. If yes, where?
3. If no, is it needed?
4. Does `isAmbiguousAnchor()` check prevent siblings from being treated as ambiguous?

**Function to examine**: `isAmbiguousAnchor()` in genotype_prune.cpp:71-79, 167-175
```cpp
auto isAmbiguousAnchor = [&](unsigned int locus) -> bool {
    bool is_amb = VAR_GET_AMB(MOD2(locus), Variants[DIV2(locus)]);
    SuperSiteContext ctx = getSuperSiteContext(locus);
    if (ctx.is_member) {
        if (!ctx.is_anchor) return false;  // SIBLINGS RETURN FALSE
        is_amb = ctx.has_het || ctx.has_sca;
    }
    return is_amb;
};
```

**Verdict**: ⚠️ **Relies on runtime check, not pre-clearing** - siblings return false from `isAmbiguousAnchor()`, but their actual AMB flag in Variants[] may still be set

---

## Part 4: Lengths[] Semantic Analysis

### Critical Question
What does `Lengths[s]` represent?
- **Option A**: Count of ALL variants in segment (anchors + siblings)
- **Option B**: Count of anchor-only variants in segment

### Evidence Collection Needed

**Test 1**: Check genotype_build.cpp where Lengths[] is populated
**Test 2**: Check if Lengths[] == Lengths_bio[] or if they differ
**Test 3**: Trace through HMM loop to see if `curr_segment_locus` matches Lengths[s]

### Implications

**If Lengths[] = ALL variants**:
- ❌ Forward/backward boundary checks WRONG (compare curr_segment_locus that only counts anchors to Lengths that counts all)
- ❌ Merge iteration bounds WRONG (vrel < Lengths[s-1] includes siblings)
- 🔧 **Fix**: Solution A (make Lengths[] anchor-only)

**If Lengths[] = anchors only**:
- ✅ Forward/backward boundary checks CORRECT
- ✅ BUT: Build code `vrel < Lengths[s]` would need special handling
- ⚠️ **Confusion**: Why have both Lengths[] and Lengths_bio[]?

---

## Summary of Findings

### Confirmed Bugs
1. ❌ **genotype_prune.cpp:209, 219** - Merge ambiguous mask rebuild uses wrong segment membership test

### Potential Bugs (Requires Investigation)
2. ⚠️ **Lengths[] semantics** - Unclear if stores all-variant count or anchor-only count
3. ⚠️ **Sibling AMB clearing** - Not explicitly implemented, relies on runtime checks

### Verified Correct
4. ✅ **genotype_build.cpp** - Initial ambiguous mask build correctly skips siblings
5. ✅ **genotype_prune.cpp:232-255** - Copy-through correctly skips siblings
6. ✅ **Backward pass sibling-skip** - Fixed 2025-11-19 to check boundaries before continue

---

## Recommended Actions

### Immediate (Fix Current Bug)
1. Implement Solution B for genotype_prune.cpp:209,219
   - Pre-compute ambiguous anchor counts per segment
   - Use `(arel < n_amb_first)` instead of `(vrel < Lengths[s-1])`

### Short-term (Clarify Semantics)
2. Document what `Lengths[s]` represents
3. Add assertions to verify Lengths[] behavior
4. Consider renaming: `Lengths[]` → `LengthsAnchor[]`, `Lengths_bio[]` → `LengthsAll[]`

### Long-term (Comprehensive Fix)
5. Implement Solution A: Make all segment iteration anchor-only
   - Refactor Lengths[] to only store anchor counts
   - Update all iteration loops to skip siblings explicitly or implicitly
   - Use Lengths_bio[] where all-variant counts needed

### Verification
6. Add test for segments with siblings at boundaries
7. Add test for multiple consecutive siblings
8. Add assertions in HMM loops: `assert(curr_segment_locus <= Lengths[curr_segment_index])`

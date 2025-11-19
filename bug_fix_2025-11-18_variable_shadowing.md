# Bug Fix: Variable Shadowing in Backward Pass

**Date:** 2025-11-18
**Bug Report:** `bug_report_2025-11-18_backward_yt_zero.md`
**Test:** `test_supersite_expansion_epochs`
**Status:** PARTIAL FIX - HMM probabilities now match, but sampling still diverges

## Root Cause

**Variable shadowing bug** in the backward pass loop at `phase_common/src/models/haplotype_segment_single.cpp:775`

### The Bug

The backward loop declared LOCAL variables `yt` and `nt`:
```cpp
for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
    ...
    float yt = 0.0f, nt = 1.0f;  // ← LOCAL variables shadow MEMBER variables!
```

These local variables shadowed the member variables declared in the header file (line 101):
```cpp
class haplotype_segment_single {
    ...
    float nt, yt;  // Member variables
```

### Why This Caused yt=0

1. The backward loop computed yt correctly in the LOCAL variable: `yt = M.getBackwardTransProb(...)`
2. But `SS_RUN_HOM()`, `SS_RUN_AMB()`, etc. accessed the MEMBER variable `this->yt`
3. The member variable was never updated during the backward pass, remaining at its initialized value of 0
4. This caused all supersite anchors to use yt=0 (no transition probability) in backward pass

### Why Forward Pass Worked

The forward pass (line 499) correctly used the member variables:
```cpp
yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(...);  // No local declaration
```

## The Fix

**File:** `phase_common/src/models/haplotype_segment_single.cpp`

### Change 1: Remove Local Variable Declaration

**Line 775 (REMOVED):**
```cpp
float yt = 0.0f, nt = 1.0f;
```

### Change 2: Initialize Member Variables Before Loop

**Lines 768-769 (ADDED):**
```cpp
yt = 0.0f;
nt = 1.0f;
```

### Complete Context

**Before:**
```cpp
// Flag: backward pass always starts with initialization
bool need_init = true;

for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
    curr_rel_locus = curr_abs_locus - locus_first;
    curr_rel_missing = curr_abs_missing - missing_first;
    const int prev_before = prev_abs_locus;
    char rare_allele = M.rare_allele[curr_abs_locus];
    bool update_prev_locus = true;
    float yt = 0.0f, nt = 1.0f;  // ← BUG: Local variables
```

**After:**
```cpp
// Flag: backward pass always starts with initialization
bool need_init = true;
yt = 0.0f;  // ← Initialize member variables
nt = 1.0f;

for (curr_abs_locus = locus_last ; curr_abs_locus >= locus_first ; curr_abs_locus--) {
    curr_rel_locus = curr_abs_locus - locus_first;
    curr_rel_missing = curr_abs_missing - missing_first;
    const int prev_before = prev_abs_locus;
    char rare_allele = M.rare_allele[curr_abs_locus];
    bool update_prev_locus = true;
    // ← No local yt/nt declaration, uses member variables
```

## Verification

### Before Fix

**Supersite locus 6 backward:**
```
[SS_RUN_HOM_TRACE] locus=6
  Input probSumT=80.4800109863
  Transition: yt=0.0000000000 nt=1.0000000000  ← BUG!
  Output probSumT=0.9940953255
```

**Biallelic locus 3 backward:**
```
[BIAL_RUN_HOM_TRACE] locus=3
  Input probSumT=80.4800109863
  Transition: yt=0.0152668785 nt=0.9847331047  ← Correct
  Output probSumT=0.9885177016
```

**Divergence:** 0.994095 / 0.988518 = 1.0056x (0.56% error)

### After Fix

**Supersite locus 6 backward:**
```
[SS_RUN_HOM_TRACE] locus=6
  Input probSumT=80.4800109863
  Transition: yt=0.0152668795 nt=0.9847331047  ✓ FIXED!
  Output probSumT=0.9885177016
```

**Biallelic locus 3 backward:**
```
[BIAL_RUN_HOM_TRACE] locus=3
  Input probSumT=80.4800109863
  Transition: yt=0.0152668785 nt=0.9847331047
  Output probSumT=0.9885177016
```

**Result:** ✓ IDENTICAL backward probabilities!

### Full Backward Probability Comparison (Burn3)

| Anchor Pair | Biallelic betaSumT | Supersite betaSumT (Before) | Supersite betaSumT (After) |
|-------------|-------------------|----------------------------|---------------------------|
| BI4 ≡ SS8 | 80.480011 | 80.480011 | 80.480011 ✓ |
| BI3 ≡ SS6 | 0.988518 | 0.994095 ❌ | 0.988518 ✓ FIXED |
| BI2 ≡ SS4 | 0.698952 | 0.702958 ❌ | 0.698952 ✓ FIXED |
| BI1 ≡ SS2 | 0.505000 | 0.505000 | 0.505000 ✓ |
| BI0 ≡ SS0 | 0.505000 | 0.505000 | 0.505000 ✓ |

## Remaining Issue

Despite the fix, the test still fails:
```
Anchor mismatch at locus 0: bial (0|1) vs supersite (0|0)
Anchor haplotypes diverged during burn3
```

### Analysis

- ✓ Forward pass probabilities: IDENTICAL
- ✓ Backward pass probabilities: IDENTICAL (after fix)
- ✓ Panel haplotype selection: IDENTICAL
- ❌ Final sampled haplotypes: DIFFERENT

### Hypothesis

The divergence is likely occurring in:
1. **Sampling step**: The HMM sampling algorithm may have non-deterministic behavior or numerical precision issues
2. **Random number generation**: Different sampling draws despite identical probabilities
3. **Another subtle bug**: Possibly in how probabilities are normalized or how sampling weights are computed

### Next Investigation Steps

1. Add logging to the sampling step to capture:
   - The exact posterior probabilities being sampled from
   - The random number draws used for sampling
   - The conditioning haplotype indices selected

2. Check if there are floating-point precision differences that affect sampling:
   - Compare probSumH values at higher precision (more decimal places)
   - Check if normalization produces slightly different weights

3. Verify the sampling algorithm is deterministic given identical inputs:
   - Check random seed initialization
   - Verify sampling uses the same probability arrays

## Impact

The variable shadowing bug caused:
- **0.56% divergence** in backward probabilities at each affected anchor
- **Compounding errors** through the backward pass
- **Different final haplotype sampling** due to accumulated probability differences

With this fix, the HMM forward and backward probabilities are now numerically identical between biallelic and supersite representations, eliminating the primary source of divergence.

## Files Modified

1. `phase_common/src/models/haplotype_segment_single.cpp`
   - Line 775: Removed local variable declaration
   - Lines 768-769: Added member variable initialization

## Diagnostic Code Added (for debugging)

Lines 810-826 in `haplotype_segment_single.cpp`: YT_DEBUG logging

This can be removed or kept for future debugging.

## Additional Fixes in This Session

1. **Enhanced trace logging** (4 locations in `haplotype_segment_single.h`):
   - Removed donor count limits (k < 4 or k < 10) to log all conditioning haplotypes
   - Lines 599, 724, 1137, 1329

These changes can be kept as they improve debugging capabilities without affecting performance (only active when SHAPEIT5_TEST_TRACE=1).

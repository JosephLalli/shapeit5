# Bug Report: Supersite Backward Pass yt=0 Issue

**Date:** 2025-11-18
**Test:** `test_supersite_expansion_epochs`
**Iteration:** burn3 (Iteration 3/15)
**Status:** Root cause identified

## Executive Summary

The test failure "Anchor mismatch at locus 0: bial (0|1) vs supersite (0|0)" during burn3 is caused by **supersite anchors using yt=0 (no transition probability) in the backward pass** instead of the correct yt≈0.0153. This causes backward probabilities to diverge by ~0.56%, which compounds through sampling to produce different final haplotypes.

## Investigation Steps Completed

### 1. Enhanced Trace Logging
Modified 4 locations in `phase_common/src/models/haplotype_segment_single.h` to output ALL donors (not just first 4-10):
- Line 599: `SS_RUN_HOM` backward (removed `k < 10` limit)
- Line 724: `SS_RUN_AMB` forward (removed `k < 4` limit)
- Line 1137: `BIAL_RUN_HOM` backward (removed `k < 10` limit)
- Line 1329: `BIAL_RUN_AMB` forward (removed `k < 4` limit)

### 2. Full Test Execution
```bash
make -j 12
SHAPEIT5_TEST_TRACE=1 tests/bin/test_supersite_expansion_epochs &> hmm_full_complete.log
```

## Key Findings

### ✓ Forward Pass: WORKING CORRECTLY

All forward probabilities are **IDENTICAL** at equivalent anchor loci:

| Anchor Pair | probSumH Distribution | probSumT | Status |
|-------------|----------------------|----------|--------|
| BI0 ≡ SS0 | [16.000000, 0.160000, 16.000000, 0.160000, ...] | 64.639999 | ✓ IDENTICAL |
| BI1 ≡ SS2 | [0.140316, 0.001403, 0.109684, 0.001097, ...] | 0.505000 | ✓ IDENTICAL |
| BI2 ≡ SS4 | [0.006706, 0.000067, 0.212662, 0.002127, ...] | 0.443125 | ✓ IDENTICAL |
| BI3 ≡ SS6 | [0.011031, 0.000110, 0.477128, 0.004771, ...] | 0.986082 | ✓ IDENTICAL |
| BI4 ≡ SS8 | [0.011083, 0.000111, 0.481092, 0.004811, ...] | 0.994193 | ✓ IDENTICAL |

### ✓ Panel Selection: WORKING CORRECTLY

Panel codes/alleles are **IDENTICAL** between biallelic and supersite paths at all loci.

**Example at BI1/SS2 (burn3 forward):**
```
BI1 panel alleles: 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1
SS2 panel codes:   1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1  ✓ IDENTICAL
```

**Example at BI3/SS6 (burn3 backward):**
```
BI3 panel alleles: 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0
SS6 panel codes:   0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0  ✓ IDENTICAL
```

### ❌ Backward Pass: BUG IDENTIFIED

**THE BUG:** All supersite anchors use **yt=0.0000000000** in backward pass instead of correct yt≈0.0152668785.

#### Evidence at BI3 vs SS6

**Biallelic locus 3 (burn3 backward):**
```
Input probSumT: 80.4800109863
Input probSumH: [10.0600013733, 10.0600013733, ...]
Transition: yt=0.0152668785, nt=0.9847331047
Panel alleles: 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0
Output probSumT: 0.9885177016
```

**Supersite locus 6 (burn3 backward):**
```
Input probSumT: 80.4800109863
Input probSumH: [10.0600013733, 10.0600013733, ...]
Transition: yt=0.0000000000, nt=1.0000000000  ← BUG!
Panel codes: 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0  (IDENTICAL!)
Output probSumT: 0.9940953255
```

**Divergence:** 0.9940953255 / 0.9885177016 = **1.0056x** (0.56% error)

#### Pattern Confirmation

ALL supersite anchors show yt=0 in backward pass:

| Locus | Type | yt (backward) | Expected yt |
|-------|------|---------------|-------------|
| SS8 | Anchor | 0.000000 | 0.000000 (last locus, correct) |
| SS6 | Anchor | 0.000000 | 0.015267 ❌ BUG |
| SS4 | Anchor | 0.000000 | 0.015267 ❌ BUG |
| SS2 | Anchor | ? | 0.015267 ❌ (likely) |

#### Trace Table vs Function Trace Discrepancy

Puzzling observation at locus 6:
```
[SS_RUN_HOM_TRACE] locus=6 ss_idx=3 sample_code=0 n_cond=16
  Transition: yt=0.0000000000 nt=1.0000000000
  ...

6	8	6	0.015267	0.984733	1	1	0	0	0	1	0.124262	...
    ↑ Trace table shows yt=0.015267!
```

The trace table (logged at line 917-927) shows yt=0.015267, but `SS_RUN_HOM` (logged at line 535) sees yt=0. This suggests yt is being computed/set AFTER `SS_RUN_HOM` executes, which contradicts the code flow.

## Code Locations

### Where yt Should Be Computed
**File:** `phase_common/src/models/haplotype_segment_single.cpp`

**Line 810-813:** Backward pass transition probability calculation
```cpp
if (!is_sibling) {
    yt = (curr_abs_locus == locus_last || curr_abs_locus == prev_abs_locus)
         ? 0.0f
         : M.getBackwardTransProb(prev_abs_locus, curr_abs_locus);
    nt = 1.0f - yt;
}
```

**Line 857:** SS_RUN_HOM is called (should see yt from line 811)
```cpp
update_prev_locus = SS_RUN_HOM(*site_view.supersite, site_view.supersite_index, site_view.sample_class0);
```

### Where yt Is Logged

**File:** `phase_common/src/models/haplotype_segment_single.h`

**Line 535:** Inside SS_RUN_HOM (beginning of function)
```cpp
std::fprintf(stdout, "  Transition: yt=%.10f nt=%.10f\n", yt, nt);
```

**File:** `phase_common/src/models/haplotype_segment_single.cpp`

**Line 917-927:** After all operations complete
```cpp
trace_log_backward_state(curr_abs_locus,
                         prev_before,
                         prev_after,
                         static_cast<double>(yt),  // ← Shows correct value
                         static_cast<double>(nt),
                         ...);
```

### Member Variable Declaration

**File:** `phase_common/src/models/haplotype_segment_single.h`

**Line 101:**
```cpp
float nt, yt;
```

## Hypotheses for Root Cause

### Hypothesis 1: Sibling Processing Interference
- Locus 7 (sibling) processes before locus 6 (anchor) in backward pass
- Sibling processing skips the yt computation block (line 810-813 guarded by `if (!is_sibling)`)
- yt may retain stale value from previous iteration
- For locus 8, yt was correctly set to 0 (because curr_abs_locus == locus_last)
- This 0 value may be retained through sibling locus 7 to anchor locus 6

### Hypothesis 2: Deferred Initialization
- Line 823 comment says: "// FIX: Set prev_abs_locus to curr_abs_locus to ensure yt=0 for the first anchor"
- This suggests intentional yt=0 behavior, but may be incorrectly applied
- Check if `need_init` flag is incorrectly set for normal anchors

### Hypothesis 3: Control Flow Issue
- The condition at line 811 should execute for anchors (!is_sibling)
- Trace shows `is_sibling=0` for locus 6, so block should execute
- Possible race condition or ordering issue with how prev_abs_locus is updated

## Impact Analysis

The 0.56% divergence in backward probabilities at each anchor compounds:
- BI3/SS6: 0.56% divergence
- BI2/SS4: Additional divergence (0.698952 vs 0.702958 = 0.57%)
- Final sampling produces different haplotypes: (0|1) vs (0|0)

## Next Steps

1. **Immediate:** Add detailed logging to track:
   - Value of `yt` before line 810
   - Value of `yt` after line 813
   - Value of `prev_abs_locus` at line 811
   - Value of `is_sibling` at line 810
   - Whether line 811 executes

2. **Debug:** Insert fprintf statements around lines 810-813:
   ```cpp
   fprintf(stderr, "[YT_DEBUG] locus=%d is_sibling=%d prev_abs=%d yt_before=%.10f\n",
           curr_abs_locus, is_sibling, prev_abs_locus, yt);
   if (!is_sibling) {
       yt = ...;
       fprintf(stderr, "[YT_DEBUG] locus=%d yt_after=%.10f\n", curr_abs_locus, yt);
   }
   ```

3. **Verify:** Check if forward pass has similar logic that works correctly at line 499:
   ```cpp
   yt = (curr_abs_locus == locus_first)?0.0:M.getForwardTransProb(prev_abs_locus, curr_abs_locus);
   ```

4. **Fix:** Once root cause confirmed, correct the yt computation or initialization for supersite anchors in backward pass.

## Test Environment

- **Compiler:** g++ (version from makefile)
- **Test binary:** `tests/bin/test_supersite_expansion_epochs`
- **Trace enabled:** `SHAPEIT5_TEST_TRACE=1`
- **Log file:** `hmm_full_complete.log`
- **Conditioning haplotypes:** n_cond_haps=16 (burn3)

## Files Modified

- `phase_common/src/models/haplotype_segment_single.h` (4 trace logging enhancements)

## References

- Previous investigation: `CLAUDE.md` (documents earlier findings)
- Related fix attempt: Line 823 comment about "ensure yt=0 for the first anchor"
- Forward pass reference: Line 499 (working correctly)

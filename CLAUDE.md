# SHAPEIT5 Supersite HMM Debugging Summary

## Problem Statement

Test `test_supersite_expansion_epochs` is failing with divergent results between biallelic and supersite (multiallelic) representations of the same genetic variants. The test runs 15 iterations of MCMC phasing and expects identical results from both representations.

## Diagnostic Commands

**Run test with full trace:**
```bash
SHAPEIT5_TEST_TRACE=1 tests/bin/test_supersite_expansion_epochs &> hmm_trace.log
```

**Check backward probabilities at burn3:**
```bash
awk '/Iteration 3\/15/,/Iteration 4\/15/' hmm_trace.log | grep "^[0-9]" | tail -20
```

**Check if yt is now correct:**
```bash
awk '/Iteration 3\/15/,/Iteration 4\/15/' hmm_trace.log | grep "Transition: yt=" | head -10
```

**Extract burn3 backward for both samples:**
```bash
awk '/Iteration 3\/15/,/Iteration 4\/15/' hmm_trace.log | awk '/# Backward.*bial/,/SET_FIRST_TRANS/' | grep "^[0-9]"
awk '/Iteration 3\/15/,/Iteration 4\/15/' hmm_trace.log | awk '/# Backward.*supersite/,/SET_FIRST_TRANS/' | grep "^[0-9]"
```

## Files Modified

1. **Enhanced debugging (can be kept):**
   - `phase_common/src/models/haplotype_segment_single.h` (lines 599, 724, 1137, 1329)
   - `phase_common/src/models/haplotype_segment_single.cpp` (lines 810-826)

2. **Documentation:**
   - `bug_report_2025-11-18_backward_yt_zero.md` - Initial investigation
   - `bug_fix_2025-11-18_variable_shadowing.md` - Fix documentation

## Test Configuration

- **Test:** `tests/bin/test_supersite_expansion_epochs`
- **PBWT depth:** 16 (but iteration 1 uses 32)
- **Iterations:** 15 (burn1, burn2, burn3, ...)
- **Failure point:** burn3 (iteration 3)
- **Trace enabled:** `SHAPEIT5_TEST_TRACE=1`

---

## Scenario 8 (8x Replication) Probability Redistribution Bug (2025-11-26)

### Problem Statement

Test scenario 8 (40 supersites, 80 total variants, 8 segments) fails at iteration 9 (burn7) with divergent ambiguous masks:

```
Iteration 8/15 [prune2]
[ITER_TRACE] amb_state post-prune2: bial(seg=2,amb0=0xaa) supersite(seg=2,amb0=0x6a)

Iteration 9/15 [burn7]
Anchor mismatch at locus 0: bial (1|0) vs supersite (0|1)
```

**Status:**
- ✅ Scenarios 1, 2, 4 (5, 10, 20 supersites) PASS
- ❌ Scenario 8 (40 supersites, 8 segments) FAILS

### Root Cause Investigation (2025-11-26)

#### Timeline of Divergence

Added comprehensive `BWD_PROB_DETAIL` logging to track backward emission probabilities (`prob[]`) at every locus:

```cpp
fprintf(stderr, "[BWD_PROB_DETAIL] sample=%s locus=%d bio_anchor=%d is_sib=%d seg=%d seg_locus=%d probSumT=%.15f\n",
        G->name.c_str(), curr_abs_locus, bio_anchor_idx, (int)is_sibling,
        curr_segment_index, curr_segment_locus, (double)probSumT);
fprintf(stderr, "  prob[0][0-7]= ");
for (int i = 0; i < std::min(8, (int)prob.size()); i++) {
    fprintf(stderr, "%.8e ", (double)prob[i]);
}
```

**Timeline Analysis:**

| Iteration | Phase | Segments | Divergence? | Details |
|-----------|-------|----------|-------------|---------|
| 5 | burn5 | 8 | ❌ NO | All prob[] values identical |
| 6 | prune1 | 8→3 | ❌ NO | After merge: `amb0=0xa6` for BOTH |
| 7 | burn6 | 3 | ✅ **YES** | **100x divergence first appears** |
| 8 | prune2 | 3→2 | ✅ YES | Masks diverge: `0xaa` vs `0x6a` |
| 9 | burn7 | 2 | ✅ YES | Test fails with haplotype mismatch |

**Critical Finding:** Divergence appears AFTER prune1 completes, during the first backward pass of burn6 (iteration 7). This indicates that **prune1's segment merging creates an inconsistent state** that manifests as divergent backward probabilities in scenario 8 (8 segments) but not in scenarios 1-4 (1-4 segments).

#### Divergence Details

**Location:** bio_anchor=36, prob[0] lane 2

**Backward emission probabilities at bio_anchor=36 during burn6:**
- Biallelic locus=36, seg=2, seg_locus=6: **probSumT=0.505000114440918**
  ```
  prob[0][0-7]= 3.90625000e-03 3.90624991e-05 3.90625000e-03 3.90624991e-05
                3.90624991e-05 3.90624991e-05 3.90625000e-03 3.90625000e-03
  ```

- Supersite locus=72, seg=2, seg_locus=12: **probSumT=0.505000114440918** (IDENTICAL!)
  ```
  prob[0][0-7]= 3.90625000e-03 3.90624991e-05 3.90624991e-05 3.90625000e-03
                3.90624991e-05 3.90625000e-03 3.90625000e-03 3.90624991e-05
  ```

**Analysis:**
- Total probability `probSumT` is **identical** (0.5050001144)
- But `prob[0][2]` differs by **100x**: biallelic=3.906e-03 vs supersite=3.906e-05
- Also `prob[0][3]` is swapped: biallelic=3.906e-05 vs supersite=3.906e-03
- The prob[] distribution pattern is **completely different** despite same total

**Pattern visualization:**
```
Biallelic:  [hi, lo, HI, lo, lo, lo, hi, hi]   ← bit 2 is HIGH
Supersite:  [hi, lo, LO, hi, lo, hi, hi, lo]   ← bit 2 is low, bit 3 is HIGH
```

### Root Cause Hypothesis

This is a **probability redistribution bug**, NOT the 2x accumulation bug seen in scenario 4. The characteristics:

1. **Total probability preserved:** probSumT is identical between biallelic and supersite
2. **Distribution incorrect:** Probability mass is redistributed to different haplotype lanes
3. **Triggered by pruning:** Appears after prune1 merges 8→3 segments
4. **Affects merged segments:** bio_anchor=36 is in segment 2 (a merged segment)
5. **Sibling-related:** Supersite has `seg_locus=12` (includes siblings) vs biallelic `seg_locus=6`

**Hypothesis:** When prune1 rebuilds ambiguous masks for merged segments using `performMerges()` (genotype_prune.cpp:274-305), the segment membership logic (`vrel < Lengths[s-1]`) correctly identifies which variants belong to which segment, BUT there may be an additional bug in how the **probabilities** are distributed across haplotype lanes when:
- The merged segment contains multiple supersites (each with siblings)
- Segment boundaries don't align with supersite anchor/sibling patterns
- The `Lengths[]` array (which includes siblings) is used inconsistently

### Segment Structure at bio_anchor=36

After prune1 merges 8→3 segments:
- **Segment 2** (where bio_anchor=36 resides):
  - Biallelic: `Lengths[2]` includes positions creating `seg_locus=6` for anchor 36
  - Supersite: `Lengths[2]` includes positions + siblings creating `seg_locus=12` for anchor 36
  - The 2x difference in `seg_locus` suggests siblings are counted in Lengths[]

**Critical observation:** Despite segment structure differences, the HMM should produce identical prob[] distributions for anchors. The fact that probSumT matches but individual lanes differ suggests a lane-assignment or indexing bug.

### Next Investigation Steps

1. **Trace prune1 merge operation for segment 2:**
   - Which original segments were merged to create segment 2?
   - How were Ambiguous masks rebuilt during the merge?
   - Check `performMerges()` logging for segment 2 in scenario 8

2. **Compare lane assignments before/after prune1:**
   - Check if conditioning haplotype indices are correct
   - Verify that merged haplotypes map correctly to lanes
   - Look for off-by-one errors in lane indexing

3. **Check sibling handling in merged segments:**
   - Are siblings correctly skipped in ambiguous mask rebuilding?
   - Does `vrel < Lengths[s-1]` work correctly when Lengths[] includes siblings?
   - Verify that `arel` (ambiguous anchor counter) increments correctly

4. **Test hypothesis: Force identical Lengths[]:**
   - Temporarily modify code to use `Lengths_bio[]` instead of `Lengths[]` in segment membership check
   - If this fixes the divergence, confirms that sibling counting in Lengths[] is the issue

### Files to Investigate

1. **genotype_prune.cpp:274-305** - `performMerges()` ambiguous mask rebuilding
   - The recent fix (2025-11-25) changed `arel < n_amb_first` to `vrel < Lengths[s-1]`
   - This fixed scenarios 1-4 but scenario 8 still fails
   - May need additional logic for multi-segment merges with siblings

2. **haplotype_segment_single.cpp:786-1100** - Backward pass with merged segments
   - Check if merged segments use correct Lengths[] values
   - Verify sibling handling doesn't affect prob[] distribution

3. **genotype_build.cpp:86-220** - Initial segment creation and Lengths[] calculation
   - Confirm Lengths[] includes siblings while Lengths_bio[] counts only anchors
   - Verify this is consistent across all segment operations

### Diagnostic Commands

**Extract prune1 merge details:**
```bash
SHAPEIT5_TEST_TRACE=1 timeout 30 tests/bin/test_supersite_expansion_epochs 2>&1 | \
  awk '/Scenario: repeat_factor=8/,0' | \
  awk '/Iteration 6\/15/,/Iteration 7\/15/' | \
  grep "PRUNE_DEBUG\|performMerges" > /tmp/prune1_merge.log
```

**Compare prob[] at bio_anchor=36 before/after prune1:**
```bash
for iter in 5 6 7; do
  echo "=== Iteration $iter ==="
  awk "/Iteration $iter\/15/,/Iteration $((iter+1))\/15/" /tmp/scenario8_bwd_detail2.log | \
    grep "BWD_PROB_DETAIL.*bio_anchor=36" -A1
done
```

### Status

- ✅ Identified that scenario 8 has a DIFFERENT bug than scenario 4 (100x redistribution vs 2x accumulation)
- ✅ Pinpointed divergence first appears in burn6 (after prune1)
- ✅ Confirmed probSumT matches but prob[] distribution differs
- ✅ Identified bio_anchor=36 in merged segment 2 as divergence point
- ❌ Root cause in prune1 segment merging not yet identified
- 🔍 Next: Investigate performMerges() lane assignment logic for multi-segment merges with siblings

---

## AMB Indexing Investigation - Hypothesis Disproven (2025-11-26)

### Investigation Goal

Verify whether siblings are incorrectly marked as AMB, causing the `a` counter in make() to increment incorrectly when indexing the Ambiguous[] array.

### Key Findings

**1. Siblings are correctly excluded from AMB counting**

`window_set.cpp:78-86` defines `is_locus_ambiguous()` which explicitly filters siblings:

```cpp
auto is_locus_ambiguous = [&](unsigned int locus) -> bool {
    if (ctx.is_member && !ctx.is_anchor) return false;  // SIBLINGS EXCLUDED
    return ctx.has_het || ctx.has_sca;
};
```

**2. Verified with WindowAmbLocus trace**

Scenario 1 supersite_sample shows siblings correctly excluded:

```
[WindowAmbLocus] locus=0 is_amb=1  ← anchor (het)  ✓
[WindowAmbLocus] locus=1 is_amb=0  ← sibling ✓ CORRECTLY EXCLUDED
[WindowAmbLocus] locus=2 is_amb=1  ← anchor (het)  ✓
[WindowAmbLocus] locus=3 is_amb=0  ← sibling ✓ CORRECTLY EXCLUDED
```

**3. genotype::build() correctly skips siblings** (genotype_build.cpp:318-319):

```cpp
if (ctx.is_member && !ctx.is_anchor) continue;  // Skip siblings when counting AMB
```

### Conclusion

**AMB indexing hypothesis DISPROVEN.** Siblings are correctly excluded throughout. The 100x divergence in scenario 8 must have a different root cause.

---

**Document Status:** 2025-11-26 - AMB indexing hypothesis disproven

---

## Corrected Timeline - Bug is in prune2, NOT burn6 (2025-11-26)

### Investigation Correction

Previous investigation in CLAUDE.md incorrectly identified burn6 as the divergence point. **New findings show the bug is in prune2 (iteration 8).**

### Verified Timeline for Scenario 8

```
burn1-5:  ✓ Haplotypes IDENTICAL (amb0=0xaa for both)
prune1:   ✓ Haplotypes IDENTICAL after merge (amb0=0xa6 for both)
burn6:    ✓ Haplotypes IDENTICAL (amb0=0xa6 for both)
prune2:   ❌ Haplotypes DIVERGE (bial=0xaa, supersite=0x6a)
burn7:    ❌ Test FAILS with anchor mismatch at locus 0
```

### Panel vs Test Sample Analysis

**Panel haplotypes after prune1:** IDENTICAL ✓
- Verified: First 10 panel samples match across all anchors
- Conclusion: Panel is not the source of divergence

**Test sample haplotypes:**
- Through burn6: IDENTICAL
- After prune2: DIVERGED

### Ambiguous Mask Divergence Detail

**Post-prune2 divergence:**
- Biallelic: `amb0 = 0xaa = 0b10101010`
- Supersite: `amb0 = 0x6a = 0b01101010`
- **Bit difference:** Bit 7 (MSB) - set in biallelic, clear in supersite

### Key Question

Why does prune1 work correctly but prune2 fails? Possible differences:
1. Different segment structures going into merge
2. Different number of segments being merged (prune1: 6→3, prune2: 3→2)
3. Edge case in performMerges() that only triggers in second pruning

### Next Investigation

Add detailed logging to prune2 to trace:
1. Input segment structure (lengths, amb counts)
2. Which segments are being merged
3. Ambiguous mask rebuilding step-by-step
4. Compare biallelic vs supersite at each merge step

---

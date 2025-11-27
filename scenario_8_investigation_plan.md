# Scenario 8 Investigation Plan

## Current Status

**Fixed:** Scenarios 1, 2, 4 (5, 10, 20 supersites) - All PASS ✅
**Failing:** Scenario 8 (40 supersites, 80 variants, 8 segments) - FAILS at iteration 9 (burn7) ❌

## Failure Details

```
Iteration 8/15 [prune2]
[ITER_TRACE] amb_state post-prune2: bial(seg=2,amb0=0xaa) supersite(seg=2,amb0=0x6a)

Iteration 9/15 [burn7]
Anchor mismatch at locus 0: bial (1|0) vs supersite (0|1)
```

### Divergence Point
- **When:** After prune2 (iteration 8), manifests in burn7 (iteration 9)
- **What:** Ambiguous masks diverge: `0xaa` (biallelic) vs `0x6a` (supersite)
- **Bit difference:** Bit 7 (MSB) - set in biallelic, clear in supersite
- **Impact:** Different haplotype selection → test failure

## Investigation Steps

### 1. Verify Fix Application in Scenario 8

Run scenario 8 with detailed trace to confirm the Lengths_bio fix is being applied:

```bash
SHAPEIT5_TEST_TRACE=1 timeout 60 tests/bin/test_supersite_expansion_epochs 2>&1 | \
  awk '/Scenario: repeat_factor=8/,0' > /tmp/scenario8_trace.log
```

Check:
- Are segment boundaries correct for 8-segment structure?
- Is `curr_segment_locus` incrementing only for anchors?
- Are siblings correctly excluded from segment locus counting?

### 2. Compare Segment Structures

Extract segment information for both representations:

```bash
grep -E "Lengths|n_segments|seg=" /tmp/scenario8_trace.log | head -50
```

Verify:
- Do both representations have the same number of segments?
- Do segment boundaries align on biological anchors?
- Are Lengths[] and Lengths_bio[] consistent?

### 3. Trace prune2 Merge Operation

The divergence happens after prune2. Investigate the merge operation:

```bash
SHAPEIT5_TEST_TRACE=1 timeout 60 tests/bin/test_supersite_expansion_epochs 2>&1 | \
  awk '/Scenario: repeat_factor=8/,0' | \
  awk '/Iteration 8\/15/,/Iteration 9\/15/' > /tmp/prune2_detail.log
```

Check:
- How many segments are merged? (likely 3→2 based on "seg=2" in output)
- Are ambiguous masks correctly rebuilt during merge?
- Does `performMerges()` handle siblings correctly?

### 4. Backward Probability Analysis

According to CLAUDE.md, the 100x divergence pattern suggests a probability redistribution bug:

```bash
# Extract backward probabilities at bio_anchor=36 (known divergence point from docs)
grep "BWD_PROB_DETAIL.*bio_anchor=36" /tmp/scenario8_trace.log -A1
```

Look for:
- Total probability (probSumT) - should be identical
- Individual lane probabilities (prob[0][0-7]) - may show redistribution
- Pattern: `[hi, lo, HI, lo, lo, lo, hi, hi]` vs `[hi, lo, LO, hi, lo, hi, hi, lo]`

### 5. Check for Additional Bugs

Scenario 8 is more complex (8 segments vs 1-4 for passing scenarios). Possible additional issues:

**A. Multi-segment window handling**
- Windows spanning multiple segments may have edge cases
- Check window boundaries: `[WindowAmbRange]` traces

**B. Segment indexing after multiple prunes**
- Prune1 (6→3 segments), Prune2 (3→2 segments)
- Each prune calls `performMerges()` - check for cumulative errors

**C. Transition probability buffer**
- More segments = more transitions
- Verify `curr_abs_transition` doesn't go negative (checked with asserts)

**D. Alpha storage indexing**
- `Alpha[rel_seg]` access pattern with merged segments
- Check if relative segment index calculation is correct

### 6. Hypothesis: Sibling Positioning in Merged Segments

From CLAUDE.md:
> bio_anchor=36 is in segment 2 (a merged segment)
> Biallelic: seg_locus=6
> Supersite: seg_locus=12 (includes siblings)

Even with Lengths_bio fix, there might be an issue where:
1. Segment 2 is created by merging multiple original segments
2. The merged segment has complex sibling patterns
3. `curr_segment_locus` tracking gets out of sync during merge

**Test this hypothesis:**
```bash
# Look for segment 2 structure after prune2
grep "seg=2" /tmp/prune2_detail.log | head -20
```

### 7. Double-Check the Fix

Review the changes made to ensure they're complete:

```bash
git diff HEAD~1 phase_common/src/models/haplotype_segment_single.cpp | grep -E "Lengths_bio|is_sibling"
```

Verify:
- All 8 locations use `Lengths_bio`
- All `curr_segment_locus++` have sibling checks
- All `curr_segment_locus--` have sibling checks
- No remaining `Lengths[...]` in segment boundary logic

### 8. Potential Missing Fix

Check if there are other places where `Lengths[]` is used for segment logic:

```bash
grep -n "Lengths\[" phase_common/src/models/haplotype_segment_single.cpp | \
  grep -v "Lengths_bio"
```

Look for uses of `Lengths[]` that should be `Lengths_bio[]` but were missed.

## Next Steps

1. Run step 1-3 to gather detailed traces
2. Analyze traces to identify divergence point
3. Compare with scenarios 1-4 behavior (which pass)
4. If issue is in `performMerges()`, check genotype_prune.cpp:274-305
5. If issue is in segment locus tracking, add more detailed logging
6. Consider adding scenario-specific assertions to catch the exact moment of divergence

## Success Criteria

- Scenario 8 completes all 15 iterations without mismatch
- Biallelic and supersite ambiguous masks remain identical through all iterations
- All backward probabilities match between representations

## Related Files

- `phase_common/src/models/haplotype_segment_single.cpp` - Forward/backward HMM
- `phase_common/src/objects/genotype/genotype_prune.cpp` - Segment merging logic
- `phase_common/src/objects/genotype/genotype_build.cpp` - Segment creation
- `CLAUDE.md` - Existing investigation notes
- `tests/src/test_supersite_expansion_epochs.cpp` - Test implementation

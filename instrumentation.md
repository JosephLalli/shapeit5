# Instrumentation and Debugging Guide

This document provides a comprehensive guide to debugging SHAPEIT5, particularly for diagnosing divergences between biallelic and supersite (multiallelic) HMM pathways.

## Environment Variables

### Primary Trace Controls

| Variable | Values | Default | Purpose |
|----------|--------|---------|---------|
| `SHAPEIT5_TEST_TRACE` | `0` (off) / `1` (on) | `0` | Enables most trace logging including HMM internals, backward/forward passes, emission calculations |
| `SHAPEIT5_ENABLE_OUTER_PRODUCT` | `0` (off) / `1` (on) | `0` | **Enables** outer product at segment boundaries (disabled by default due to divergence bugs) |
| `SHAPEIT5_SUPERSITE_GUARDS` | `0` (off) / `1` (on) | `1` | Enables supersite invariant guards (mutual exclusivity, h0/h1 vs bits, c0/c1 compatibility, packed code parity, SC normalization). |
| `SHAPEIT5_SUPERDEBUG_INVARIANTS` | `0` (off) / `1` (on) | `0` | When guards are on, emit concise invariant violations to stderr (init, projection, panel update). |
| `SHAPEIT5_SUPERDEBUG_SAMPLENAME` | string | unset | If set, filters invariant checks/logging to this sample. |
| `SHAPEIT5_SUPERDEBUG_BP` | integer bp | unset | If set, filters invariant checks/logging to supersite anchors at this bp. |
| `SHAPEIT5_SUPERSITE_GUARDS` | `0` (off) / `1` (on) | `1` | Enables supersite invariant guards (mutual exclusivity, h0/h1 vs bits, c0/c1 compatibility, packed code parity). |
| `SHAPEIT5_SUPERDEBUG_INVARIANTS` | `0` (off) / `1` (on) | `0` | When guards are on, emit concise invariant violations to stderr (init, projection, panel update). |
| `SHAPEIT5_SUPERDEBUG_SAMPLENAME` | string | unset | If set, invariant checks/logging are filtered to this sample. |
| `SHAPEIT5_SUPERDEBUG_BP` | integer bp | unset | If set, invariant checks/logging are filtered to supersite anchors at this bp. |

### Usage Example
```bash
# Basic trace logging
SHAPEIT5_TEST_TRACE=1 timeout 20 tests/bin/test_supersite_expansion_epochs &> trace.log

# Enable outer product for comparison (not recommended due to bugs)
SHAPEIT5_ENABLE_OUTER_PRODUCT=1 SHAPEIT5_TEST_TRACE=1 tests/bin/test_supersite_expansion_epochs
```

## Trace Log Tags

### HMM Backward Pass

| Tag | Location | Data | When Enabled |
|-----|----------|------|--------------|
| `[BWD_PROB_TRACE]` | `haplotype_segment_single.cpp:1001` | `locus`, `is_sib`, `probSumT`, `head_sum`, `prob[0]` | `SHAPEIT5_TEST_TRACE=1` |
| `[YT_DEBUG]` | `haplotype_segment_single.cpp:878-892` | `locus`, `is_sibling`, `is_anchor`, `prev_abs`, `yt` (transition prob) | `SHAPEIT5_TEST_TRACE=1` |
| `[SIBLING_SKIP_TRACE]` | `haplotype_segment_single.cpp:991` | `locus`, `prev_loc_in/out`, `kind`, `is_sib`, `update_prev`, `yt` | `SHAPEIT5_TEST_TRACE=1` |

### Segment Boundary Collapse Functions

| Tag | Location | Data | Purpose |
|-----|----------|------|---------|
| `[BIAL_COLLAPSE_HOM_ENTER]` | `haplotype_segment_single.h:1242` | `locus`, `probSumT_before`, `yt`, `nt`, `n_cond` | Entry to biallelic segment boundary |
| `[BIAL_COLLAPSE_HOM_FLAGS]` | `haplotype_segment_single.h:1259` | `use_outer`, `rel_prev_seg`, `nt_val`, `tFreq_val` | Biallelic collapse algorithm choice |
| `[BIAL_COLLAPSE_HOM_EXIT]` | `haplotype_segment_single.h:1282` | `locus`, `probSumT_after`, `probSumH[0-3]` | Exit from biallelic collapse |
| `[SS_COLLAPSE_HOM_ENTER]` | `haplotype_segment_single.h:797` | `locus`, `ss_idx`, `probSumT_before`, `yt`, `nt`, `n_cond` | Entry to supersite segment boundary |
| `[SS_COLLAPSE_HOM_FLAGS]` | `haplotype_segment_single.h:820` | `use_outer`, `rel_prev_seg`, `nt_val`, `tFreq_val` | Supersite collapse algorithm choice |
| `[SS_COLLAPSE_HOM_PROBSUMK]` | `haplotype_segment_single.h:826` | First 4 `probSumK[]` values | Per-haplotype sums from previous segment |
| `[SS_COLLAPSE_HOM_EXIT]` | `haplotype_segment_single.h:852` | `locus`, `probSumT_after`, `probSumH[0-3]` | Exit from supersite collapse |

### HMM Emission Functions

| Tag | Location | Data | Purpose |
|-----|----------|------|---------|
| `[BIAL_RUN_HOM_TRACE]` | Various | `locus`, `sample_allele`, `n_cond` | Biallelic homozygous emission |
| `[SS_RUN_HOM_TRACE]` | Various | `locus`, `ss_idx`, `sample_code`, `n_cond` | Supersite homozygous emission |
| `[BIAL_RUN_AMB_TRACE]` | Various | `locus`, `amb_code` | Biallelic ambiguous (heterozygous) emission |
| `[SS_RUN_AMB_TRACE]` | Various | `locus`, `ss_idx`, `amb_mask` | Supersite ambiguous emission |
| `[SupersiteEmit]` | `haplotype_segment_single.cpp:831` | `locus`, `emit`, `allele`/`anchor_class`, `lane_class` | Emission classification and lane assignments |
| `[BIAL_LANE_DEBUG]` | Various | `amb_code`, per-lane `g0`/`g1` emission values | Detailed lane emission probabilities |
| `[SS_LANE_DEBUG]` | Various | `amb_mask`, per-lane expected classes | Detailed supersite lane assignments |

### Transition Probability Storage

| Tag | Location | Data | Purpose |
|-----|----------|------|---------|
| `[TRANS_HAP]` | `haplotype_segment_single.h:~1587` | `locus`, `seg`, `prev_total`, `yt`, `nt`, `sumHProbs`, `prob[0-3]` | Haplotype transition probabilities at segment boundaries |
| `[SET_OTHER_TRANS]` | `haplotype_segment_single.h:~1309` | `locus`, `bio_anchor`, `seg`, `is_ss_anchor`, `sumHProbs`, `sumDProbs` | Stores transition probs for segment boundary |
| `[SET_FIRST_TRANS_DEBUG]` | Various | `use_outer`, `probSumT` | First segment initialization |

### Sampling and Path Selection

| Tag | Location | Data | Purpose |
|-----|----------|------|---------|
| `[SAMPLE_PATH]` | Various | `sample`, `direction`, `n_segments` | Start of sampling process |
| `[SAMPLE_FWD]` | Various | `sample`, `segment`, `sumProbs`, per-diplotype `prob`/`norm`/`cum` | Forward sampling details with cumulative probabilities |
| `[SAMPLE_FWD_PICK]` | Various | `sample`, `seg`, `selected_idx`, `dipcode`, `hap0`, `hap1` | Selected diplotype and haplotypes |
| `[SAMPLE_BWD_PICK]` | Various | Same as FWD_PICK | Backward sampling selection |
| `[SAMPLE_DEBUG]` | Various | `sample`, `direction`, `dip_mask`, `DipSampled0` | Sampling summary |
| `[MAKE_SS_AMB]` | Various | `locus`, `ss_idx`, `sampled_lanes`, `amb_code`, allele classes | Ambiguous mask construction for supersites |
| `[SupersiteProject]` | Various | `sample`, `ss_idx`, `locus`, `h0`, `h1` | Projection of supersite to constituent variants |

### Iteration and Test Progress

| Tag | Location | Data | Purpose |
|-----|----------|------|---------|
| `Iteration X/15 [name]` | Test output | Iteration number and phase name | Marks iteration boundaries (burn1-9, prune1-5, main1) |
| `Scenario: repeat_factor=X` | Test output | Scenario identifier | Marks test scenario boundaries |
| `[ITER_TRACE] amb_state` | Various | Iteration phase, segment, ambiguous mask values | Ambiguous mask state after each iteration |
| `Anchor mismatch` | Test output | `locus`, biallelic haplotype, supersite haplotype | Test failure with divergence location |

## Test Structure

### test_supersite_expansion_epochs

**Scenarios** (by repeat_factor):
- `repeat_factor=1`: 5 supersites (10 total variants with siblings)
- `repeat_factor=2`: 10 supersites (20 total variants)
- `repeat_factor=4`: 20 supersites (40 total variants) ← Most common failure point
- `repeat_factor=8`: 40 supersites (80 total variants, multi-segment)

**Iterations** (15 total):
1. `burn1` - burn-in iteration 1
2. `burn2` - burn-in iteration 2
3. `burn3` - burn-in iteration 3
4. `burn4` - burn-in iteration 4
5. `burn5` - burn-in iteration 5
6. `prune1` - first pruning iteration
7. `burn6` - burn-in after prune
8. `prune2` - second pruning
9. `burn7` - burn-in after prune
10. `prune3` - third pruning
11. `burn8` - burn-in after prune
12. `prune4` - fourth pruning
13. `burn9` - burn-in after prune
14. `prune5` - fifth pruning
15. `main1` - final main iteration

**Sample names**:
- `bial_sample` - uses biallelic representation (each variant separate)
- `supersite_sample` - uses supersite representation (multiallelic grouped)

**Locus numbering**:
- Biallelic: sequential (0, 1, 2, ..., N-1)
- Supersite: includes siblings (0, 1=sib, 2, 3=sib, 4, ..., 2N-1)
  - **Mapping**: biallelic locus X = supersite locus 2*X (anchor only, not siblings)

**Segment structure** (example for repeat_factor=4):
- 3 segments total (indices 0, 1, 2)
- Segment lengths differ between biallelic and supersite due to siblings
- Boundaries align at same biological anchors

## Common Analysis Workflows

### 1. Capture Full Trace
```bash
# Proper stderr redirection (use &> or put 2>&1 at END)
SHAPEIT5_TEST_TRACE=1 timeout 20 tests/bin/test_supersite_expansion_epochs &> /tmp/trace.log

# WRONG: 2>&1 > file  (redirects stderr to OLD stdout)
# RIGHT: &> file      (redirects both)
# RIGHT: > file 2>&1  (redirects stdout, then stderr to stdout)
```

### 2. Isolate Specific Scenario and Iteration
```bash
# Extract scenario 4, iteration 5
awk '/Scenario: repeat_factor=4/,0' /tmp/trace.log | awk '/Iteration 5\/15/,/Iteration 6\/15/' > /tmp/burn5.log

# Extract all of scenario 4
awk '/Scenario: repeat_factor=4/,/Scenario: repeat_factor=8/' /tmp/trace.log > /tmp/scenario4.log
```

### 3. Find First Backward Probability Divergence
```bash
# Extract BWD_PROB_TRACE for both samples
awk '/Iteration 5\/15/,/Iteration 6\/15/' /tmp/trace.log | grep "BWD_PROB_TRACE" > /tmp/bwd_traces.txt

# First ~20 lines are biallelic, next ~40 are supersite (with siblings)
# Compare probSumT values for matching biological anchors:
#   biallelic locus X ↔ supersite locus 2*X (is_sib=0)

# Example: check if locus 14 diverges
grep "locus=14.*is_sib=0" /tmp/bwd_traces.txt
```

### 4. Compare Segment Boundary Handling
```bash
# Get COLLAPSE traces for specific locus
grep "COLLAPSE_HOM.*locus=14" /tmp/trace.log

# Check use_outer flag differences
grep "COLLAPSE_HOM_FLAGS" /tmp/trace.log | grep "locus=14\|locus=28"

# Expected: biallelic L14 ↔ supersite L28 (same biological anchor 14)
```

### 5. Check Emission Classifications
```bash
# Find emission types for a locus
grep "SupersiteEmit.*locus=14" /tmp/trace.log

# Check lane assignments
grep "lane_class" /tmp/trace.log | grep "locus=14"
```

### 6. Track Iteration Progress and Divergence
```bash
# Find when divergence first occurs
grep "ITER_TRACE\|Anchor mismatch" /tmp/trace.log

# Compare ambiguous masks between samples
grep "ITER_TRACE.*amb_state" /tmp/trace.log | grep "prune1"
```

### 7. Side-by-Side Comparison
```bash
# Extract biallelic backward for iteration 5
awk '/# Backward trace - sample=bial_sample/,/# Backward trace - sample=supersite/' /tmp/trace.log | \
  grep "BWD_PROB_TRACE" > /tmp/bial_bwd.txt

# Extract supersite backward for iteration 5
awk '/# Backward trace - sample=supersite/,/SAMPLE_PATH/' /tmp/trace.log | \
  grep "BWD_PROB_TRACE" > /tmp/super_bwd.txt

# Compare with paste (align by biological anchor, accounting for siblings)
# Manual inspection needed due to locus number differences
```

### 8. Performance Profiling
```bash
# Count trace lines by type
grep -o "\[.*\]" /tmp/trace.log | sort | uniq -c | sort -rn

# Find most common trace tags
grep -oP '^\[\K[A-Z_]+' /tmp/trace.log | sort | uniq -c | sort -rn
```

## Key Debugging Concepts

### use_outer Flag
Controls segment boundary collapse algorithm:
- **`use_outer=0`**: Standard transition formula with `nt/probSumT` normalization
  - Formula: `prob = probSumK[k] * (nt/probSumT) + (yt/n_cond_haps)`
  - Used when outer product disabled or not available
- **`use_outer=1`**: Outer product mixing
  - Formula: `prob = col_mix * (row_stay * probSumK[k] + row_switch)`
  - Used when `supersites_enabled_flag=true` and conditions met

**Bug**: Biallelic has `supersites_enabled_flag=false` → always `use_outer=0`
Supersite has `supersites_enabled_flag=true` → can use `use_outer=1`
This causes divergence at segment boundaries!

**Default Fix**: Outer product is **disabled by default** (2025-11-26 fix). Both biallelic and supersite now use `use_outer=0` unless `SHAPEIT5_ENABLE_OUTER_PRODUCT=1` is explicitly set.

### Normalization Issues

When `use_outer=0`, the normalization factor `nt/probSumT` can **amplify** instead of normalize:
- If `probSumT < nt` (e.g., 0.505 < 0.985), then `nt/probSumT > 1` (e.g., 1.95)
- This multiplies probabilities by 1.95, causing `probSumT` to grow to 8.0 (HAP_NUMBER)
- Root cause: dividing by the **old** `probSumT` from previous locus

### Segment Boundaries vs Regular Loci

| Context | Function Called | Behavior |
|---------|----------------|----------|
| Regular locus (not at boundary) | `RUN_HOM`, `RUN_AMB`, `RUN_MIS` | Standard HMM update |
| Segment boundary (last locus) | `COLLAPSE_HOM`, `COLLAPSE_AMB`, `COLLAPSE_MIS` | Combines segments using transition probs or outer product |
| Sibling locus | `RUN_SIB`, `COLLAPSE_SIB`, or skipped | Bookkeeping only, no HMM math |

**Detection**: `pending_collapse = (curr_segment_locus == G->Lengths[curr_segment_index] - 1)`

### Supersite Anchors vs Siblings

- **Anchor**: First variant in a supersite group (is_sibling=0)
  - Represents the full supersite for HMM calculations
  - HMM probabilities computed here
- **Sibling**: Additional variants in supersite (is_sibling=1)
  - Inherit probabilities from anchor
  - Minimal processing (bookkeeping only)
  - Skip most HMM calculations with `continue` statements

**Identification**:
- Check `is_sibling` flag in traces
- Anchors have `locus == ss.global_site_id`
- Siblings have `locus != ss.global_site_id` but `ss_idx >= 0`

## Troubleshooting

### No trace output
- Check: `SHAPEIT5_TEST_TRACE=1` is set
- Check: Redirecting stderr correctly (`&>` or `2>&1`)
- Check: Trace function is gated by `supersite_trace_enabled()` or similar

### Can't find specific iteration
- Use: `grep -n "Iteration X/15"` to find line numbers
- Use: `awk '/Iteration X/,/Iteration Y/'` to extract range

### Locus numbers don't match between samples
- Remember: supersite locus = 2 × biallelic locus (for anchors)
- Filter supersite by `is_sib=0` to see only anchors
- Siblings alternate with anchors in supersite traces

### Trace file too large
- Use `timeout` to limit test duration
- Filter early: `grep "specific_tag" | head -1000`
- Process in chunks with `sed -n 'N,M p'`

## Quick Reference: Common Grep Patterns

```bash
# Find all segment boundary collapses
grep "COLLAPSE_HOM_ENTER"

# Find divergences in use_outer flag
grep "COLLAPSE_HOM_FLAGS"

# Compare probSumT at specific locus
grep "locus=14" | grep "probSumT"

# Check emission type
grep "SupersiteEmit.*locus=14"

# Find iteration failures
grep "Anchor mismatch\|assertion\|FAIL"

# Extract all backward traces
grep "BWD_PROB_TRACE"

# Filter by sample name
awk '/sample=bial_sample/,/sample=supersite/'

# Get transition storage events
grep "SET_OTHER_TRANS\|TRANS_HAP"
```

## Related Documentation

- `CLAUDE.md`: Historical debugging notes and bug fixes
- `gemini_theory.md`: Root cause analysis and hypotheses
- Test source: `tests/src/test_supersite_expansion_epochs.cpp`
- HMM implementation: `phase_common/src/models/haplotype_segment_single.{h,cpp}`

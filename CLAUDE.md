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
   - `.AGENT_markdowns/bug_report_2025-11-18_backward_yt_zero.md` - Initial investigation
   - `.AGENT_markdowns/bug_fix_2025-11-18_variable_shadowing.md` - Fix documentation

## Test Configuration

- **Test:** `tests/bin/test_supersite_expansion_epochs`
- **PBWT depth:** 16 (but iteration 1 uses 32)
- **Iterations:** 15 (burn1, burn2, burn3, ...)
- **Failure point:** burn3 (iteration 3)
- **Trace enabled:** `SHAPEIT5_TEST_TRACE=1`

---

## K-Inflation in Missing Data (2025-12-06)

### Problem Statement

Real-world integration tests on chr22:18-25mb region show massive accuracy degradation in supersite mode when missing data is present:

**Switch Error Rates:**
- Biallelic mode: 0.39% (0.05% pure switch)
- Supersite mode: 26.55% (7.5% pure switch)
- **68x worse performance**

**Key Finding:** Runs WITHOUT missing data do NOT show K-inflation. The bug is specific to missing data handling in supersite mode.

### Evidence

Sample-level comparison shows every sample has 50-60x more errors in supersite mode:
- HG00438: 0.41% → 24.30%
- HG01891: 0.69% → 29.89% (worst case)

Test data: `test/tmp/chr22.1KGP.18-25mb.phase_common.*.42.full_unphased.*`

### Root Cause Hypotheses

#### Hypothesis 1: Missing Index (`curr_abs_missing`) Mismatch (Most Likely)

**Location:** `haplotype_segment_single.cpp` lines 721-760

In supersite mode, sibling loci skip `AlphaMissing` caching:
```cpp
if (is_sibling && supersites_enabled_flag) {
    // Do NOT advance curr_abs_missing / curr_rel_missing
}
```

If the backward pass doesn't use identical skip logic, `AlphaMissing[idx]` will retrieve wrong forward probabilities, causing K-inflation through incorrect probability aggregation.

#### Hypothesis 2: `n_cond_haps` Semantic Mismatch

**Location:** `haplotype_segment_single.cpp` line 350

```cpp
n_cond_haps = idxH.size();  // Set once at construction
```

If conditioning set size differs between biallelic and supersite (due to different panel construction), all missing data normalization will be incorrect:
- RUN_MIS: `yt / (n_cond_haps * probSumT)`
- COLLAPSE_MIS: `yt / n_cond_haps`

#### Hypothesis 3: Supersite Missing Imputation Class Aggregation

**Location:** `haplotype_segment_single.cpp` lines 1172-1200

The backward pass aggregates Alpha*Beta across donor codes at missing supersites:
```cpp
for (int k = 0; k < n_cond_haps; ++k) {
    uint8_t dcode = get_donor_code_at_ss(ss_idx, k);
    for (int h = 0; h < HAP_NUMBER; ++h) {
        class_sum[h][dcode] += AlphaMissing[idx][k*HAP_NUMBER+h] * prob[k*HAP_NUMBER+h];
    }
}
```

If donor codes are incorrect or class aggregation is wrong, posteriors will be systematically biased.

### Debugging Instrumentation

Extensive tracing already in place via environment variables:

| Variable | Traces |
|----------|--------|
| `SHAPEIT5_TEST_TRACE=1` | Supersite emissions, probSumK, missing cache |
| `SHAPEIT5_SUPERDEBUG_SAMPLENAME=HG01891` | Focus on specific sample |
| `SHAPEIT5_SUPERDEBUG_BP=19000000` | Focus on genomic position |
| `SHAPEIT5_IMPUTE_DEEP=1` | Deep missing imputation traces |

**Key trace patterns:**
- `[FWD_MIS_TRACE]` - Forward missing cache (line 735)
- `[FWD_MIS_MAP]` - Missing index mapping (line 752)
- `[FWD_MIS_SKIP]` - Sibling skip logic (line 727)
- `[BWD_MIS_MAP]` - Backward retrieval (line 1154)
- `[SS_ALPHA_BETA_TRACE]` - Supersite Alpha*Beta (line 1184)
- `[IMPUTE_BIAL_TRACE]` - Biallelic imputation (line 1251)

### Investigation Commands

```bash
# Run biallelic with trace
SHAPEIT5_TEST_TRACE=1 \
./phase_common/bin/phase_common \
  --input 1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.filtered.bcf \
  --filter-maf 0.001 \
  --region chr22:18000000-25000000 \
  --map test/info/chr22.gmap.gz \
  --output test/tmp/debug_bial.bcf \
  --seed 42 --thread 1 \
  2>&1 | tee test/tmp/trace_bial_full.log

# Run supersite with trace
SHAPEIT5_TEST_TRACE=1 \
./phase_common/bin/phase_common \
  --input 1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.filtered.bcf \
  --filter-maf 0.001 \
  --region chr22:18000000-25000000 \
  --map test/info/chr22.gmap.gz \
  --output test/tmp/debug_super.bcf \
  --seed 42 --thread 1 --enable-supersites \
  2>&1 | tee test/tmp/trace_super_full.log

# Compare missing data handling
grep -E "FWD_MIS|BWD_MIS|n_missing" test/tmp/trace_bial_full.log > test/tmp/mis_bial.txt
grep -E "FWD_MIS|BWD_MIS|n_missing|SS_ALPHA_BETA" test/tmp/trace_super_full.log > test/tmp/mis_super.txt
diff test/tmp/mis_bial.txt test/tmp/mis_super.txt | head -100
```

### Critical Files

| File | Lines | Purpose |
|------|-------|---------|
| `phase_common/src/models/haplotype_segment_single.cpp` | 721-760 | Forward missing cache + sibling skip |
| `phase_common/src/models/haplotype_segment_single.cpp` | 1154-1291 | Backward missing retrieval |
| `phase_common/src/models/haplotype_segment_single.h` | 802-819 | SS_RUN_MIS |
| `phase_common/src/models/haplotype_segment_single.h` | 1030-1055 | SS_COLLAPSE_MIS |
| `phase_common/src/models/haplotype_segment_single.h` | 1629-1633 | SUMK() |

### Next Steps

1. Run full trace on 18-25mb region (both modes)
2. Compare missing index advancement between modes
3. Verify `n_cond_haps` is identical
4. Analyze Alpha*Beta aggregation at missing supersites
5. Identify first divergence point

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

## Supersite Imputation Produces Impossible Genotypes (2025-12-11)

### Problem Statement

When supersites are enabled, missing data imputation can produce **biologically impossible genotypes** where multiple mutually exclusive ALT alleles are assigned to the same haplotype at a multiallelic position.

### Evidence

**Test command:** `test/scripts/phase.chr22.wgs.sh 42 without_missing`

**Results:**
- Main algorithm: 9,142 incompatible sample-positions (0.29% of het positions at multiallelic sites)
- Supersite algorithm: 3,820 incompatible sample-positions (0.12%)

While supersites reduced incompatibility by 58%, **there should be exactly 0 incompatible haplotypes** with supersites enabled. Every remaining case is a bug.

### Example Case

**Position:** chr22:19010790
**Sample:** HG00351
**Status:** Missing genotype that was imputed

**Supersite output shows:**
```
A>G       GT=1|0
A>ATGTG   GT=1|0
ATG>A     GT=1|0
A>ATG     GT=0|0
ATGTGTG>A GT=0|1
ATGTG>A   GT=0|0
```

**Problem:** Haplotype 0 has THREE ALT alleles (G, ATGTG, and the deletion ATG>A). This is biologically impossible - each haplotype can only have one allele at any position.

### Scope of Problem

**Incompatibility breakdown by position complexity:**
- 2 variants: 219 incompatible cases
- 3 variants: 228 incompatible cases
- 4 variants: 291 incompatible cases
- 5 variants: 831 incompatible cases
- 6 variants: 2,251 incompatible cases (59% of all errors)

**Worst positions:** All have 5-6 variants
- pos=21637703: 101 incompatible samples, 6 variants
- pos=24196546: 97 incompatible samples, 6 variants
- pos=19164186: 92 incompatible samples, 5 variants

### Root Cause Hypothesis

The supersite imputation logic is not correctly constraining imputed genotypes to valid supersite classes. When a site is missing, the imputation should:

1. Compute posterior probability for each **supersite class** (not each biallelic variant independently)
2. Select a single class per haplotype
3. Expand that class back to biallelic genotypes

If step 1 or 2 treats the biallelic variants independently, impossible combinations can be imputed.

### Key Files to Investigate

| File | Lines | Purpose |
|------|-------|---------|
| `haplotype_segment_single.cpp` | 1172-1291 | Backward pass missing imputation |
| `haplotype_segment_single.h` | SS_RUN_MIS, SS_COLLAPSE_MIS | Supersite missing handling macros |
| `genotype_build.cpp` | | How imputed genotypes are written |

### Debugging Scripts

```bash
# Count incompatible haplotypes
python3 scripts/count_incompatible_haplotypes.py <phased.bcf> [region]

# Detailed analysis of incompatible sites
python3 scripts/debug_incompatible_sites.py <phased.bcf> [region] [max_examples]
```

### Next Steps

1. Trace the imputation path for sample HG00351 at position 19010790
2. Verify supersite class assignments during backward pass
3. Check if imputation is incorrectly treating biallelic loci independently
4. Ensure genotype writing respects supersite mutual exclusivity

---

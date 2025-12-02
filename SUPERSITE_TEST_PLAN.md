# Supersite Unit Test Gaps and Proposed Additions

This document summarizes where the current supersite test suite is strong, where gaps remain, and proposes a small set of focused tests to close those gaps.

The intent is **not** to add a huge new harness, but to introduce a few surgical tests that:

- Directly exercise multi-ALT (`ALT2`, `ALT3`, …) behavior in the HMM.
- Validate lane-level AMB broadcasting semantics at supersite anchors.
- Tie PBWT conditioning behavior more explicitly to multiallelic genotype similarity.
- Extend the existing epoch tests beyond the “dummy-alt (all 0/0 siblings)” scenario.

## 1. Existing Coverage (High Level)

Current supersite-related tests already cover:

- **Builder / Accessors / Encoding**
  - Grouping split records into `SuperSite` (`test_supersite_builder.cpp`, `test_supersite_vs_biallelic_simple.cpp`).
  - 4-bit panel code packing/unpacking and offsets (`test_supersite_unpack.cpp`, `test_supersite_accessor.cpp`).
  - Sample allele-code inference and HOM/AMB/MIS classification microcases (`test_supersite_accessor.cpp`, `test_supersite_class_microcases.cpp`).
- **Emissions and Class Equality**
  - Scalar and AVX2 precompute for supersite emissions, including non-REF/ALT1 codes (`test_supersite_emissions.cpp`, `test_supersite_emissions_real.cpp`).
  - Emission parity / match-count comparison against biallelic logic (`test_emission_logic_differences.cpp`).
- **PBWT and K-States**
  - PBWT conditioning parity between biallelic and supersite layouts (`test_pbwt_conditioning_parity.cpp`).
  - K-state sizes and contents tracking through multi-epoch MCMC (`test_supersite_expansion_epochs.cpp`).
- **Anchor Gating / Siblings / Boundary Behavior**
  - Only anchors run DP; siblings are true no-ops (`test_supersite_anchor_gating.cpp`, `test_supersite_boundary_conditions.cpp`).
  - Anchor encoding updated correctly for heterozygous supersites (`test_supersite_anchor_encoding.cpp`).
- **Phase 3 Multivariant Imputation**
  - SC buffer normalization and shape (`test_supersite_backward_sc.cpp`).
  - SC → projection → mutual exclusivity across splits (`test_supersite_backward_projection.cpp`, `test_missing_multiallelic_multinomial.cpp`).
- **Epoch Parity (Dummy-Alt Scenario)**
  - Biallelic vs supersite parity on a repeated micro dataset where extra splits are dummy 0/0 (`test_supersite_expansion_parity.cpp`, `test_supersite_expansion_epochs*.cpp`).

The gaps below are intentionally defined relative to this baseline.

## 2. Gap A – HMM Microcases for Real Multi-ALT Anchors

**Problem:** We do not currently have a single, analytic HMM microcase that covers:

- A supersite anchor where the sample is truly multiallelic (e.g. `ALT1|ALT2` or `ALT2|ALT3`).
- Per-lane AMB broadcasting semantics (`amb_code`) for those multi-ALT classes.
- Forward + backward behavior with a ground-truth probability we can compute by hand.

Most existing HMM tests either:

- Use dummy-alt siblings (superfluous splits with sample 0/0), or
- Focus on float/double parity / underflow safety rather than semantic correctness for multi-ALT anchors.

### Proposed Test: `test_supersite_hmm_multialt_microcase.cpp`

**Goal:** White-box verify that a single-window HMM with a truly multi-ALT supersite anchor produces the expected forward and backward probabilities, including lane-level emissions.

**Sketch:**

- **Variant map / supersite:**
  - Build a single supersite with 2 splits at the same `(chr,bp)` and `n_alts = 2` or `3`.
  - Example: `A→C` (ALT1) and `A→G` (ALT2) at one locus.
- **Conditioning panel:**
  - Use `K = 2` or `4` donors with simple, distinct codes:
    - hap0: class 0 (REF)
    - hap1: class 1 (ALT1)
    - hap2: class 2 (ALT2)
    - hap3: class 0 (REF)  (for redundancy)
- **Sample genotype:**
  - Case 1: `ALT1|ALT1` at the supersite (HOM).
  - Case 2: `ALT1|ALT2` at the supersite (AMB).
  - Case 3: `ALT2|ALT2` for contrast.
  - Construct the per-split bits so that `getSampleSuperSiteAlleleCode()` returns the correct `(c0, c1)` pair in each case.
- **HMM parameters and window:**
  - Use a single-segment, single-window setup:
    - `n_variants = supersite.var_count` but only `global_site_id` is a DP locus.
    - `t` and `nt` set so transitions are trivial (e.g. no transitions, or a single known `t`).
  - Choose `ed` and `ee` small enough that probabilities can be computed analytically.
- **Assertions:**
  - After `forward()`:
    - Check `prob[k*HAP_NUMBER + lane]` for each donor `k` and lane `h` matches the emission we expect from:
      - `lane_class[h]` derived from `(c0,c1)` and `amb_code`.
      - class equality mask: `emit = 1.0` if donor code == `lane_class[h]`, else `ed/ee`.
    - Check `probSumH[h]` and `probSumT` equals the analytically computed sum over donors and lanes.
  - After `backward()`:
    - With `t=0` (no transitions), the posteriors at the anchor should match the normalized emissions (forward-only case).
    - Optionally assert SC entries (if we call the multivariant imputation path) against a hand-computed `P(class_c | data, hap)`.

**Value:** This directly tests the “correct matching logic for non REF/ALT1 genotypes” and the lane-wise AMB broadcasting semantics at a supersite anchor, independent of PBWT or window segmentation.

## 3. Gap B – Lane-Level AMB Broadcasting Parity (Bial vs Supersite)

**Problem:** Existing emission parity tests focus on **match counts** across donors and lanes but do not explicitly verify that, for a given heterozygous multi-ALT pattern:

- Each lane’s “wants c0 vs c1” pattern (from `amb_code`) is preserved, and
- Supersite anchors behave like biallelic anchors when both are representing the same underlying diploid site.

### Proposed Test: `test_supersite_amb_lane_parity.cpp`

**Goal:** Verify that the 8-lane AMB semantics at supersite anchors mirror the biallelic AMB semantics when both are encoding the same underlying biological site.

**Sketch:**

- Build a tiny **paired dataset**:
  - Bial path: a single biallelic variant with REF/ALT.
  - Supersite path: the same site represented as a 2-split supersite where only one split carries an ALT (the other is dummy REF).
  - Ensure the sample genotype and conditioning panel yield an identical effective “two-class” problem.
- For both datasets:
  - Build `haplotype_segment_single` objects.
  - Run `forward()` to the AMB site.
- White-box inspect:
  - For each donor `k` and lane `h`, compare the per-lane emissions:
    - Extract biallelic `prob` entries immediately after the AMB emission.
    - Extract supersite `prob` entries at the corresponding anchor.
  - Assert they are identical (within tolerance) across all `k, h`.

**Value:** This test reuses the existing AVX2 emission logic but explicitly checks the lane semantics, not just aggregate match counts.

## 4. Gap C – PBWT Selection Driven by Multi-ALT Similarity

**Problem:** `test_pbwt_conditioning_parity.cpp` shows that supersite vs bial layouts produce *identical* conditioning sets when the underlying variation is identical. It does **not** assert that PBWT actually reflects multi-ALT similarity relative to a focal sample in a situation where different allele classes matter.

Because PBWT still operates on the split-bit biallelic matrix, the key questions are:

- Are we correctly masking siblings (`applySupersiteAnchorMask` / anchor redirect) so we don’t double-count multi-ALT sites?
- Does the PBWT neighbor set for anchors **prefer donors carrying the same allele class** (ALT1 vs ALT2 vs REF) in synthetic scenarios where this should be strongly favored?

### Proposed Test: `test_pbwt_multiallelic_similarity.cpp`

**Goal:** Construct a small synthetic panel where donors differ only at one multiallelic supersite, and verify that PBWT conditioning sets at that anchor are enriched for donors with the same allele code as the target haplotype.

**Sketch:**

- Build a supersite dataset:
  - One supersite with 3 splits (ALT1, ALT2, ALT3), plus a few flanking biallelic variants to give PBWT some context.
  - Panel haplotypes:
    - A block of donors carrying ALT1 only.
    - A block with ALT2 only.
    - A block with ALT3 only.
    - A block REF-only.
- Build a genotype set with one sample:
  - Case A: hap0 carries ALT1, hap1 REF at the supersite.
  - Case B: hap0 carries ALT2, hap1 REF.
- For each case:
  - Run PBWT selection (`conditioning_set::initialize` + `select`) with supersites enabled and proper anchor masking / redirect.
  - Extract the neighbor set at the supersite anchor locus using the same helper as in `test_pbwt_conditioning_parity.cpp`.
- Assertions:
  - At the anchor, donors carrying the same class as the sample’s haplotype should dominate the neighbor set (e.g. >80% of donors in K come from the matching block, given a strong block structure).
  - Switching from Case A (ALT1) to Case B (ALT2) should flip which donor block is enriched.

**Value:** This validates that PBWT selection “reflects the observed genotypes” in multiallelic contexts, not just that bial vs supersite layouts are identical.

## 5. Gap D – Epoch Harness with Non-Dummy Multi-ALT Anchors

**Problem:** The current 15-epoch tests (`test_supersite_expansion_epochs*.cpp`) intentionally use dummy 0/0 siblings and simple bial-like anchors. This is perfect for parity, but it does not exercise:

- Repeated sampling of multi-ALT anchors over burn-in / prune / main.
- The interaction between Phase 3 SC-based imputation, segment pruning, and supersite projection in a scenario where the sample is actually multiallelic.

### Proposed Test: `test_supersite_expansion_epochs_multialt.cpp`

**Goal:** Mirror the existing epoch harness, but with genuine multiallelic anchors, and assert stability of key invariants over epochs.

**Sketch:**

- Start from the existing `MiniContext` machinery in `test_supersite_expansion_epochs.cpp`.
- Modify the supersite scenario so that:
  - Some anchors are HET between distinct ALT classes (e.g. `ALT1|ALT2`, `ALT2|ALT3`).
  - The panel patterns ensure that each ALT class has a distinct block of donor haplotypes.
- Run a shorter schedule (e.g. 5 burn-in + 2 prune + 5 main) for tractability.
- At each epoch:
  - Check mutual exclusivity at all supersite anchors.
  - Track the per-anchor `(c0,c1)` class codes inferred via `getSampleSuperSiteAlleleCode()` and ensure they are consistent with the projected hap bits.
  - Optionally, check that the empirical distribution of sampled classes at an anchor reflects the underlying donor composition (e.g. ALT1 more frequent than ALT2 when there are more ALT1 donors).

**Value:** Extends the “epochs tests” from the dummy-alt world into a realistic multiallelic setting, validating that Phase 3 behavior is stable across the full MCMC schedule.

## 6. Gap E – Corner-Case Stress for Multi-ALT + Missing + Low Transitions

**Problem:** We have targeted tests for:

- Underflow behavior (`test_supersite_transition_underflow.cpp`).
- Fully missing supersites (`test_supersite_backward_sc.cpp`, `test_missing_multiallelic_multinomial.cpp`).

But we do not explicitly test the combination:

- Multi-ALT supersites,
- Many donors (`K` reasonably large),
- Very small `yt` (rare transitions),
- Mixtures of missing and observed supersites in the same window.

### Proposed Test: `test_supersite_multialt_underflow_stress.cpp`

**Goal:** Ensure that, in a window with alternating observed and missing multi-ALT supersites under low transition rates, the HMM:

- avoids numerical underflow (or gracefully escalates to double precision), and
- maintains well-normalized SC distributions at missing anchors.

**Sketch:**

- Build a variant map with 4–6 supersites in a single window:
  - Pattern: observed → missing → observed → missing.
  - Each supersite has 3–4 ALTs with different panel frequencies.
- Use `K ≈ 64` donors, with `yt` very small and `nt` near 1.
- Mark every second supersite as fully missing in the sample; others are multi-ALT HET.
- Run `forward()` and `backward()` in `haplotype_segment_single`, letting the implementation fall back to double if needed.
- Assertions:
  - No negative underflow codes returned.
  - For each missing supersite, SC per-haplotype sums to 1 within tolerance.
  - `probSumT` remains finite and non-zero across the window.

**Value:** Stress-tests the numerical regime where supersites are dense, multi-ALT, and transitions are rare, which is close to the worst-case for accumulation error.

## 7. Suggested Implementation Order

To keep iteration manageable:

1. **Start with Gap A** – `test_supersite_hmm_multialt_microcase.cpp`  
   Gives a precise, analytic anchor for multi-ALT semantics and will immediately catch any emission or classification mistakes.
2. **Add Gap B** – `test_supersite_amb_lane_parity.cpp`  
   Reuses the same micro harness but ties supersite AMB behavior directly to the well-tested bial HMM.
3. **Add Gap C** – `test_pbwt_multiallelic_similarity.cpp`  
   Extends PBWT coverage from pure parity to phenotype-guided similarity.
4. **Extend epochs** – Gap D (`test_supersite_expansion_epochs_multialt.cpp`)  
   Once the microcases are solid, drive them through a realistic epoch schedule.
5. **Optional stress test** – Gap E (`test_supersite_multialt_underflow_stress.cpp`)  
   Only if we observe numerical fragility in WGS-scale logs or want extra confidence.

These tests collectively address the original concerns:

- **PBWT selection that reflects observed genotypes** → Gap C.  
- **Accurate broadcasting** → Gaps A and B.  
- **Correct matching logic for non REF/ALT1 genotypes** → Gaps A, C, D, and E.


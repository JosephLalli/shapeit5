# Supersite Instrumentation and Unit Test Plan

This document is a concrete, step‑by‑step plan for adding:

- Targeted **instrumentation** to catch supersite‑related state corruption in real runs.
- Focused **unit tests** to enforce those invariants on small synthetic datasets.

The goal is to explain *exactly* what to add, where to add it, and how to run it, without needing to rediscover the design from previous notes.

The plan is organized in roughly the order we should implement it.

---

## 0. Goals and Scope

**Primary goals**

- Detect and localize **structural supersite bugs** that:
  - Corrupt haplotype bits (`G`, `H_opt_hap`, `H_opt_var`) around supersites.
  - Break the mapping between supersite class space (`c0/c1`, `SC`, `h0/h1`) and the actual split‑record bits.
  - Propagate errors to non‑supersite loci (as seen in `chr22.wgs.sh` plots).
- Make supersite invariants **observable and testable**:
  - At runtime in WGS tests (via env‑guarded instrumentation).
  - In unit tests on tiny, analytically understandable scenarios.

**Out of scope**

- Fixes to any specific bug; this plan only adds diagnostics and tests.
- Changes to the high‑level supersite algorithm (we assume the intended design is correct and want to enforce it).

---

## 1. New / Extended Environment Controls (status)

We will reuse existing flags where possible and add a small number of new ones:

1. **Existing flags to reuse**
   - `SHAPEIT5_SUPERSITE_GUARDS` (already described in AGENTS):
     - Default `=1`. Implemented to cover:
       - Mutual exclusivity (hap bits across siblings).
       - Hap bits vs sampled `h0/h1`.
       - Hap bits vs immutable `c0/c1`.
       - Packed donor code parity at supersite build.
       - SC row finiteness/normalization.
   - `SHAPEIT5_SUPERDEBUG_SAMPLENAME`:
     - Implemented: filters invariant checks/logging to this sample.
   - `SHAPEIT5_SUPERDEBUG_BP`:
     - Implemented: filters invariant checks/logging to supersite anchors at this bp.

2. **New flag (if needed)**
   - `SHAPEIT5_SUPERDEBUG_INVARIANTS` (optional; default `0`):
     - Implemented: when `1`, emit concise invariant violations to stderr (init, projection, panel update), still gated by sample/BP filters.

**Action items**

- [ ] Document these flags in `instrumentation.md` (linking to this plan).
- [ ] Make sure they are read once (e.g. in a small `SupersiteDebugConfig` helper) and cached, not parsed repeatedly in hot loops.

---

## 2. Core Helper: Supersite Consistency Checker (implemented)

### 2.1. Helper location and API

**Where**

- Add a small header/source pair within `phase_common`, e.g.:
  - `phase_common/supersite_debug.h`
  - `phase_common/supersite_debug.cpp`
- Keep this module *header‑only for hot helpers* if needed, but putting implementations in `.cpp` is fine because usage will be sparse and/or guarded.

**Core API**

Implemented in `phase_common/src/objects/supersite_debug.{h,cpp}`:
- `SupersiteDebugConfig::from_env()` reads guards/verbose/sample/BP filters; `enabled_for_sample()` applies filters.
- `SupersiteInvariantViolation` carries ss_idx/global_site_id/var range/message.
- `class_from_hap_bits(...)` and `check_supersite_consistency_for_sample(...)` enforce:
  - Mutual exclusivity (at most one ALT per hap across siblings).
  - Hap bits vs sampled `h0/h1`.
  - Hap bits vs immutable `c0/c1` snapshot.
  - Optional violation reporting; no abort unless caller asserts.

### 2.2. “Class from bits” utility

Add an internal utility to recover a supersite class from actual hap bits:

```cpp
// Returns 0=REF, 1..n_alts for ALT_i, or 255 for invalid (“multi-ALT” or inconsistent).
uint8_t class_from_hap_bits(
    const genotype& sample_g,
    const SuperSite& ss,
    const std::vector<int>& super_site_var_index,
    int hap_index   // 0 or 1
);
```

**Intended behavior**

- Scan all member variants:
  - Let `alt_for_hap = 0` initially.
  - For each `ai in [0, ss.n_alts)`:
    - `v_abs = super_site_var_index[ss.var_start + ai]`.
    - If hap has ALT at `v_abs`:
      - If `alt_for_hap == 0`: set `alt_for_hap = ai + 1`.
      - Else: return 255 (more than one ALT across siblings → invalid).
- Return:
  - `0` if `alt_for_hap == 0` and no REF/ALT inconsistency is observed.
  - `alt_for_hap` if exactly one ALT is set across siblings.
  - `255` for impossible patterns (e.g. incompatible with supersite semantics).

### 2.3. Invariants to check

Inside `check_supersite_consistency`, enforce the following for the *debug sample(s)* and supersites in scope:

1. **Mutual exclusivity across siblings**
   - For each supersite `ss` and each hap (`h=0,1`):
     - `class_from_hap_bits(...)` must be in `{0..n_alts}` (not `255`).

2. **Anchor vs sampled class (`h0/h1`)**
   - For anchors (where `locus == ss.global_site_id`):
     - Retrieve the stored sampled `h0/h1` classes for the current epoch (we may need a small accessor for this).
     - Require: `class_from_hap_bits(h) == hX` (where X is 0 or 1).

3. **`c0/c1` compatibility**
   - For non‑missing anchors:
     - `c0` and `c1` should be immutable snapshot codes.
     - Check:
       - `class_from_hap_bits(0) ∈ {0, c0, c1}`
       - `class_from_hap_bits(1) ∈ {0, c0, c1}`
   - For sites where both haplotypes are fully missing:
     - It is sufficient to check that no ALT bits are spuriously set in the splits.

4. **Sibling alignment**
   - For all siblings (variants sharing a supersite id but not the anchor locus):
     - Confirm that hap bits are consistent with the anchor class:
       - If `hX == 0`: all siblings must be REF for hap X.
       - If `hX == ai+1`: exactly one sibling, the `ai` index, is ALT for hap X; all others REF.

5. **Panel consistency (optional, but recommended)**
   - On the donor side (`H_opt_hap`):
     - Run `class_from_hap_bits` like logic over panel haplotypes to verify that:
       - The packed donor codes (`packed_allele_codes`) match the decoded bits.

### 2.4. Failure handling

When an invariant fails:

- If `guards_enabled`:
  - For tests: prefer to `assert(false)` or `throw std::runtime_error` with a concise message.
  - For production binaries: consider logging once per run and either `abort()` or continue depending on severity (we can iterate later).
- If `verbose`:
  - Print a concise line with:
    - `sample`, `chr`, `bp`, `ss_idx`, `global_site_id`
    - `class_from_hap_bits` values
    - Expected `h0/h1` and `c0/c1`

---

## 3. Instrumentation Hook Points in the Pipeline

This section describes **where** to call `check_supersite_consistency` (and lightweight variants) to catch corruption early.

### 3.1. After supersite build / metadata snapshot

**File(s):**

- `phase_common/super_site_builder.*`
- `phaser_initialise.cpp` (or wherever supersite metadata is attached to `genotype_set`).

**Hook**

- After `buildSuperSites()` and after `c0/c1` snapshot:
  - For debug sample(s) and relevant supersites, call a **lightweight** consistency check that only validates:
    - `class_from_hap_bits` is compatible with `c0/c1`.
    - No multi‑ALT patterns in initial hap bits.

**Reason**

- Ensures the snapshot itself is coherent before any HMM work or projection.

### 3.2. After `make(...)` and `projectSupersites()` in sampling

**File(s):**

- `genotype_management.cpp` (or wherever `make()` and `projectSupersites()` live).

**Hook**

- Immediately after `projectSupersites()` completes for a sample in an iteration:
  - Call `check_supersite_consistency(...)` for that sample (guarded).

**Reason**

- This is the **prime suspect** region: projection from sampled classes to bits.
- Catch bugs before haplotypes are written back to the panel and spread in PBWT.

### 3.3. After `H.updateHaplotypes(...)` and `H.transposeHaplotypes_H2V()`

**File(s):**

- `conditioning_set_update.cpp` or wherever these methods are implemented.

**Hook**

- After the panel update for a given sample:
  - Re‑check the supersite invariants, but now against `H_opt_hap` / `H_opt_var`.

**Reason**

- Confirms that panel storage faithfully reflects the sample’s supersite state and that no packing/unpacking bug is corrupting donor codes.

### 3.4. HMM anchor emission debug (forward / backward)

**File(s):**

- `haplotype_segment_single.cpp` / `.h`

**Hook (partial)**

- Supersite anchor emissions now log lane classes and donor codes under `SHAPEIT5_TEST_TRACE`:
  - `SS_RUN_AMB_CLASSES` (lane_class per lane) and `SS_RUN_AMB_CODES` (first donor codes).
  - `SS_RUN_HOM_CLASSES` logs per-lane class for HOM anchors.
- Remaining: add class-equality mask logging if needed.

**Reason**

- Verifies that emissions use the correct classes and that lane semantics match the `Ambiguous` mask at anchors.

### 3.5. Multivariant SC and sampling debug

**File(s):**

- Supersite‑aware backward pass, e.g. `haplotype_segment_single.cpp` in `IMPUTE_SUPERSITE_MULTIVARIANT`.
- Sampling code that turns SC into `h0/h1`.

**Hook (implemented for SC; sampling log pending)**

- SC rows are guard-checked for finiteness/normalization in both precisions.
- For missing supersite anchors, SC rows are logged per hap when tracing.
- Logging of the sampled class choices is still pending.

**Reason**

- Ensures multivariant imputation is numerically sensible and semantically consistent with projection.

---

## 4. Bial vs Supersite Equivalence Harness

The error plots strongly suggest that supersite runs are corrupting state in ways that propagate beyond supersite loci. We want a **small, deterministic harness** to compare:

- Run A: biallelic path (supersites disabled).
- Run B: supersite path (supersites enabled).

### 4.1. Harness design

**New test file**

- `tests/src/test_supersite_bial_equivalence_window.cpp`

**Dataset**

- Tiny synthetic VCF/BCF region (single chromosome, ~50–100 variants) with:
  - 2–3 supersites (anchors + siblings), each with 2–3 ALTs.
  - Flanking pure biallelic variants.
  - A few samples, including one “debug sample” with non‑trivial supersite genotypes (e.g. `0/1`, `0/2`, `1/2`).

### 4.2. Execution pattern

1. **Build two configurations**
   - Config A: `--enable-supersites=0`.
   - Config B: `--enable-supersites=1`.
   - Ensure identical:
     - Random seed.
     - Iteration counts (e.g. 1–2 burn‑in only).
     - Threading (single‑threaded for deterministic tests).

2. **Run phase_common once in each mode**
   - Use existing test harness machinery (many tests already run small in‑memory contexts).
   - Capture final haplotypes for all samples and loci into simple arrays.

### 4.3. Assertions

- At **non‑supersite loci**:
  - Hap bits must be **bit‑for‑bit identical** between runs A and B.

- At supersite loci:
  - For run B (supersites), reconstruct split bits by:
    - Reading the projected child variants’ hap bits.
    - Confirming that when you “forget” supersites and look only at split bial records, the sample’s hap bits match run A.

- If desired, allow for trivial label permutation (e.g. hap0/hap1 swap) but otherwise require exact equality.

---

## 5. Unit Tests for Supersite Invariants

We now connect the instrumentation to automated tests, so new regressions are caught early.

### 5.1. Test: Supersite projection invariants

**File**

- `tests/src/test_supersite_projection_invariants.cpp`

**Scenario**

- Construct a small genotype set with:
  - One chromosome.
  - One supersite group with 2 or 3 splits at the same `(chr,bp)`.
  - A few flanking bial variants (optional).
  - One or two samples with genotypes spanning:
    - `0/0`, `0/1`, `0/2`, `1/1`, `1/2`, `2/2` (for 2–3 ALTs).

**Procedure**

1. Build `SuperSite` metadata and snapshot `c0/c1`.
2. Run a single iteration of `phase_common` in a minimal mode where:
   - HMM runs, sampling occurs, `make()` and `projectSupersites()` are invoked.
3. With `SHAPEIT5_SUPERSITE_GUARDS=1`:
   - Call `check_supersite_consistency(...)` at the points described in Section 3.

**Assertions**

- The test should complete without any invariant violation or assertion.
- Optionally, perform explicit checks:
  - For each supersite and hap, `class_from_hap_bits` matches the expected class from the known genotype.

### 5.2. Test: Negative cases (intentional invariant breaks)

**File**

- Extend `test_supersite_projection_invariants.cpp` or create an additional `test_supersite_projection_invariants_negative.cpp`.

**Scenario**

- Use the same synthetic setup but manually corrupt:
  - Set two siblings ALT for the same hap.
  - Set a sibling ALT incompatible with `c0/c1`.

**Procedure**

- Call `check_supersite_consistency(...)` directly in the test and assert that it:
  - Returns `false` and populates `SupersiteInvariantViolation`, or
  - Triggers a controlled assertion / exception that the test expects.

---

## 6. Optional: PBWT Neighbor Tracking Around Supersites

To connect structural supersite errors to the observed WGS spikes, we may add a focused PBWT neighbor test.

### 6.1. Test: PBWT similarity reflects supersite classes

**File**

- `tests/src/test_pbwt_multiallelic_similarity.cpp`

**Scenario**

- Synthetic panel with:
  - One supersite with 3 ALT classes.
  - Donor haplotypes grouped into blocks by class (ALT1 block, ALT2 block, REF block).
  - One target sample with class `ALT1` at the anchor in one sub‑scenario, `ALT2` in another.

**Procedure**

1. Run PBWT selection with supersites enabled.
2. Extract neighbor indices at the supersite anchor for the target hap.

**Assertions**

- The K neighbors for the target hap should be dominated by donors whose class matches the target.
- Changing the target from ALT1 to ALT2 should flip which donor block dominates.

---

## 7. Running and Using the New Instrumentation

Once the above is implemented:

1. **Unit tests**
   - Build and run:
     - `make -C tests`
     - `ctest -R supersite_projection_invariants`
     - `ctest -R supersite_bial_equivalence_window`

2. **Real‑data WGS regression (e.g. `chr22.wgs.sh`)**
   - For a specific problematic region:
     - Set:
       - `SHAPEIT5_SUPERSITE_GUARDS=1`
       - `SHAPEIT5_SUPERDEBUG_SAMPLENAME=<sample>`
       - `SHAPEIT5_SUPERDEBUG_BP=<anchor_bp>`
       - Optionally `SHAPEIT5_SUPERDEBUG_INVARIANTS=1` for verbose logs.
   - Run `chr22.wgs.sh` and inspect logs for:
     - Supersite invariant violations.
     - Emission/SC summaries at the anchor.

This closes the loop between **real‑data spikes** and **small unit tests**: any structural bug that corrupts haplotype state around supersites should now manifest as a clear invariant violation in both environments.

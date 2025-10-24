# SHAPEIT5 — Multiallelic Site Constraint Enforcement (Handover)

This document explains, at a practical engineering level, what SHAPEIT5 does, why multiallelic sites are tricky, the three remedies we implemented (in rising complexity), how they’re wired into the codebase, and the main challenges we hit. It is written for a strong C++ developer without a genetics background.

## 1) What SHAPEIT5 currently does

- Phases diploid genotypes into haplotypes using:
  - PBWT (Positional Burrows–Wheeler Transform) to select a small panel of donor haplotypes per window.
  - A Li–Stephens HMM to sample each sample’s haplotypes across windows (burn‑in, pruning, main iterations).
- Two pipelines:
  - Common variants: `phase_common` (array/WGS)
  - Rare variants on a scaffold: `phase_rare`
- Variant representation is biallelic (REF vs ALT). Real multiallelic sites are split into several biallelic rows with identical (chr, pos).

## 2) Why multiallelic handling is a problem

Splitting a multiallelic site into K biallelic rows causes the HMM to treat rows independently. It can assign multiple ALT rows at the same position to the same haplotype in a sample (biologically impossible — a haplotype can carry at most one allele at a base). This can mislead downstream tools (imputation, burden tests, STR analysis).

## 3) Three solutions (in order of increasing complexity)

All three run after an HMM sampling step — they are small, local consistency fixes that do not alter the core HMM.

- Transition‑only post‑hoc flipping (simplest)
  - Use nearest phased heterozygous anchors p (left) and n (right) within the same phase set. For a violation (e.g., two ALT rows on hap0), evaluate the two ways to flip one row to hap1 using Li–Stephens stay/switch probabilities on [p→site] and [site→n]; choose the higher likelihood.
  - Pros: O(1) per violation, minimally invasive.
  - Limits: Does not consider per‑site emissions; does not use the actual donor identities; with >2 rows or missing anchors, may not fully resolve.

- Micro re‑decode (donor‑agnostic)
  - Enumerate all valid compound assignments for the K split rows (≤ (K+1)^2, pruned to ≤1 ALT per haplotype). Score each assignment with the run’s transition model (and a simple emission proxy or GL/PL if present). Apply the ML assignment and rewrite bits.
  - Pros: Handles arbitrary K rows; can fix both allele identity and haplotype assignment together.
  - Limits: Still donor‑agnostic; emissions in `phase_common` are approximate without GL/PL.

- Micro re‑decode with frozen donors (micro‑donor, most integrated)
  - Use the donor haplotypes (PBWT Kstates) selected for the sample/window. For each candidate, compute a small donor‑weighted chain score across p→site→n to reward candidates supported by donors.
  - Pros: Better ties resolution and consistency with the model actually used in the pass.
  - Limits: Current code implements a donor‑weighted chain (emissions + stay/switch) rather than a full donor‑conditioned DP; requires window/donor plumbing.

## 4) Current status in `phase_common` and `phase_rare`

### phase_common

- Feature flags:
  - `--enforce-oneallele` (enable) and `--oneallele-stats <file>`
  - `--oneallele-mode {transition|micro|micro-donor}` (select algorithm)

- Where it runs:
  - Per sample, inside worker loop after each sampling step: `phase_common/src/phaser/phaser_algorithm.cpp`
  - Final “belt‑and‑suspenders” pass after `G.solve()`: `phase_common/src/phaser/phaser_finalise.cpp`

- Key modules and functions:
  - `phase_common/src/modules/multiallelic_position_map.{h,cpp}`: build groups of split rows by (chr, pos).
  - `phase_common/src/modules/transition_scorer.{h,cpp}`: stay/switch scoring helpers.
  - `phase_common/src/modules/oneallele_enforcer.{h,cpp}`:
    - `set_enabled()`, `set_mode()`, `stats()` (API)
    - `enforce_sample(...)`: per‑sample entry point in worker; dispatches to the selected mode per position group.
    - `enforce(...)`: batch fallback (used at finalize step)
    - `enforce_group_transition(...)`: transition‑only fix at a position
    - `enforce_group_micro(...)`: micro/micro‑donor enumeration across arbitrary K rows
    - `evaluate_candidate_micro(...)`: scores a single compound assignment
    - `compute_chain_score_with_donors(...)`: donor‑weighted chain score
    - Helpers: `find_left_neighbor(...)`, `find_right_neighbor(...)`, `flip_alt_to_other_hap(...)`, `rebalance_excess_alts(...)`

### phase_rare

- Feature flags:
  - `--enforce-oneallele-rare`, `--oneallele-rare-stats`
- Where it runs:
  - After merging rare data back into the scaffold (step 3.5): `phase_rare/src/phaser/phaser_algorithm.cpp`
- Approach:
  - Rare pipeline stores per‑row phase probabilities (PP). It resolves violations by flipping the lowest‑confidence ALT rows on the offending haplotype until each hap has ≤1 ALT. This is conceptually similar to a simplified micro decision but tailored to the sparse data structure.

## 5) Challenges encountered

- Access to donors/windows outside the worker: we solved this by running enforcement per sample inside the worker loop.
- Anchor detection: respecting phase‑set/segment boundaries is critical; anchors outside the segment are ignored.
- Emissions in `phase_common`: without GL/PL we rely on a simple error model for micro scoring; accuracy improves if GL/PL are present.
- Determinism/ties: we use deterministic tie‑breaking to keep runs reproducible.
- Performance: transition‑only is essentially free; micro/micro‑donor remain cheap with small K and m (typically ≤8).

## Code map (quick reference)

- Common pipeline
  - `phase_common/src/phaser/phaser_algorithm.cpp`: per‑sample enforcement call sites
  - `phase_common/src/phaser/phaser_finalise.cpp`: final pass and stats output
  - `phase_common/src/phaser/phaser_parameters.cpp`: CLI flags and mode selection
  - `phase_common/src/modules/*`: implementation of grouping, scoring and enforcement

- Rare pipeline
  - `phase_rare/src/phaser/phaser_algorithm.cpp`: rare enforcement (step 3.5)
  - `phase_rare/src/phaser/phaser_parameters.cpp`: rare flags

## Running and validating

- Typical command (common):
  - `phase_common/bin/phase_common --input X.bcf --map GMAP.gz --region 1 --output out.bcf --enforce-oneallele --oneallele-mode micro`
- Checker: `tools/check_oneallele out.bcf` reports per‑haplotype violation counts per position.
- See `test/scripts` for integration exercises, including a multiallelic‑focused script that runs all three modes.


## WGS PAR2 Fixtures (1KGP T2T)

The WGS tests use CHM13 T2T 1KGP PAR2 inputs on chrX, with known multiallelic violations to exercise one-allele enforcement.

### Region and Map
- Region: `chrX:153929053-154248138` (PAR2)
- Map: `test/info/par2.gmap.gz`

### Inputs (BCF)
- Family: `test/wgs/target.family.1kgp_t2t.par2.bcf` (≈1803 samples)
- Unrelated: `test/wgs/target.unrelated.1kgp_t2t.par2.bcf` (≈1399 samples)
- Haploid cohort: `test/wgs/target.haploid.1kgp_t2t.par2.bcf` (≈1599 male samples)

These BCFs carry only `FORMAT:GT` (no `FORMAT:PL`). This is sufficient for phasing workflows used in tests.

### Sample Lists and Pedigree
- Family union from PED (Son, Father, Mother; drop NA/-1): `test/wgs/samples.family.1kgp_t2t.par2.txt`
- All samples in master BCF: `test/wgs/samples.all.1kgp_t2t.par2.txt`
- Unrelated complement: `test/wgs/samples.unrelated.1kgp_t2t.par2.txt`
- Haploid sample list (first column only): `test/wgs/samples.haploid.1kgp_t2t.par2.txt`
- Pedigree TSV used by scripts: `test/info/1kgp_t2t.par2.ped`

### Haploid Handling
- The haploid BCF remains diploid in representation (GT); haploid samples are declared via `--haploids info/1kgp_t2t.par2.haploid.txt` so heterozygous GTs are reset to missing within the pipeline.
- We do not require `bcftools +fixploidy` in the final setup; if ever needed, `test/wgs/par2.ploidy.txt` and `test/wgs/samples.haploid.sex.txt` illustrate ploidy control.

### Scripts
- `test/scripts/phase.wgs.family.sh`: uses family BCF + PED, region and map above
- `test/scripts/phase.wgs.unrelated.sh`: uses unrelated BCF
- `test/scripts/phase.wgs.haploid.sh`: uses haploid BCF + haploid list

All three use `--map test/info/par2.gmap.gz`, `--region chrX:153929053-154248138`, `--seed 15052011`, and `--thread 1`.

### Expectations
- Ground-truth outputs live in `test/scripts/expected/` as headerless `.vcf.gz` with corresponding `.md5` files:
  - `phase.wgs.family.{vcf.gz,md5}`
  - `phase.wgs.unrelated.{vcf.gz,md5}`
  - `phase.wgs.haploid.{vcf.gz,md5}`

### Violation Spot-Checks
- The repository includes `test/find_one_alt_per_hap_violation.py` to confirm split multiallelic positions with same-haplotype ALT duplications are present in the inputs. Running it on the three BCFs reports multiple violations in each subset, ensuring the constraint logic is exercised.

## Array Fixtures

Array tests remain unchanged functionally but now also pin `--seed 15052011` and `--thread 1`. Expected outputs are maintained as headerless `.vcf.gz` + `.md5` pairs in `test/scripts/expected/`.

## Determinism & Reproducibility

- Always run tests with `--thread 1` and a fixed `--seed` to avoid haplotype-orientation drift.
- Validate with `assert_same_md5` (not raw VCF diffs) to avoid spurious header/metadata mismatches.
- Expected files should be regenerated using known-good binaries and then committed.

This architecture ensures reliable regression testing and eliminates false failures from metadata differences, providing a solid foundation for SHAPEIT5 development.

## Handover: Where We Left Off and Next Steps

This section captures the precise development state in the codebase and the next engineering steps to complete multiallelic one-allele enforcement as described above.

### Current Implementation Snapshot

- Common pipeline (phase_common)
  - Implemented
    - Grouping of split-biallelic rows by position: `phase_common/src/modules/multiallelic_position_map.{h,cpp}` (built at init in `phase_common/src/phaser/phaser_initialise.cpp:137`).
    - Transition-only enforcement with left/right heterozygous anchors per phase segment: `phase_common/src/modules/oneallele_enforcer.{h,cpp}` using `transition_scorer.{h,cpp}`.
    - Enforcement runs once per MCMC iteration after sampling, before updating haplotypes: `phase_common/src/phaser/phaser_algorithm.cpp:139`.
    - Flags/stats: `--enforce-oneallele` and `--oneallele-stats` parsed and reported at finalize.
  - Not yet implemented
    - Mode selector `--oneallele-mode {transition|micro|micro-donor}` (docs mention it, code does not yet parse/use it).
    - Micro and micro-donor enumeration/scoring paths (`enforce_group_micro`, donor-weighted scoring, etc.).
    - Per-sample enforcement inside the worker loop to access donor Kstates; current enforcement runs outside workers post-iteration.
    - Final “belt-and-suspenders” enforcement pass after `G.solve()` just before writing output (finalize currently only prints stats).

- Rare pipeline (phase_rare)
  - Implemented
    - Enforcement after merging rare variants onto the scaffold (step 3.5) that flips the lowest-PP ALT on an offending hap until each hap has ≤1 ALT: `phase_rare/src/phaser/phaser_algorithm.cpp`.
    - Flags `--enforce-oneallele-rare` and `--oneallele-rare-stats` with finalize reporting.
  - Not yet implemented
    - Micro-style enumeration for rare sites (optional improvement); current approach is PP-based flipping only.

- Tests and tools
  - Unit: `test/unit/test_oneallele_enforcer.cpp` exercises transition scoring and resolving a basic violation.
  - Integration: WGS PAR2 and array scripts under `test/scripts`, plus `tools/check_oneallele` verifier.
  - Note: `test/scripts/phase.oneallele.array.family.sh` is missing a trailing line-continuation after `--enforce-oneallele`, which likely breaks the next `--seed` line; add a `\` if needed.

### Next Steps (Concrete Plan)

1) Add mode selection (common)
- Parse `--oneallele-mode` with default `transition` in `phase_common/src/phaser/phaser_parameters.cpp`.
- Extend `OneAlleleEnforcer` API with `set_mode()` and branch by mode.

2) Implement micro re-decode (donor-agnostic)
- Add an enumeration path in `oneallele_enforcer.{h,cpp}` that, for each multiallelic position and sample, explores all valid assignments across K split rows (≤ (K+1)^2, ≤1 ALT per hap), scores with transition terms (and GL/PL if present), uses deterministic tie-breaking, and applies the ML assignment.
- Keep O(K) memory and short-circuit when already valid.

3) Implement micro-donor scoring
- Plumb donor/window context from worker to enforcer. Move enforcement into the per-sample worker loop (right after `sample(...)` in `phase_common/src/phaser/phaser_algorithm.cpp::phaseWindow`) so Kstates are available.
- Add a donor-weighted chain score for p→site→n that rewards candidates supported by donors (conditioning size `m` already available via `set_conditioning_size`).

4) Final pass at finalize (common)
- After `G.solve()` and before writing output in `phase_common/src/phaser/phaser_finalise.cpp`, run a last `enforce(...)` call as a belt-and-suspenders sweep.

5) Stats and knobs
- Expand `OneAlleleStats` with per-mode counters; optionally expose `--oneallele-min-cm` to control the minimum genetic distance used by transition scoring (`set_min_distance_cm`).

6) Rare pipeline optional upgrade
- Replace the lowest-PP flipping loop with the same micro enumeration used in common (when K>2 or ties make flipping ambiguous), retaining current flags and stats.

7) Tests and fixtures
- Add unit tests for micro and micro-donor paths (K>2, missing anchors, ties).
- Add integration runs that exercise all three modes on the PAR2 fixtures (`--thread 1`, fixed `--seed`).
- Regenerate `.md5` expectations with known-good binaries after landing the new modes.

This encapsulates the practical “where we left off” and an ordered plan to complete multiallelic constraint enforcement, aligned with the architecture described above.

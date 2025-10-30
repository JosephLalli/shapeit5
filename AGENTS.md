## Project
SHAPEIT5 — `phase_common` super-site integration

## Source of Truth
- This AGENTS.md mirrors `.github/copilot-instructions.md`, which is the authoritative, detailed guide distilled from `SUPERSITE_CONVERSATION_SUMMARY.md`.
- If anything diverges, prefer `.github/copilot-instructions.md` and update this file to match.

## Objectives (unchanged)
- Expose optional super-site phasing for multi-allelic loci (`--enable-supersites`) while keeping default behaviour untouched.
- Maintain the original biallelic code paths (and performance characteristics) when supersites are disabled.

## Current Status (Oct 29, 2025)
- **Phase 3 Multinomial Imputation: COMPLETE** ✅
  - Native multinomial imputation for missing supersites implemented end-to-end
  - Guarantees mutual exclusivity by design (exactly one ALT per haplotype across all splits)
  - Clean compilation, ready for integration testing
- Data model / builder: `buildSuperSites` partitions split records, fills `super_sites`, `locus_to_super_idx`, `super_site_var_index`, and packs panel haplotype codes. Conditioning haplotypes are encoded with 4-bit classes (two per byte) as planned.
- HMM integration: `haplotype_segment_{single,double}` accept optional supersite state and compute multinomial posteriors in backward pass via `IMPUTE_SUPERSITE_MULTINOMIAL()`.
- Genotype projection: `genotype::make()` samples from multinomials and projects to split records, enforcing mutual exclusivity.
- Unit test coverage:
  - `test_supersite_emissions` (float/double kernels)
  - `test_supersite_accessor` (per-hap sample codes)
  - `test_supersite_unpack` (4-bit decode)
  - `test_supersite_builder` (grouping/panel packing)
  - `test_supersite_hmm` (forward path regression)
  - `test_supersite_hmm_states` (INIT/RUN/COLLAPSE micro-harness)
- Build / tooling: tests have a standalone `tests/makefile` mirroring top-level flags and linking supersite objects.

## Recent Work (Phase 3 Implementation)
**Problem Identified:** Split multi-allelic records were being imputed independently, allowing multiple ALTs on the same haplotype (e.g., both split1=1|0 and split2=1|0 at chr:pos, implying hap0 carries both ALT1 and ALT2 simultaneously). This violates biological reality.

**Solution Implemented:** Phase 3 native multinomial imputation treats supersites as first-class multi-allelic loci throughout the HMM pipeline:

1. **HMM Backward Pass** (`haplotype_segment_{single,double}.{h,cpp}`):
   - Added `IMPUTE_SUPERSITE_MULTINOMIAL()` function that computes P(class_c | Alpha, Beta) for all classes c ∈ {REF, ALT1, ..., ALTn}
   - Uses C SIMD accumulators (C ≤ 16) following same Alpha×Beta pattern as biallelic IMPUTE
   - Stores multinomial posteriors in SC buffer: `SC[offset + hap*C + c]`
   - Updated `backward()` signature to accept `SC` and `anchor_has_missing` pointers
   - Added conditional logic to call multinomial imputation at anchors, skip siblings

2. **Genotype Projection** (`genotype_header.h`, `genotype_managment.cpp`):
   - Added supersite context pointers to genotype class (super_sites, locus_to_super_idx, etc.)
   - Implemented `setSuperSiteContext()` to pass HMM multinomial posteriors
   - Updated `make()` to sample one class per haplotype from multinomial distribution
   - Projects sampled class to split records: if class=ALTi, set split_i=ALT, all others=REF
   - **Mathematical guarantee:** Exactly one split set to ALT per haplotype (mutual exclusivity)

3. **Integration** (`phaser_algorithm.cpp`, `compute_job.{h,cpp}`):
   - `compute_job` allocates SC buffer and populates `anchor_has_missing` per window
   - `phaser_algorithm` passes SC to backward() and calls `setSuperSiteContext()` before sample()
   - Works across all iteration types (BURN, PRUN, MAIN) and precision paths (float/double)

**Key Design Choice:** Phase 3 was chosen over Phase 1 (reconstruction) because:
- Phase 1: Compute biallelic posteriors, then reconstruct multinomial in genotype::make()
  - Requires complex normalization, potential for numerical instability
  - Post-hoc constraints may not perfectly enforce mutual exclusivity
- Phase 3: Compute multinomial natively in HMM
  - Single source of truth for posteriors
  - Sampling from multinomial mathematically guarantees mutual exclusivity
  - Cleaner separation of concerns (HMM computes probabilities, genotype samples)

## Outstanding Issues
1. **Integration testing needed**
   - Run `test/scripts/phase.wgs.unrelated.sh --enable-supersites` to verify correctness
   - Validate SC buffer contents (multinomials sum to ~1.0)
   - Check output BCF for mutual exclusivity violations (should be zero)
   - Compare phasing accuracy vs. biallelic mode
2. Unit test harness
   - `tests/bin/test_supersite_hmm_states` may need update for Phase 3 validation
   - Fixtures under `tests/data/` still needed for smoke tests (placeholders exist)
3. Build stability
   - `make -C tests` succeeds after setting `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib`. Ensure CI/exported instructions include this.

## Next Steps
1. **Integration testing** (priority 1): Run full phasing pipeline with --enable-supersites on test data
2. Add validation checks:
   - Assert Σ_c SC[offset+hap*C+c] ≈ 1.0 in IMPUTE_SUPERSITE_MULTINOMIAL
   - Log sampled classes in genotype::make() for debugging
   - Add mutual exclusivity validator in output writer
3. Profile performance: measure overhead of multinomial computation vs. biallelic
4. Author integration smoke test harness (`tests/run_supersite_smoke.sh`) with known multi-allelic VCF fixtures
5. Optional optimizations:
   - SIMD stores in IMPUTE_SUPERSITE_MULTINOMIAL (currently scalar extraction)
   - Separate SuperProbMissing storage if needed (currently SC serves this purpose)

## Architecture & Data Flow (Essentials)
- HMM cores: `phase_common/src/models/haplotype_segment_{single,double}` implement forward/backward with AVX2 on 8-lane vectors (HAP_NUMBER=8) using `aligned_vector32<T>` arrays for Alpha/Beta and sums.
- States: INIT/RUN/COLLAPSE per genotype class (HOM/AMB/MIS), mirrored between float/double implementations.
- Supersite check: Always guard supersite logic with `if (super_sites && locus_to_super_idx)` to preserve biallelic defaults.
- Precision fallback: Forward/backward return non-zero on underflow; caller (`phaser_algorithm.cpp`) retries with double precision.

## Supersite Data Model
- Built by `phase_common/src/objects/super_site_builder.{h,cpp}` via `buildSuperSites()` once per iteration after PBWT selection.
- 4‑bit allele codes (0=REF, 1..15=ALT) per conditioning haplotype; two codes per byte in a packed array; `panel_offset` is the starting byte for a supersite.
- `locus_to_super_idx[v]` maps variant index to supersite index or -1 when not in a supersite.
- Member variants are referenced via flat `super_site_var_index` with `var_start/var_count` (CSR-style). Large sites chunk when `n_alts>15`.
- **Phase 3 additions:**
  - `SuperSite.class_prob_offset`: byte offset into SC buffer for this supersite's multinomial posteriors
  - `SuperSite.n_classes`: number of classes (C = 1 + n_alts), cached for convenience
  - `compute_job.SC`: CurrentSuperClassPosteriors buffer storing P(class_c | hap) per missing supersite
  - `compute_job.anchor_has_missing`: boolean vector indicating which supersites have all members missing for this sample
- Runtime lookups:
  - Panel donor code: `unpackSuperSiteCode(panel_codes, ss.panel_offset, hap_idx)`
  - Sample code: `getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap)`
  - Multinomial posterior: `SC[class_prob_offset + hap*C + c]` = P(class_c | hap)
- Anchor gating: Run DP only at `SuperSite.global_site_id`; sibling split records skip emission/transition to avoid double counting.

## Supersite Coding Patterns (Mirror Biallelic)
- Function taxonomy: Extract `SS_INIT_*`, `SS_RUN_*`, `SS_COLLAPSE_*` with signatures/returns matching biallelic functions; dispatch from existing entry points after centralized classification.
- Centralized classification:
  - Compute `c0`/`c1` via `getSampleSuperSiteAlleleCode` for both haplotypes.
  - Route to MIS if both missing, HOM if `c0==c1`, else AMB.
- AMB lane semantics:
  - Build per-lane want_c0/want_c1 masks from `amb_code`; blend emissions per lane accordingly.
  - Do not split semantics by low/high 128-bit halves; half-splitting is allowed only for packing.
- Anchor gating: If `curr_abs_locus != ss.global_site_id`, treat as sibling and skip DP while maintaining bookkeeping (e.g., missing counters if applicable).
- Naming consistency: Match biallelic locals (`_prob`, `_emit`, `_nt`, `_tFreq`, `_mismatch`, `_factor`) and accumulators (`probSumH`, `probSumT`, `probSumK`).

## Performance & Memory
- Use `aligned_vector32<T>` for HMM arrays and supersite scratch buffers; avoid `std::vector` inside hot loops.
- Pre-allocate as class members: `ss_emissions`, `ss_emissions_h1`, and (optionally) cached donor codes `ss_cond_codes` with a parallel `ss_cached` bitset.
- Precompute emission vectors once per supersite event (INIT/RUN/COLLAPSE) rather than per function body; reuse the same AVX broadcast pattern as biallelic (`_factor`, `_tFreq`, `_nt`, `_mismatch`).

## Error Handling & Globals
- No exceptions in hot paths; use assertions for invariants. Underflow is signalled by return code and retried in double.
- Global singletons: exactly one `.cpp` defines `_DECLARE_TOOLBOX_HERE` to instantiate `rng`, `stb`, `alg`, `vrb`, `tac`; others include externs from `common/src/utils/otools.h`. In tests, `tests/src/test_toolbox.cpp` owns the declaration.
- Commit tracing: Makefile injects `__COMMIT_ID__` and `__COMMIT_DATE__` via `-D` flags; binaries report via `--version`.

## Debugging & Logging
- Debug builds: `make -C phase_common debug` adds `-g` and reduces optimization for easier stepping; run `gdb` or `valgrind` on binaries/tests.
- Logging: use `vrb.bullet()` and `vrb.title()` for structured console output. Inspect `Alpha`, `AlphaSum` (32‑byte aligned) in debugger for HMM state.

## Common Pitfalls (Do/Don’t)
- Don’t allocate inside `RUN_AMB` (or any hot path); hoist buffers to class members to avoid allocator churn and `posix_memalign` crashes.
- Don’t access supersite state without null/size checks (`super_sites`, `locus_to_super_idx`).
- Don’t perform unaligned AVX2 loads; ensure 32‑byte alignment for `_mm256_load_ps` sources.
- Do gate DP to anchors to prevent double counting at sibling splits.
- Do preserve AMB lane semantics from biallelic code; avoid half-lane logic divergence.

## Known Bugs (from current sprint)
**All Phase 3 blockers resolved as of Oct 29, 2025:**
- ✅ INIT_AMB overlapping store: Fixed by using single combined vector from 128-bit halves
- ✅ AMB lane mask semantics: Not applicable to Phase 3 (uses multinomial, not split biallelic)
- ✅ Per-split classification: Superseded by multinomial classification in backward pass
- ✅ Double counting at sibling loci: Fixed via anchor gating in backward loop
- ✅ HOM masquerading as AMB: Not applicable (Phase 3 classifies at supersite level)
- ✅ Repeated `ss_cond_codes` unpacking: Still occurs but not critical for Phase 3; can optimize later

**Remaining validation tasks:**
1. Integration testing on real data with --enable-supersites
2. Verify multinomials sum to 1.0 (add assertion if needed)
3. Validate mutual exclusivity in output BCF
4. Performance profiling (multinomial overhead vs. biallelic)

## Prioritized Work Queue
**Phase 3 Implementation Complete - Ready for Testing:**
1. ✅ Data structures (SuperSite.class_prob_offset, compute_job.SC, anchor_has_missing)
2. ✅ Backward pass multinomial computation (IMPUTE_SUPERSITE_MULTINOMIAL)
3. ✅ Genotype projection (sample from multinomial, project to splits)
4. ✅ Integration (phaser_algorithm passes SC, calls setSuperSiteContext)
5. ✅ Clean compilation (both single and double precision)

**Next Actions (Testing & Validation):**
1. Run integration tests: `bash test/scripts/phase.wgs.unrelated.sh --enable-supersites`
2. Add multinomial validation: assert Σ_c SC[...+c] ≈ 1.0 in IMPUTE_SUPERSITE_MULTINOMIAL
3. Add mutual exclusivity checker in output writer or post-processing script
4. Profile performance: compare runtime with/without supersites
5. Optional: SIMD optimization for IMPUTE_SUPERSITE_MULTINOMIAL (currently scalar stores)
6. Optional: Cache ss_cond_codes per supersite (10-20% speedup in SS-heavy windows)
7. Add smoke test fixtures under `tests/data/` with known multi-allelic VCFs

## Key Files
- `phase_common/src/models/haplotype_segment_single.{h,cpp}`: Float HMM core + supersite integration
- `phase_common/src/models/haplotype_segment_double.{h,cpp}`: Double HMM core + supersite integration
- `phase_common/src/models/super_site_macros.h`: Emission macro templates (INIT/RUN/COLLAPSE)
- `phase_common/src/objects/super_site_builder.{h,cpp}`: Multi‑allelic grouping + panel encoding
- `phase_common/src/phaser/phaser_algorithm.cpp`: Window segmentation, threading, precision fallback
- `common/src/utils/otools.h`: `aligned_vector32`, toolbox externs
- `tests/makefile`: Unit test build rules, external object linking

## Testing Guidance
- Always run `make clean && make -C tests -j8` (project policy) followed by `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib tests/bin/<name>` to execute unit tests.
- The new HMM state harness depends on the supersite builder; ensure `tests/src/test_toolbox.cpp` continues to define `_DECLARE_TOOLBOX_HERE` so globals are instantiated exactly once.

## Conventions & References
- Gate supersite logic by checking `super_sites && locus_to_super_idx` to preserve baseline behaviour.
- Restrict new buffers to `aligned_vector32<T>` unless stack or STL containers are intentionally used (document justification when deviating).
- Keep supersite-related fixtures synthetic for now; real data will live under `tests/data/` once smoke tests are added.
- External refs: HTSlib docs, Intel AVX2 intrinsics guide, Li & Stephens (2003).


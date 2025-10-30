## Project
SHAPEIT5 — `phase_common` super-site integration

## Source of Truth
- This AGENTS.md mirrors `.github/copilot-instructions.md`, which is the authoritative, detailed guide distilled from `SUPERSITE_CONVERSATION_SUMMARY.md`.
- If anything diverges, prefer `.github/copilot-instructions.md` and update this file to match.

## Objectives (unchanged)
- Expose optional super-site phasing for multi-allelic loci (`--enable-supersites`) while keeping default behaviour untouched.
- Maintain the original biallelic code paths (and performance characteristics) when supersites are disabled.

## Current Status (Oct 29)
- Data model / builder: `buildSuperSites` partitions split records, fills `super_sites`, `locus_to_super_idx`, `super_site_var_index`, and packs panel haplotype codes. Conditioning haplotypes are encoded with 4-bit classes (two per byte) as planned.
- HMM integration: `haplotype_segment_{single,double}` accept optional supersite state, derive sample codes on demand, and emit multiallelic probabilities for INIT / RUN / COLLAPSE / MIS. AVX2 double path matches expectations; float path is functionally correct but still relies on scalar shuffles for AMB lanes.
- Unit test coverage:
  - `test_supersite_emissions` (float/double kernels)
  - `test_supersite_accessor` (per-hap sample codes)
  - `test_supersite_unpack` (4-bit decode)
  - `test_supersite_builder` (grouping/panel packing)
  - `test_supersite_hmm` (forward path regression)
  - `test_supersite_hmm_states` (INIT/RUN/COLLAPSE micro-harness; currently forward-only after latest crash)
- Build / tooling: tests have a standalone `tests/makefile` mirroring top-level flags and linking supersite objects.

## Recent Work
- Investigated float RUN_AMB double-counting; refactored supersite branches to split lower/upper lanes explicitly instead of using `_mm256_blend_ps`.
- Added deterministic supersite HMM state harness comparing double-precision buffers against recorded expectations and ensuring the float path stays finite.
- Began wiring backward-path validation; currently disabled because `haplotype_segment_single::RUN_AMB` still allocates per-iteration temporaries (switch to `std::vector<float>` exposed a `posix_memalign` crash during backward).

## Outstanding Issues
1. Float RUN_AMB supersite lane fix
   - Forward supersite AMB handling now splits lanes correctly, but backward still uses heap allocations (`std::vector<float>` inside the inner loop). Pre-allocate scratch buffers on the object to avoid heap churn and re-enable backward testing.
2. Unit test harness
   - `tests/bin/test_supersite_hmm_states` fails when backward is run; re-enable assertions once allocator fix lands.
   - Fixtures under `tests/data/` still needed for smoke tests (placeholders exist).
3. Build stability
   - `make -C tests` succeeds after setting `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib`. Ensure CI/exported instructions include this.

## Next Steps
1. Hoist supersite emission scratch buffers (`ss_emissions`, `ss_emissions_h1`) into persistent `aligned_vector32` members so both forward and backward avoid heap churn and are ready for SIMD again.
2. Re-enable and harden backward-path comparisons in `test_supersite_hmm_states` once allocator churn is resolved.
3. Author the integration smoke test harness (`tests/run_supersite_smoke.sh`) once VCF fixtures are assembled.
4. Profile float supersite paths to decide whether an AVX2 implementation is warranted or if scalar remains adequate.

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
- Runtime lookups:
  - Panel donor code: `unpackSuperSiteCode(panel_codes, ss.panel_offset, hap_idx)`
  - Sample code: `getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap)`
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
1. INIT_AMB overlapping store: two 256‑bit stores corrupt next block. Fix by composing a single combined vector from 128‑bit halves before store.
2. AMB lane mask semantics broken when splitting halves. Fix with per‑lane want_c0/want_c1 masks from `amb_code`.
3. Per‑split classification used in supersite paths. Fix by centralized supersite classification of both haplotypes before dispatch.
4. Double counting at sibling loci. Fix with anchor gating (`curr_abs_locus == ss.global_site_id`).
5. HOM masquerading as AMB when `c0==c1`. Fix by routing to HOM in classification.
6. Repeated `ss_cond_codes` unpacking in hot paths. Fix via per‑supersite caching and invalidation on PBWT changes.

## Prioritized Work Queue
1. Fix INIT_AMB overlapping store (correctness).
2. Fix AMB lane mask semantics (correctness).
3. Implement anchor gating (correctness/perf).
4. Centralize supersite classification (correctness/maintainability).
5. Extract `SS_*` functions to mirror biallelic taxonomy (maintainability/reviewability).
6. Add `ss_cond_codes` caching (performance in SS‑heavy windows).
7. Re-enable backward pass assertions in `test_supersite_hmm_states` after allocator fixes.
8. Add smoke test fixtures under `tests/data/` and `tests/run_supersite_smoke.sh`.

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


## Project
SHAPEIT5 — `phase_common` super-site integration

## Objectives (unchanged)
- Expose optional super-site phasing for multi-allelic loci (`--enable-supersites`) while keeping default behaviour untouched.
- Maintain the original biallelic code paths (and performance characteristics) when supersites are disabled.

## Current Status (Oct 29)
- **Data model / builder**: `buildSuperSites` partitions split records, fills `super_sites`, `locus_to_super_idx`, `super_site_var_index`, and packs panel haplotype codes. Conditioning haplotypes are encoded with 4-bit classes (two per byte) as planned.
- **HMM integration**: `haplotype_segment_{single,double}` accept optional supersite state, derive sample codes on demand, and emit multiallelic probabilities for INIT / RUN / COLLAPSE / MIS. AVX2 double path matches expectations; float path is functionally correct but still relies on scalar shuffles for AMB lanes.
- **Unit test coverage**:  
  - `test_supersite_emissions` (float/double kernels)  
  - `test_supersite_accessor` (per-hap sample codes)  
  - `test_supersite_unpack` (4-bit decode)  
  - `test_supersite_builder` (grouping/panel packing)  
  - `test_supersite_hmm` (forward path regression)  
  - `test_supersite_hmm_states` (new micro-harness anchoring INIT/RUN/COLLAPSE for supersites; currently exercises `forward()` only after the latest crash)
- **Build / tooling**: tests now own a standalone makefile (`tests/makefile`) that mirrors the top-level flags and links against the supersite objects.

## Recent Work
- Investigated float RUN_AMB double-counting; refactored supersite branches to split lower/upper lanes explicitly instead of relying on `_mm256_blend_ps`.
- Added deterministic supersite HMM state harness comparing double-precision buffers against recorded expectations and ensuring the float path stays finite.
- Began wiring backward-path validation, but it is temporarily disabled because `haplotype_segment_single::RUN_AMB` still allocates per-iteration temporaries on the hot path (switching to `std::vector<float>` exposed a `posix_memalign` crash during backward).

## Outstanding Issues
1. **Float RUN_AMB supersite lane fix**  
   - Forward supersite AMB handling now splits lanes correctly, but backward still uses heap allocations (`std::vector<float>` inside the inner loop) as a stop-gap. Need to pre-allocate scratch buffers on the object to avoid repeated allocations and re-enable backward testing.
2. **Unit test harness**  
   - `tests/bin/test_supersite_hmm_states` currently fails when backward is run; re-enable assertions once the allocator fix lands.  
   - We still need fixtures under `tests/data/` for future smoke tests (placeholder directories exist).
3. **Build stability**  
   - `make -C tests` succeeds after running with `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib`. Ensure CI/exported instructions include this requirement.

## Next Steps
1. Hoist supersite emission scratch buffers (`ss_emissions`, `ss_emissions_h1`) into persistent `aligned_vector32` members so both forward and backward avoid heap churn and ARE ready for SIMD again.
2. Re-enable and harden backward-path comparisons in `test_supersite_hmm_states` once allocator churn is resolved.
3. Author the integration smoke test harness (`tests/run_supersite_smoke.sh`) once VCF fixtures are assembled.
4. Profile float supersite paths to decide whether an AVX2 implementation is warranted or if scalar remains adequate.

## Testing Guidance
- Always run `make clean && make -C tests -j8` (project policy) followed by `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib tests/bin/<name>` to execute unit tests.
- The new HMM state harness depends on the supersite builder; ensure `tests/src/test_toolbox.cpp` continues to define `_DECLARE_TOOLBOX_HERE` so globals are instantiated exactly once.

## Notes / Conventions
- Continue to gate supersite logic by checking `super_sites && locus_to_super_idx` to preserve baseline behaviour.
- Restrict new buffers to `aligned_vector32<T>` unless stack or STL containers are intentionally used (document justification when deviating).
- Keep supersite-related fixtures synthetic for now; real data will live under `tests/data/` once smoke tests are added.

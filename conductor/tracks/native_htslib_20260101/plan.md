# Track Plan: Refactor `xcf_reader` to Native HTSlib

## Phase 1: Preparation and C++ Style Guide
- [x] Task: Create a specific C++ style guide based on the codebase conventions
- [x] Task: Analyze `xcf_reader` usage and internal implementation to list all features being removed (binary/sparse support).

## Phase 2: Refactor Scanning Loop
- [~] Task: Create a reproduction test case (smoke test) that runs the current `phase_common` on a small dataset.
- [ ] Task: Refactor `genotype_reader_scaning.cpp` to use `bcf_srs_t` directly.
    - Remove `xcf_reader` initialization.
    - Implement `bcf_sr_init` and `bcf_sr_add_reader`.
    - Rewrite the `while` loop using `bcf_sr_next_line`.
    - Replace metadata accessors (chr, pos, ref, alt) with direct `bcf1_t` struct access.
- [ ] Task: Verify compilation and run the smoke test to ensure scanning (variant counting) still works.

## Phase 3: Refactor Reading Loop
- [ ] Task: Refactor `genotype_reader_reading.cpp` to use `bcf_srs_t` directly.
    - Re-implement reader initialization (similar to scanning).
    - Replace `XR.readRecord` with `bcf_get_genotypes`.
    - **Critical:** Implement manual memory management for the genotype buffer ( `dst` and `ndst` for `bcf_get_genotypes`).
    - Remove all `if (type == RECORD_BINARY_...)` branches.
- [ ] Task: Verify compilation and run the smoke test to ensure full phasing execution works.

## Phase 4: Cleanup and Verification
- [ ] Task: Remove `xcf_reader` includes from `genotype_reader_header.h` and other modified files.
- [ ] Task: Run a regression check on a biallelic dataset to ensure output identity (basic parity check).
- [ ] Task: Conductor - User Manual Verification 'Cleanup and Verification' (Protocol in workflow.md)

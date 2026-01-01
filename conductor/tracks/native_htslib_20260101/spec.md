# Track Specification: Refactor `xcf_reader` to Native HTSlib

## 1. Objective
To modernize and simplify the input/output layer of `phase_common` by removing the custom `xcf_reader` wrapper and replacing it with direct calls to `htslib`'s synced reader (`bcf_srs_t`). This refactoring removes the legacy support for the custom "binary/sparse" (`.bin`) format, paving the way for native multiallelic support by using standard VCF/BCF handling.

## 2. Scope
- **Target Project:** `phase_common` (initially).
- **Files to Modify:**
    - `phase_common/src/io/genotype_reader/genotype_reader_scaning.cpp`
    - `phase_common/src/io/genotype_reader/genotype_reader_reading.cpp`
    - `phase_common/src/io/genotype_reader/genotype_reader_header.h` (to remove include)
- **Files to Delete/Deprecate:**
    - Usage of `xcf_reader` class.
    - Logic related to `RECORD_BINARY_GENOTYPE`, `RECORD_BINARY_HAPLOTYPE`, `RECORD_SPARSE_HAPLOTYPE`.

## 3. Requirements

### 3.1 Functional Requirements
- **Direct HTSlib Usage:** The scanning and reading loops must use `bcf_sr_init`, `bcf_sr_add_reader`, `bcf_sr_next_line`, and `bcf_get_genotypes` directly.
- **Multiallelic Awareness:** The new implementation must correctly identify multiallelic sites using `bcf1_t->n_allele`.
- **Memory Management:** The refactored code must properly manage memory for genotype buffers (allocating and freeing `int32_t *` buffers used by `bcf_get_genotypes`), replicating the safety that `xcf_reader` might have provided.
- **Feature Parity (VCF/BCF):** Reading of standard VCF and BCF files must remain functional and identical in behavior to the previous version.
- **Feature Removal:** Explicitly remove support for reading the custom `.bin` + `.vcf` XCF format.

### 3.2 Non-Functional Requirements
- **Performance:** The overhead of direct `htslib` calls should be negligible or better than the wrapper.
- **Code Style:** Strict adherence to the existing project style (naming conventions, indentation).

## 4. Verification Plan
- **Compilation:** The project must compile without errors using `make`.
- **Smoke Test:** Run the `phase_common` binary on a small test dataset (e.g., from `test/` directory) and verify it runs to completion without crashing.
- **Output Validation:** Ensure the output phased file is valid BCF/VCF.

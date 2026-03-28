# Analysis of xcf_reader Removal

## Features being Removed
The following features provided by the `xcf_reader` wrapper will be removed in favor of standard HTSlib VCF/BCF support:

1. **XCF Binary Format Support:**
   - Reading `.bin` files associated with VCFs.
   - Parsing `INFO/SEEK` fields to jump to binary data offsets.
   - `RECORD_BINARY_GENOTYPE`: 2-bit packed genotype reading.
   - `RECORD_BINARY_HAPLOTYPE`: 1-bit packed haplotype reading.
   - `RECORD_SPARSE_HAPLOTYPE`: Sparse haplotype reading.

2. **Custom Pedigree Reading (for Binary files):**
   - Reading `.fam` files for sample metadata when using binary files (HTSlib reads samples from VCF header).

## Refactoring Strategy

### Scanning Phase (`genotype_reader_scaning.cpp`)
- **Current:** Uses `XR.nextRecord()` to iterate and `XR.typeRecord()` to validate formats.
- **New:** Use `bcf_sr_next_line(sr)` to iterate.
- **Action:** Remove all checks for `RECORD_BINARY_*`. Assume standard BCF/VCF input.

### Reading Phase (`genotype_reader_reading.cpp`)
- **Current:** Uses `XR.readRecord()` which dispatches to either `bcf_get_genotypes` or `bin_fds.read`.
- **New:** Call `bcf_get_genotypes` directly.
- **Action:** Delete the `else if (main_type == RECORD_BINARY_GENOTYPE)` blocks.
- **Action:** Delete the `else if (ref_type == RECORD_BINARY_HAPLOTYPE)` and `SPARSE` blocks.

## Memory Management
- `xcf_reader` managed its own buffers.
- **New:** We must allocate `int32_t *` buffers (`main_buffer`, `ref_buffer`, `scaf_buffer`) and pass them to `bcf_get_genotypes`, handling re-allocation (though `bcf_get_genotypes` handles realloc if size passes, we just need to track the size).

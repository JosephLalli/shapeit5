# SHAPEIT5 Unit Test Suite

## Overview

This directory contains comprehensive unit tests for SHAPEIT5, with a particular focus on validating the Phase 3 multinomial imputation implementation for missing multiallelic sites.

## Test Organization

### Core Infrastructure

- **`test_toolbox.cpp`** - Global singleton instantiation (`rng`, `vrb`, `tac`, etc.)
  - Must be compiled with `-D_DECLARE_TOOLBOX_HERE` to instantiate globals exactly once
  - All other test files extern-declare these globals

### Supersite Multiallelic Support Tests

#### 1. test_supersite_emissions.cpp
**Purpose**: Validate emission probability computation for supersite HMM

**Coverage:**
- Match emission = 1.0 (donor haplotype matches sample)
- Mismatch emission = ed/ee (error rate ratio)
- Missing emission = 1.0 (uninformative, imputed via flanking sites)
- Double precision emission consistency

**Key Assertions:**
- Emission values match expected HMM parameters
- Float and double precision produce consistent results

---

#### 2. test_supersite_unpack.cpp
**Purpose**: Validate 4-bit allele code unpacking from packed panel buffer

**Coverage:**
- Single code unpacking: `unpackSuperSiteCode(buffer, offset, hap_idx)`
- Batch unpacking: `unpackSuperSiteCodesBatch()` for multiple haplotypes
- Vectorized unpacking: `unpackSuperSiteCodesVectorized_PEXT()` using PEXT instruction
- Offset handling: Verify correct byte offset calculation for multi-supersite buffers
- REF (code 0) vs ALT1-15 (codes 1-15) differentiation

**Key Assertions:**
- Packing format: 2 codes per byte (lower 4 bits = even hap, upper 4 bits = odd hap)
- All unpacking methods produce identical results
- Offset arithmetic correctly handles multi-supersite panels

---

#### 3. test_supersite_accessor.cpp
**Purpose**: Validate SuperSite data structures and accessor functions

**Coverage:**
- SuperSite structure creation and field initialization
- Panel haplotype code unpacking via `unpackSuperSiteCode()`
- Sample genotype code inference via `getSampleSuperSiteAlleleCode()`
  - Returns 0 for REF
  - Returns 1-n_alts for ALT1..ALTn (first matching split)
  - Returns `SUPERSITE_CODE_MISSING` for missing data
- Multi-allelic supersites (2-15 ALTs)
- aligned_vector32 32-byte alignment verification

**Key Assertions:**
- SuperSite fields correctly populated by builder
- Code inference matches split record genotypes
- Missing data properly detected across all splits
- Memory alignment meets AVX2 requirements

---

#### 4. test_supersite_builder.cpp
**Purpose**: Validate multiallelic variant grouping and panel encoding

**Coverage:**
- Basic 2-split supersite grouping by chr:bp
- Locus-to-supersite index mapping (`locus_to_super_idx[v] → ss_idx`)
- Variant index array construction (CSR-style flat array)
- Panel code packing:
  - 4-bit codes per conditioning haplotype
  - 0 = REF, 1-15 = which ALT carried
  - 2 codes per byte (packed)
- Multiple supersites at different genomic positions
- No-supersites case (all biallelic variants)

**Key Assertions:**
- All split records at same position grouped into one SuperSite
- Panel codes correctly inferred from split record ALT patterns
- `locus_to_super_idx[v] = -1` for non-supersite variants
- Anchor variant is first split record in supersite

---

#### 5. test_missing_multiallelic_multinomial.cpp
**Purpose**: Validate Phase 3 multinomial imputation for missing supersites

**Coverage:**
- 3-split multiallelic site creation (4 classes: REF, ALT1, ALT2, ALT3)
- Conditioning panel with diverse allele distribution
- Missing data detection across all splits of supersite
- Supersite context setting for genotype sampling
- Structure validation for multinomial posterior computation

**Key Assertions:**
- All splits at same position grouped correctly
- Conditioning haplotypes encoded with correct class codes (0=REF, 1-3=ALT1-3)
- `getSampleSuperSiteAlleleCode()` returns `SUPERSITE_CODE_MISSING` when all splits missing
- `setSuperSiteContext()` properly initializes genotype for multinomial sampling

**Note**: End-to-end multinomial sampling (calling `genotype::make()`) requires full HMM forward/backward pass to populate SC buffer. This is tested in integration tests.

---

### Legacy Tests (Expected to Fail)

- **`test_supersite_hmm.cpp`** - Tests old HMM implementation (pre-Phase 3)
- **`test_supersite_hmm_states.cpp`** - Tests old state transitions (pre-Phase 3)

These tests validate pre-Phase 3 biallelic supersite logic. They are retained for reference but skip execution in `make test-run` as the implementation has evolved to Phase 3 multinomial imputation.

---

## Building and Running Tests

### Prerequisites

- C++17 compiler with AVX2 support
- Boost libraries: `program_options`, `iostreams`, `serialization`
- HTSlib 1.12+
- Git submodules initialized (`git submodule update --init --recursive`)

### Build All Tests

```bash
make -C tests clean
make -C tests -j$(nproc)
```

### Run All Tests

```bash
# Run all passing tests (skips legacy tests)
make -C tests test-run

# Or manually with proper library path
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_supersite_emissions
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_supersite_accessor
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_supersite_unpack
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_supersite_builder
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_missing_multiallelic_multinomial
```

### Individual Test Execution

```bash
# Example: Run emission tests
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_supersite_emissions

# Example: Run multinomial imputation validation
LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH tests/bin/test_missing_multiallelic_multinomial
```

---

## Test Design Principles

### 1. Minimal Dependencies
- Tests use synthetic data (no external VCF fixtures required for unit tests)
- Variants created programmatically via `variant` constructor
- Conditioning panels built in-memory via `conditioning_set` API

### 2. Toolbox Singleton Pattern
- **Exactly one** `.cpp` defines `_DECLARE_TOOLBOX_HERE` (`test_toolbox.cpp`)
- All other tests `#include "../../common/src/utils/otools.h"` to extern-declare globals
- Avoids multiple definition errors in linker

### 3. Private Member Access
- Tests use `#define private public` / `#define protected public` to access HMM internals
- Required for white-box validation of Alpha/Beta arrays, state transitions, etc.
- Undefined immediately after includes to avoid leakage

### 4. Alignment Verification
- Supersite tests verify 32-byte alignment for `aligned_vector32<T>` members
- Critical for AVX2 `_mm256_load_ps()` instructions (segfault on misalignment)
- Example: `assert((reinterpret_cast<uintptr_t>(vec.data()) % 32) == 0)`

---

## Phase 3 Multinomial Imputation Summary

### Problem
Independent biallelic imputation of split multiallelic records can produce multiple ALTs on the same haplotype (e.g., both split1=1|0 and split2=1|0 at chr:pos, implying hap0 carries ALT1 and ALT2 simultaneously).

### Solution
Native multinomial imputation:

1. **HMM Backward Pass** (`IMPUTE_SUPERSITE_MULTINOMIAL`):
   - Computes P(class_c | Alpha, Beta) for all classes c ∈ {REF, ALT1, ..., ALTn}
   - Uses C SIMD accumulators (C ≤ 16): `sum[code] += Alpha[k] × Beta[k] / AlphaSum`
   - Stores normalized posteriors in SC buffer: `SC[offset + hap*C + c]`

2. **Genotype Sampling** (`genotype::make`):
   - Samples one class per haplotype from multinomial distribution
   - Projects to splits: if class=ALTi, set split_i=ALT, all others=REF
   - **Mathematical guarantee**: Exactly one split set to ALT per haplotype

### Validation Strategy

**Unit Tests** (this directory):
- Data structure correctness (SuperSite, allele codes, mappings)
- Code unpacking accuracy (4-bit packed format)
- Missing data detection across splits
- Supersite context setting

**Integration Tests** (`test/scripts/`):
- End-to-end phasing with `--enable-supersites`
- BCF output validation (no multiple ALTs per haplotype)
- Phasing accuracy vs. biallelic mode
- Performance profiling

---

## Test Coverage Metrics

| Component | Unit Tests | Integration Tests |
|-----------|-----------|-------------------|
| Supersite builder | ✅ 6 scenarios | ✅ phase.*.sh |
| 4-bit code packing/unpacking | ✅ 5 scenarios | ✅ phase.*.sh |
| Accessor functions | ✅ 6 scenarios | ✅ phase.*.sh |
| Missing data detection | ✅ 2 scenarios | ✅ phase.wgs.*.sh |
| Multinomial structure | ✅ 5 scenarios | ⚠️ Pending |
| Mutual exclusivity | ⚠️ Requires HMM | ⚠️ Pending |

**Legend:**
- ✅ Implemented and passing
- ⚠️ Pending implementation
- ❌ Known failure

---

## Future Enhancements

### High Priority
1. **Full HMM integration test** - End-to-end forward/backward with multinomial imputation
2. **Mutual exclusivity validator** - Post-processing script to check BCF outputs
3. **Smoke test fixtures** - Real VCF examples under `tests/data/` with known ground truth

### Medium Priority
4. **Performance benchmarks** - Multinomial overhead vs. biallelic (target: <10% slowdown)
5. **Numerical stability tests** - Edge cases (zero probabilities, underflow, overflow)
6. **SIMD optimization tests** - Validate AVX2 vectorization of multinomial computation

### Low Priority
7. **Fuzzing harness** - Random multiallelic site generation for stress testing
8. **Comparison vs. Phase 1** - Validate Phase 3 produces same/better results as reconstruction approach

---

## Troubleshooting

### Common Build Issues

**Problem**: `boost/program_options.hpp: No such file or directory`
```bash
# Solution: Install boost dev packages
sudo apt install -y libboost-program-options-dev libboost-iostreams-dev libboost-serialization-dev
```

**Problem**: `htslib/hts.h: No such file or directory`
```bash
# Solution: Build and install htslib
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xf htslib-1.16.tar.bz2
cd htslib-1.16
make -j$(nproc)
sudo make install
```

**Problem**: `cannot find -lhts` at link time
```bash
# Solution: Update library cache
sudo ldconfig
```

### Common Runtime Issues

**Problem**: `error while loading shared libraries: libhts.so.3`
```bash
# Solution: Add /usr/local/lib to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
# Or use the test-run target which sets this automatically
make -C tests test-run
```

**Problem**: Segmentation fault in test
```bash
# Likely causes:
# 1. Unaligned AVX2 load (check aligned_vector32 usage)
# 2. Null pointer dereference (check supersite context setup)
# 3. Buffer overflow (check SC buffer sizing in multinomial tests)

# Debug with:
gdb tests/bin/test_name
(gdb) run
(gdb) bt  # backtrace on crash
```

---

## References

- **AGENTS.md** - Detailed implementation notes for Phase 3 multinomial imputation
- **SUPERSITE_CONVERSATION_SUMMARY.md** - Evolution of supersite design decisions
- **.github/copilot-instructions.md** - Authoritative coding guidelines for supersites
- **test/scripts/** - Integration test suite for end-to-end validation

---

## Contact

For questions about test design or to report test failures, please open an issue on the SHAPEIT5 GitHub repository.

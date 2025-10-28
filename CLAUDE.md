# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with SHAPEIT5 source code.

## Project Overview

SHAPEIT5 is a state-of-the-art haplotype phasing software designed for large-scale genomic datasets. It implements novel PBWT-based algorithms for accurately phasing both common and rare variants across hundreds of thousands of samples.

**Key features**:
- Two-stage phasing: common variants (MAF ≥0.1%) followed by rare variant integration
- PBWT data structures for efficient identity-by-descent detection
- Scalable to 150,000+ whole genome sequences (validated on UK Biobank)
- Switch error rates <5% for variants present in 1/100,000 samples
- Supports pedigree information, haploid regions, and reference scaffolding

## Architecture Overview

### Core Components
- **`phase_common/`**: Phases common variants (typically SNP array data)
- **`phase_rare/`**: Phases rare variants onto common variant scaffold  
- **`switch/`**: Computes switch error rates and genotyping accuracy
- **`ligate/`**: Combines multiple phased chunks into chromosome-length files
- **`simulate/`**: Generates simulated haplotype datasets for testing

### Key Source Directories
- **`phase_common/src/`**: Main phasing algorithms and PBWT implementation
- **`test/`**: Comprehensive test suite with optimized datasets
- **`static_bins/`**: Pre-compiled static executables for deployment
- **`resources/`**: Genetic maps and reference data

## Development Guidelines

### Building SHAPEIT5
```bash
make all              # Build all components (dynamic linking)
make static_exe       # Build static executables 
make clean           # Clean build artifacts

# Build specific component
cd phase_common && make
cd phase_common && make static_exe
```

### Testing System
SHAPEIT5 uses a **headerless VCF.gz + MD5 validation system** for deterministic testing:

```bash
# Run all tests
cd test && for s in scripts/phase*.sh; do bash "$s"; done

# Run specific test scenarios
cd test/scripts
./phase.array.reference.sh    # Array data with reference panel
./phase.wgs.family.sh         # WGS data with pedigree information
./phase.array.haploid.sh      # Haploid region phasing
```

**Important**: All tests use:
- **Single-threaded execution** (`--thread 1`) for reproducibility
- **Optimized datasets**: 5k samples, 181-10k variants for fast execution
- **MD5 checksum validation** to avoid header formatting issues

### Test Validation Functions
Located in `test/scripts/lib/test_utils.sh`:
- `extract_variants()`: Converts BCF to headerless VCF.gz
- `assert_same_md5()`: Validates output against expected MD5 checksums
- `generate_headerless_vcf_gz()`: Creates ground truth files

### Common Development Tasks

#### Adding New Functionality
1. **Modify source code** in relevant component (e.g., `phase_common/src/`)
2. **Update command-line options** in main.cpp if needed
3. **Add tests** following existing patterns in `test/scripts/`
4. **Regenerate expected outputs** using published SHAPEIT5 binaries
5. **Verify all tests pass** before committing

#### Debugging and Profiling
```bash
# Debug build with symbols
cd phase_common && make debug

# No optimization for detailed debugging  
cd phase_common && make CXXFLAG="-g -O0"

# Memory profiling available through debug builds
```

#### Working with Test Data
- **Never modify test input files** (array/*.bcf, wgs/*.bcf) directly
- **Use published SHAPEIT5 binaries** to generate expected outputs
- **Preserve deterministic testing** by using MD5 validation
- **Keep test execution fast** (target: 10-20 seconds per test)

### Code Style and Conventions
- **C++17 codebase** with AVX2/FMA optimizations
- **Dependencies**: HTSlib for VCF/BCF I/O, Boost for utilities
- **Naming**: Use existing patterns in codebase
- **Comments**: Focus on algorithmic logic and complex calculations
- **Error handling**: Use proper exception handling and logging

### Performance Considerations
- **Release builds use -O3 optimization** for production performance
- **PBWT algorithms are memory-intensive**: consider chunking for large datasets
- **Multi-threading available** but use single-thread for testing reproducibility
- **Static linking preferred** for deployment to avoid dependency issues

### Algorithmic Background
SHAPEIT5 implements advanced population genetics algorithms:
- **PBWT (Positional Burrows-Wheeler Transform)**: Efficient haplotype matching
- **HMM (Hidden Markov Model)**: Statistical phasing with recombination maps
- **IBD detection**: Identity-by-descent for rare variant phasing
- **Sparse representation**: Memory-efficient storage for rare variants

### Documentation and Resources
- **Main documentation**: `docs/` directory with Jekyll-based site
- **Example workflows**: `example/` directory with sample commands
- **Version tracking**: `versions/` directory with changelog
- **External tools**: `external/` directory with dependencies

### Troubleshooting Common Issues
1. **Compilation errors**: Check HTSlib and Boost installations
2. **Test failures**: Verify single-threaded execution and input file integrity
3. **Memory issues**: Consider chunking large chromosomes or increasing available RAM
4. **Performance problems**: Profile with debug builds and check PBWT parameters

### Contributing Guidelines
- **Follow existing code patterns** and architectural decisions
- **Maintain backward compatibility** for command-line interfaces
- **Add comprehensive tests** for new functionality
- **Update documentation** for user-visible changes
- **Use meaningful commit messages** describing changes and rationale

This framework provides the foundation for reliable, scalable haplotype phasing in population genomics applications.
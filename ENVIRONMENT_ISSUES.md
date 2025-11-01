# Environment Configuration Issues & Solutions

## Problem Summary

The SHAPEIT5 test suite was failing to build in the Copilot environment due to library path mismatches and missing dependencies.

## Root Causes

### 1. Library Path Mismatch

**Issue**: The main makefile (`common/makefile_common.mk`) expects libraries in Homebrew paths:
- HTSlib: `~/.linuxbrew/include` and `~/.linuxbrew/lib`
- Boost: `~/.linuxbrew/include` and `~/.linuxbrew/lib`

**Reality**: The Copilot environment has libraries in system paths:
- HTSlib: `/usr/local/include` and `/usr/local/lib`
- Boost: `/usr/include` and `/usr/lib/x86_64-linux-gnu`

### 2. Missing Dependencies

**Issue**: The following were not installed:
- Boost development libraries (`libboost-program-options-dev`, `libboost-iostreams-dev`, `libboost-serialization-dev`)
- HTSlib 1.16 (required version for BCF/VCF support)
- xcftools submodule (not initialized)

### 3. Runtime Linking

**Issue**: Even after successful compilation, tests failed at runtime because `LD_LIBRARY_PATH` didn't include `/usr/local/lib` where HTSlib was installed.

## Comparison: CI Workflow vs Makefile Expectations

| Component | CI Workflow (.github/workflows/build.yml) | Makefile Expectation | Copilot Reality |
|-----------|-------------------------------------------|----------------------|-----------------|
| HTSlib location | `/usr/local` (built from source) | `~/.linuxbrew` | `/usr/local` |
| Boost location | `/usr/include` (apt install) | `~/.linuxbrew` | `/usr/include` |
| Runtime LD_LIBRARY_PATH | Not explicitly set | `~/.linuxbrew/lib` | Needs `/usr/local/lib` |

## Solution Options Considered

### Option 1: Modify Main Makefile (❌ Rejected)
**Pros**: One build system for all environments  
**Cons**: Would break existing development workflows, affects production code

### Option 2: Use Environment Variables (❌ Rejected)
**Pros**: No code changes  
**Cons**: Brittle, easy to forget, not documented in repository

### Option 3: Create Copilot-Specific Makefiles (✅ Implemented)
**Pros**:
- Isolates Copilot configuration from production
- Clearly documented in repository
- Easy for future Copilot agents to find and use

**Cons**:
- Slight duplication of makefile logic

## Implemented Solution

### Created Files:

1. **`common/makefile_copilot.mk`** - Overrides library paths for Copilot environment
   ```makefile
   HTSLIB_INC = /usr/local/include
   HTSLIB_LIB = /usr/local/lib
   BOOST_INC = /usr/include
   BOOST_LIB_IO = /usr/lib/x86_64-linux-gnu/libboost_iostreams.so
   BOOST_LIB_PO = /usr/lib/x86_64-linux-gnu/libboost_program_options.so
   ```

2. **`tests/makefile.copilot`** - Test build configuration for Copilot
   ```makefile
   include ../common/makefile_copilot.mk
   include makefile  # Then include regular test makefile
   ```

3. **`.github/workflows/copilot.yml`** - CI workflow for Copilot agents
   - Installs Boost via apt
   - Builds HTSlib 1.16 from source
   - Initializes xcftools submodule
   - Sets LD_LIBRARY_PATH correctly

### Dependency Installation Steps:

```bash
# Install Boost
sudo apt-get update
sudo apt-get install -y \
    libboost-program-options-dev \
    libboost-iostreams-dev \
    libboost-serialization-dev

# Build HTSlib 1.16
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/usr/local
make -j$(nproc)
sudo make install

# Initialize xcftools submodule
git submodule update --init --recursive
```

## Build Instructions for Copilot Agents

### Standard Build (Production)
```bash
make -C tests clean
make -C tests -j$(nproc)
make -C tests test-run
```

### Copilot Build
```bash
cd tests
make -f makefile.copilot clean
make -f makefile.copilot -j$(nproc)
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
make -f makefile.copilot test-run
```

## Impact Assessment

**Who is affected**:
- ✅ Copilot agents (now have working build system)
- ❌ Regular developers (no change - still use standard makefile)
- ❌ CI/CD (no change - GitHub Actions unchanged)

**Risk**: Low
- Copilot-specific files are additive only
- No modifications to production build system
- Clearly documented and isolated

## Verification

All tests now compile and run successfully in Copilot environment:
- ✅ test_supersite_emissions: PASS
- ✅ test_supersite_accessor: PASS
- ✅ test_supersite_unpack: PASS
- ✅ test_supersite_builder: PASS
- ✅ test_missing_multiallelic_multinomial: PASS

## Future Recommendations

1. **Standardize CI and development environments** - Consider containerization (Docker) to ensure consistency
2. **Document environment requirements** - Add README section explaining library requirements
3. **Automated environment detection** - Makefile could auto-detect library locations
4. **Unified build system** - Consider CMake for cross-platform builds

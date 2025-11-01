# Copilot Agent Instructions for SHAPEIT5 Testing

## Quick Start

If you're a Copilot agent working on SHAPEIT5, follow these instructions to set up and run tests.

## Environment Differences

The standard SHAPEIT5 build expects libraries in Homebrew paths (`~/.linuxbrew`), but the Copilot environment has them in system paths. We've created Copilot-specific makefiles to handle this.

| Library | Standard Location | Copilot Location |
|---------|-------------------|------------------|
| HTSlib headers | `~/.linuxbrew/include` | `/usr/local/include` |
| HTSlib library | `~/.linuxbrew/lib` | `/usr/local/lib` |
| Boost headers | `~/.linuxbrew/include` | `/usr/include` |
| Boost libraries | `~/.linuxbrew/lib` | `/usr/lib/x86_64-linux-gnu` |

## Initial Setup (First Time Only)

### 1. Install Dependencies

```bash
# Install Boost development libraries
sudo apt-get update
sudo apt-get install -y \
    libboost-program-options-dev \
    libboost-iostreams-dev \
    libboost-serialization-dev \
    build-essential \
    wget \
    bzip2

# Build and install HTSlib 1.16
cd /tmp
wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
tar -xjf htslib-1.16.tar.bz2
cd htslib-1.16
./configure --prefix=/usr/local
make -j$(nproc)
sudo make install
cd /home/runner/work/shapeit5/shapeit5

# Initialize xcftools submodule
git submodule update --init --recursive
```

### 2. Verify Installation

```bash
# Check HTSlib
ls /usr/local/include/htslib/  # Should show header files
ls /usr/local/lib/libhts.*     # Should show library files

# Check Boost
dpkg -l | grep libboost        # Should show installed Boost packages
```

## Building Tests

### Using Copilot Makefile (Recommended)

```bash
cd /home/runner/work/shapeit5/shapeit5/tests

# Clean previous builds
make -f makefile.copilot clean

# Build all tests
make -f makefile.copilot -j$(nproc)
```

### Build Output

Successful build creates binaries in `tests/bin/`:
- `test_supersite_emissions_real`
- `test_supersite_accessor`
- `test_supersite_unpack`
- `test_supersite_builder`
- `test_missing_multiallelic_multinomial`
- `test_supersite_vs_biallelic_simple`

## Running Tests

### Run All Tests

```bash
cd /home/runner/work/shapeit5/shapeit5/tests
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
make -f makefile.copilot test-run
```

### Run Individual Tests

```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

# Run specific test
./bin/test_supersite_emissions_real
./bin/test_supersite_accessor
./bin/test_supersite_unpack
./bin/test_supersite_builder
./bin/test_missing_multiallelic_multinomial
./bin/test_supersite_vs_biallelic_simple
```

## Troubleshooting

### Issue: "error while loading shared libraries: libhts.so.3"

**Cause**: LD_LIBRARY_PATH not set  
**Solution**:
```bash
export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
```

### Issue: "fatal error: htslib/vcf.h: No such file or directory"

**Cause**: HTSlib not installed  
**Solution**: Follow "Install Dependencies" section above

### Issue: "cannot find -lboost_program_options"

**Cause**: Boost libraries not installed  
**Solution**:
```bash
sudo apt-get install libboost-program-options-dev libboost-iostreams-dev libboost-serialization-dev
```

### Issue: Tests fail with segmentation fault

**Causes**:
1. LD_LIBRARY_PATH not set (see first issue)
2. Boost version mismatch

**Debug**:
```bash
# Check which libraries are loaded
ldd bin/test_supersite_accessor

# Should show:
#   libhts.so.3 => /usr/local/lib/libhts.so.3
#   libboost_iostreams.so.X => /usr/lib/x86_64-linux-gnu/...
```

### Issue: xcftools submodule empty

**Cause**: Submodule not initialized  
**Solution**:
```bash
git submodule update --init --recursive
```

## File Organization

```
shapeit5/
├── common/
│   ├── makefile_common.mk       # Standard makefile (uses Homebrew paths)
│   └── makefile_copilot.mk      # Copilot override (uses system paths)
├── tests/
│   ├── makefile                 # Standard test makefile
│   ├── makefile.copilot         # Copilot test makefile (USE THIS)
│   ├── src/                     # Test source files
│   └── bin/                     # Compiled test binaries
├── .github/workflows/
│   ├── build.yml                # Standard CI workflow
│   └── copilot.yml              # Copilot CI workflow (optional)
├── ENVIRONMENT_ISSUES.md        # Detailed problem analysis
├── COPILOT_INSTRUCTIONS.md      # This file
├── HMM_CALCULATION_GUIDE.md     # Beginner-friendly HMM explanation
└── TEST_RESULTS.md              # Test execution summary
```

## Key Points for Copilot Agents

1. **Always use `makefile.copilot`** for building tests in Copilot environment
2. **Always set `LD_LIBRARY_PATH`** before running tests
3. **HTSlib must be version 1.16+** (older versions lack required BCF functions)
4. **Boost and HTSlib locations differ** from standard developer environment
5. **Never modify `makefile_common.mk`** or main makefiles - use Copilot-specific files

## Expected Test Results

All tests should pass:
```
Testing REAL supersite emission computation...
  Test 1a: Sample=REF, expect matches at haps 0,3...
    OK: emissions = [1.0, 0.01, 0.01, 1.0]
  ...
All tests passed!

Testing supersite accessor functions...
  SuperSite structure: OK
  ...
All tests passed!

Testing supersite code unpacking...
  Single code unpacking: OK
  ...
All tests passed!

Testing supersite builder...
  Building supersites from variant map...
  ...
All tests passed!

Testing missing multiallelic site imputation (Phase 3)...
  Creating 3-split multiallelic site...
  ...
All tests passed!

Testing buildSuperSites behavior...
  Building supersites from 10-variant context...
  ...
✓ SUCCESS: buildSuperSites correctly distinguishes biallelic vs multiallelic
All tests passed!
```

## Next Steps After Setup

1. All tests pass? ✅ Environment configured correctly
2. Want to add new tests? Copy existing test file patterns
3. Need to debug? Use `gdb` with test binaries
4. Want to profile? Use `valgrind` or `perf`

## Contact

For issues with this setup, check:
- ENVIRONMENT_ISSUES.md (detailed problem analysis)
- TEST_RESULTS.md (what tests should output)
- HMM_CALCULATION_GUIDE.md (understanding HMM calculations)

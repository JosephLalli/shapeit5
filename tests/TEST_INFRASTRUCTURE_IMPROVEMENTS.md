# Test Infrastructure Improvements

## Summary

The SHAPEIT5 test suite has been enhanced with standardized reporting, better crash diagnostics, and improved result visualization. These improvements make it much easier to understand test results, diagnose failures, and track test outcomes over time.

## What Was Changed

### 1. Standardized Test Reporting (`test_reporting.h`)

**Location:** `tests/src/test_reporting.h`

A new header provides consistent test output formatting and automatic crash diagnostics:

**Features:**
- Standard `[TEST_PASS]` / `[TEST_FAIL]` format for machine-readable parsing
- Automatic signal handlers that report context when tests crash
- Built-in test summary generation
- Scenario and iteration tracking for complex tests

**Usage Example:**

```cpp
#include "test_reporting.h"

int main() {
    TEST_INIT("my_test_name");  // Sets up signal handlers

    TEST_START("test_case_1", "Description");
    // ... test code ...
    TEST_CHECK(condition, "test_case_1", "Failure reason if false");

    // For multi-scenario tests
    SCENARIO_START("scenario_x8", "8x expansion test");
    ITERATION(3, 15, "burn3");  // Iteration 3 of 15
    // ... test code ...
    SCENARIO_PASS("scenario_x8", "15 iterations");

    TEST_SUMMARY();  // Prints summary
    return 0;
}
```

### 2. Enhanced `run_tests.sh`

**Improvements:**

#### Better Output Parsing
- Prioritizes new `[TEST_PASS]` / `[TEST_FAIL]` format
- Maintains backward compatibility with legacy output formats
- Extracts scenario information from `[SCENARIO_PASS]` / `[SCENARIO_FAIL]`
- Captures crash context from `[TEST_CRASH]` markers

#### Failure Context Collection
When tests fail, the script now automatically captures:
- Last 20 lines of output before failure
- Crash context (what the test was doing when it crashed)
- Active `SHAPEIT5_*` environment variables
- Exit code and failure type classification

Example output:
```
[FAILURE_CONTEXT] Crash: Signal: SIGSEGV (segmentation fault); Context: scenario 8, iteration 9
[FAILURE_CONTEXT] Last output before failure:
    [ITERATION] 9/15: burn7
    Anchor mismatch at locus 0: bial (1|0) vs supersite (0|1)
[FAILURE_CONTEXT] Environment:
    SHAPEIT5_TEST_TRACE=1
```

#### Formatted Summary Table
At the end of each test run, a table shows all test results:

```
Detailed Test Results:
┌────────────────────────────────────────┬──────────┬──────────┬─────────────┬───────────────────┐
│ Test Binary                            │ Category │ Status   │ Duration    │ Details           │
├────────────────────────────────────────┼──────────┼──────────┼─────────────┼───────────────────┤
│ test_supersite_builder                 │ unit     │ ✓ PASS   │  0.045s     │ 6/6               │
│ test_expansion_epochs                  │ integr.. │ ✗ FAIL   │  0.143s     │ ✓ rf_1 ✓ rf_2 ... │
│ test_supersite_hmm                     │ unit     │ ✗ FAIL   │  0.059s     │ [crash: SIGSEGV]  │
└────────────────────────────────────────┴──────────┴──────────┴─────────────┴───────────────────┘

By Category:
  unit:         18/25 passed
  integration:  4/8 passed
  smoke:        2/2 passed
  concurrency:  3/3 passed
```

#### Automatic Test Categorization
Tests are categorized based on naming patterns:
- `*parity*`, `*expansion*` → `integration`
- `*smoke*`, `*perf*` → `smoke`
- `*threading*` → `concurrency`
- Everything else → `unit`

### 3. Updated Test Files

**Converted to new format:**
- `test_supersite_builder.cpp` - Example of simple unit test conversion
- `test_supersite_expansion_epochs.cpp` - Example of complex multi-scenario test

**Benefits:**
- Clear scenario boundaries with `[SCENARIO_START]` / `[SCENARIO_PASS]`
- Iteration tracking: `[ITERATION] 3/15: burn3`
- Automatic crash context: if test crashes at iteration 9, you'll know exactly where
- Machine-readable output for potential CI/CD integration

## How to Use

### Running Tests

```bash
# Run all tests (same as before)
cd tests
./run_tests.sh

# Extract results from latest log without re-running
./run_tests.sh --extract-only

# Extract from specific log
./run_tests.sh --extract test_logs/test_results_*.log
```

### Updating Existing Tests

To convert an existing test to use the new format:

1. Include the header:
```cpp
#include "test_reporting.h"
```

2. Initialize at start of main():
```cpp
int main() {
    TEST_INIT("test_name");
    // ... rest of test
}
```

3. Replace manual pass/fail output:
```cpp
// OLD:
std::cout << "Test 1: OK" << std::endl;

// NEW:
TEST_CHECK(condition, "test_1", "Optional failure message");
// or
TEST_START("test_1", "Description");
// ... test code ...
TEST_PASS("test_1");  // or TEST_FAIL("test_1", "reason")
```

4. Add summary at end:
```cpp
TEST_SUMMARY();
return 0;
```

### For Multi-Scenario Tests

```cpp
for (const auto& scenario : scenarios) {
    SCENARIO_START("scenario_name", "description");
    TEST_CONTEXT("scenario_name initialization");  // For crash diagnostics

    for (int i = 0; i < iterations; i++) {
        ITERATION(i+1, iterations, "stage_name");
        TEST_CONTEXT("scenario_name, iteration " + std::to_string(i));
        // ... test code ...
    }

    SCENARIO_PASS("scenario_name", "details");
    // or on failure:
    SCENARIO_FAIL("scenario_name", "failure reason");
}
```

## Benefits

### Before
```
test_supersite_hmm: FAIL (0.059s, exit_code=139, type=segfault)
  Binary: test_supersite_hmm (exit_code: 139)
  Output:
    Testing supersite HMM harness...
      Building supersites
      Running single-precision forward
    timeout: the monitored command dumped core
```

### After
```
test_supersite_hmm: FAIL (0.059s, exit_code=139, type=segfault)
  Binary: test_supersite_hmm (exit_code: 139)
  [FAILURE_CONTEXT] Crash: Signal: SIGSEGV; Binary: test_supersite_hmm; Context: scenario_x8, iteration 9/15 (burn7)
  [FAILURE_CONTEXT] Last output before failure:
    [SCENARIO_START] scenario_x8: 8x expansion test
    [ITERATION] 9/15: burn7
    Anchor mismatch at locus 0: bial (1|0) vs supersite (0|1)
```

**Key improvements:**
1. **Know exactly where it crashed:** "scenario_x8, iteration 9/15 (burn7)"
2. **See what happened before crash:** Last 20 lines of output
3. **Understand failure mode:** "Anchor mismatch at locus 0"
4. **Track by category:** See that integration tests are failing more than unit tests

## Migration Path

**Phase 1 (Current):**
- Infrastructure in place
- Two reference tests converted
- run_tests.sh handles both old and new formats

**Phase 2 (Optional):**
- Gradually convert remaining tests
- Tests can be converted one at a time
- No breaking changes to existing tests

**Phase 3 (Future):**
- Potential CI/CD integration using standardized format
- Test filtering by category (`./run_tests.sh --filter=unit`)
- Verbosity levels (`./run_tests.sh --verbose=1`)

## Files Changed

- `tests/src/test_reporting.h` - NEW: Standardized reporting utilities
- `tests/run_tests.sh` - ENHANCED: Better parsing, failure context, summary table
- `tests/src/test_supersite_builder.cpp` - CONVERTED: Example unit test
- `tests/src/test_supersite_expansion_epochs.cpp` - CONVERTED: Example integration test

## Backward Compatibility

All existing tests continue to work without modification. The enhanced `run_tests.sh` still parses old output formats (✓, "OK", "passed", etc.) while preferring the new standardized format when present.

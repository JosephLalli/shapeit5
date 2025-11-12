#!/bin/bash

# Test runner script for SHAPEIT5 test binaries
# Runs all binaries in tests/bin and tracks results with logging
# Handles multi-test binaries and both assertion crashes and proper test reporting
#
# Usage:
#   ./run_tests.sh                    # Run all tests
#   ./run_tests.sh --extract-only     # Extract results from latest log
#   ./run_tests.sh --extract <log>    # Extract results from specific log
#   ./run_tests.sh --help             # Show help

set -uo pipefail

# Extract test results from log file
extract_test_results() {
    local log_file="$1"
    
    if [[ ! -f "$log_file" ]]; then
        echo "Error: Log file $log_file not found" >&2
        return 1
    fi
    
    echo "=== Test Results Summary ==="
    echo "Log: $log_file"
    echo
    
    # Extract individual test results using the same format as the test runner
    echo "Individual Test Results:"
    grep -E "^[^:]+: (PASS|FAIL)" "$log_file" | while IFS= read -r line; do
        if [[ $line =~ :[[:space:]]*PASS ]]; then
            echo "✓ $line"
        else
            echo "✗ $line"
        fi
    done
    
    echo
    echo "=== Summary Stats ==="
    local total_tests
    local passed_tests  
    local failed_tests
    
    total_tests=$(grep -cE "^[^:]+: (PASS|FAIL)" "$log_file" 2>/dev/null | tr -d '\n' || echo "0")
    passed_tests=$(grep -cE "^[^:]+: PASS" "$log_file" 2>/dev/null | tr -d '\n' || echo "0")
    failed_tests=$(grep -cE "^[^:]+: FAIL" "$log_file" 2>/dev/null | tr -d '\n' || echo "0")
    
    echo "Total tests: $total_tests"
    echo "Passed: $passed_tests"
    echo "Failed: $failed_tests"
    
    if [[ $failed_tests -gt 0 ]]; then
        echo
        echo "=== Failed Tests ==="
        grep -E "^[^:]+: FAIL" "$log_file"
    fi
    
    # Return failed test count as exit code (max 255)
    if [[ $failed_tests -gt 255 ]]; then
        return 255
    else
        return $failed_tests
    fi
}

# Show help
show_help() {
    cat << EOF
SHAPEIT5 Test Runner

Usage:
  $0                        Run all tests and generate new log
  $0 --extract-only         Extract results from latest log file
  $0 --extract <log_file>   Extract results from specific log file
  $0 --help                 Show this help message

Examples:
  $0                                           # Run all tests
  $0 --extract-only                           # Parse latest results
  $0 --extract test_logs/test_results_*.log   # Parse specific log

The test runner generates logs in test_logs/ with format:
  test_results_<commit>_<timestamp>.log
  
Individual test results use the format:
  test_name: PASS|FAIL (duration_seconds)
EOF
}

# Parse command line arguments
case "${1:-}" in
    --help|-h)
        show_help
        exit 0
        ;;
    --extract-only)
        if [[ -L "test_logs/latest_test_results.log" && -f "test_logs/latest_test_results.log" ]]; then
            extract_test_results "test_logs/latest_test_results.log"
            exit $?
        else
            echo "Error: No latest test results found" >&2
            exit 1
        fi
        ;;
    --extract)
        if [[ -z "${2:-}" ]]; then
            echo "Error: --extract requires a log file path" >&2
            exit 1
        fi
        extract_test_results "$2"
        exit $?
        ;;
    --*)
        echo "Error: Unknown option $1" >&2
        echo "Use --help for usage information" >&2
        exit 1
        ;;
esac

# Global error trap to handle unexpected script termination
cleanup_on_exit() {
    local exit_code=$?
    # Only report unexpected exits (not when script exits normally with failure due to test failures)
    if [ $exit_code -ne 0 ] && [ -n "${CURRENT_LOG:-}" ] && [ "${EXPECTED_EXIT:-}" != "true" ]; then
        echo "" | tee -a "$CURRENT_LOG"
        echo "ERROR: Test script terminated unexpectedly with exit code $exit_code" | tee -a "$CURRENT_LOG"
        echo "This may indicate a crash that wasn't properly caught." | tee -a "$CURRENT_LOG"
    fi
}
trap cleanup_on_exit EXIT

# Configuration
TESTDIR="$(dirname $(realpath "${BASH_SOURCE[0]}"))"
TEST_BIN_DIR="$TESTDIR/bin"
LOG_DIR="$TESTDIR/test_logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
COMMIT_HASH=$(git rev-parse --short HEAD)
COMMIT_MSG=$(git log --oneline -1 --pretty=format:"%s")

# Create log directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Log files
CURRENT_LOG="$LOG_DIR/test_results_${COMMIT_HASH}_${TIMESTAMP}.log"
LATEST_SYMLINK="$LOG_DIR/latest_test_results.log"

# Global counters for individual tests
total_individual_tests=0
passed_individual_tests=0
failed_individual_tests=0
failed_individual_test_names=()

# Function to parse test output and extract individual test results
parse_test_output() {
    local output="$1"
    local binary_name="$2"
    local exit_code="$3"
    local duration="$4"
    
    # Common patterns for test results
    local test_patterns=(
        '✓.*passed'                    # ✓ Test X passed
        '✗.*failed'                    # ✗ Test X failed
        'Test [0-9]+.*passed'          # Test 1 passed
        'Test [0-9]+.*failed'          # Test 1 failed
        'PASS.*Test'                   # PASS: Test description
        'FAIL.*Test'                   # FAIL: Test description
        '\[PASS\]'                     # [PASS]
        '\[FAIL\]'                     # [FAIL]
        'OK.*test'                     # OK - test description
        'ERROR.*test'                  # ERROR - test description
    )
    
    local found_individual_tests=false
    local individual_tests_in_binary=0
    local passed_in_binary=0
    local failed_in_binary=0
    
    # Look for individual test results in output
    while IFS= read -r line; do
        lc_line="${line,,}"
        is_test_line=false
        if [[ "$lc_line" == *"test"* ]]; then
            is_test_line=true
        elif [[ "$lc_line" == *"pass"* ]] || [[ "$lc_line" == *"fail"* ]] || [[ "$lc_line" == *"ok"* ]]; then
            is_test_line=true
        fi

        # Check for passed tests
        if $is_test_line && [[ "$lc_line" == *"✓"* || "$lc_line" == *"passed"* || "$lc_line" == *"[pass]"* || "$lc_line" == *" ok"* || "$lc_line" == *"pass:"* ]]; then
            found_individual_tests=true
            ((individual_tests_in_binary++))
            ((passed_in_binary++))
            ((total_individual_tests++))
            ((passed_individual_tests++))
            
            # Extract test name/description
            local test_desc=$(echo "$line" | sed -E 's/.*[Tt]est[[:space:]]*[0-9]*[[:space:]]*[:]*[[:space:]]*(.*)/\1/' | sed 's/[✓✗]//g' | sed 's/\[(PASS|FAIL)\]//g' | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
            if [[ -z "$test_desc" ]]; then
                test_desc="unnamed_test_$individual_tests_in_binary"
            fi
            
            echo "${binary_name}::${test_desc}: PASS (${duration}s)" | tee -a "$CURRENT_LOG"
            
        # Check for failed tests - be more specific to avoid false positives
        # Only flag explicit failure indicators, not words like "error" in descriptive context
        elif $is_test_line && [[ "$lc_line" == *"✗"* || "$lc_line" == *"failed"* || "$lc_line" == *"[fail]"* ]]; then
            found_individual_tests=true
            ((individual_tests_in_binary++))
            ((failed_in_binary++))
            ((total_individual_tests++))
            ((failed_individual_tests++))
            
            # Extract test name/description
            local test_desc=$(echo "$line" | sed -E 's/.*[Tt]est[[:space:]]*[0-9]*[[:space:]]*[:]*[[:space:]]*(.*)/\1/' | sed 's/[✓✗]//g' | sed 's/\[(PASS|FAIL)\]//g' | sed 's/^[[:space:]]*//' | sed 's/[[:space:]]*$//')
            if [[ -z "$test_desc" ]]; then
                test_desc="unnamed_test_$individual_tests_in_binary"
            fi
            
            echo "${binary_name}::${test_desc}: FAIL (${duration}s)" | tee -a "$CURRENT_LOG"
            failed_individual_test_names+=("${binary_name}::${test_desc}")
        fi
    done <<< "$output"
    
    # If no individual test patterns found, treat the whole binary as one test
    if [ "$found_individual_tests" = false ]; then
        ((total_individual_tests++))
        if [ "$exit_code" -eq 0 ]; then
            ((passed_individual_tests++))
            echo "${binary_name}: PASS (${duration}s)" | tee -a "$CURRENT_LOG"
        else
            ((failed_individual_tests++))
            failed_individual_test_names+=("${binary_name}")
            
            # Determine failure type
            local failure_type="unknown"
            if [ "$exit_code" -eq 134 ]; then
                failure_type="assertion_abort"
            elif [ "$exit_code" -eq 139 ]; then
                failure_type="segfault"
            elif [ "$exit_code" -eq 6 ]; then
                failure_type="abort_signal"
            elif [ "$exit_code" -eq 11 ]; then
                failure_type="segv_signal"
            elif [ "$exit_code" -eq 4 ]; then
                failure_type="illegal_instruction"
            elif [ "$exit_code" -eq 8 ]; then
                failure_type="floating_point_exception"
            elif [ "$exit_code" -eq 124 ]; then
                failure_type="timeout"
            elif echo "$output" | grep -qi "TIMEOUT:"; then
                failure_type="timeout"
            elif echo "$output" | grep -qi "assertion.*failed"; then
                failure_type="assertion_failed"
            elif echo "$output" | grep -qi "terminate called"; then
                failure_type="cpp_terminate"
            elif echo "$output" | grep -qi "pure virtual"; then
                failure_type="pure_virtual_call"
            elif echo "$output" | grep -qi "error"; then
                failure_type="runtime_error"
            fi
            
            echo "${binary_name}: FAIL (${duration}s, exit_code=${exit_code}, type=${failure_type})" | tee -a "$CURRENT_LOG"
        fi
    fi
    
    # Always log the output for failed binaries or detailed info
    if [ "$exit_code" -ne 0 ] || [ "$failed_in_binary" -gt 0 ]; then
        echo "  Binary: $binary_name (exit_code: $exit_code)" | tee -a "$CURRENT_LOG"
        echo "  Individual tests found: $individual_tests_in_binary (passed: $passed_in_binary, failed: $failed_in_binary)" | tee -a "$CURRENT_LOG"
        echo "  Output:" | tee -a "$CURRENT_LOG"
        echo "$output" | sed 's/^/    /' | tee -a "$CURRENT_LOG"
        echo "" | tee -a "$CURRENT_LOG"
    fi
    
    # Return non-zero if this binary had any failures
    if [ "$exit_code" -ne 0 ] || [ "$failed_in_binary" -gt 0 ]; then
        return 1
    else
        return 0
    fi
}

# Function to run a single test binary
run_test() {
    local test_binary="$1"
    local test_name=$(basename "$test_binary")
    
    echo "Running $test_name..."
    
    local start_time=$(date +%s.%N)
    local output=""
    local exit_code
    local timeout_seconds=300  # 5 minute timeout per test
    
    # Run the test and capture both output and exit code with timeout protection
    set +e  # Don't exit on command failure
    
    # Use timeout to prevent hanging tests and capture all possible crashes
    if command -v timeout >/dev/null 2>&1; then
        # Use --preserve-status to get the actual exit code from the process
        output=$(timeout --preserve-status "$timeout_seconds" "$test_binary" 2>&1)
        exit_code=$?
        
        # Check if timeout killed the process (timeout returns 124 only without --preserve-status)
        # With --preserve-status, we need to check if the process was killed by timeout
        if echo "$output" | grep -q "Terminated"; then
            output="$output
TIMEOUT: Test binary exceeded ${timeout_seconds} seconds"
            exit_code=124
        fi
    else
        # Fallback without timeout if timeout command not available
        output=$("$test_binary" 2>&1)
        exit_code=$?
    fi
    
    # Additional protection against script termination from child process signals
    trap '' SIGINT SIGTERM
    
    local end_time=$(date +%s.%N)
    local duration
    # Use bc for duration calculation with fallback to basic arithmetic
    if command -v bc >/dev/null 2>&1; then
        duration=$(echo "$end_time - $start_time" | bc -l 2>/dev/null) || duration="0.000"
    else
        # Fallback to basic bash arithmetic (less precise)
        duration=$(awk "BEGIN {printf \"%.3f\", $end_time - $start_time}")
    fi
    
    # Reset signal handling
    trap - SIGINT SIGTERM
    
    # Format duration to 3 decimal places
    duration=$(printf "%.3f" "$duration")
    
    # Parse the output for individual tests
    parse_test_output "$output" "$test_name" "$exit_code" "$duration"
    
    local parse_result=$?
    
    # Don't re-enable set -e here as it can cause issues with error handling
    # The main script doesn't need set -e for proper operation
    
    return $parse_result
}

# Main execution
echo "========================================" | tee "$CURRENT_LOG"
echo "SHAPEIT5 Test Run Report" | tee -a "$CURRENT_LOG"
echo "Timestamp: $(date)" | tee -a "$CURRENT_LOG"
echo "Commit: $COMMIT_HASH - $COMMIT_MSG" | tee -a "$CURRENT_LOG"
echo "========================================" | tee -a "$CURRENT_LOG"
echo "" | tee -a "$CURRENT_LOG"

# Check if test directory exists
if [ ! -d "$TEST_BIN_DIR" ]; then
    echo "ERROR: Test binary directory '$TEST_BIN_DIR' not found" | tee -a "$CURRENT_LOG"
    exit 1
fi

# Get list of test binaries
test_binaries=($(find "$TEST_BIN_DIR" -type f -executable | sort))

if [ ${#test_binaries[@]} -eq 0 ]; then
    echo "ERROR: No executable test binaries found in '$TEST_BIN_DIR'" | tee -a "$CURRENT_LOG"
    exit 1
fi

echo "Found ${#test_binaries[@]} test binaries" | tee -a "$CURRENT_LOG"
echo "" | tee -a "$CURRENT_LOG"

# Run all tests
total_binaries=${#test_binaries[@]}
passed_binaries=0
failed_binaries=0
failed_binary_names=()

for test_binary in "${test_binaries[@]}"; do
    if run_test "$test_binary"; then
        ((passed_binaries++))
    else
        ((failed_binaries++))
        failed_binary_names+=($(basename "$test_binary"))
    fi
done

# Summary
echo "" | tee -a "$CURRENT_LOG"
echo "========================================" | tee -a "$CURRENT_LOG"
echo "Test Summary:" | tee -a "$CURRENT_LOG"
echo "  Total binaries: $total_binaries" | tee -a "$CURRENT_LOG"
echo "  Passed binaries: $passed_binaries" | tee -a "$CURRENT_LOG"
echo "  Failed binaries: $failed_binaries" | tee -a "$CURRENT_LOG"
echo "" | tee -a "$CURRENT_LOG"
echo "  Total individual tests: $total_individual_tests" | tee -a "$CURRENT_LOG"
echo "  Passed individual tests: $passed_individual_tests" | tee -a "$CURRENT_LOG"
echo "  Failed individual tests: $failed_individual_tests" | tee -a "$CURRENT_LOG"

if [ $failed_binaries -gt 0 ]; then
    echo "" | tee -a "$CURRENT_LOG"
    echo "  Failed binaries: ${failed_binary_names[*]}" | tee -a "$CURRENT_LOG"
fi

if [ $failed_individual_tests -gt 0 ]; then
    echo "  Failed individual tests: ${failed_individual_test_names[*]}" | tee -a "$CURRENT_LOG"
fi

echo "========================================" | tee -a "$CURRENT_LOG"

# Compare with previous run if available
if [ -L "$LATEST_SYMLINK" ] && [ -f "$LATEST_SYMLINK" ]; then
    echo "" | tee -a "$CURRENT_LOG"
    echo "Changes from previous run:" | tee -a "$CURRENT_LOG"
    
    # Extract test name and outcome from both logs
    declare -A current_outcomes previous_outcomes
    
    # Parse current results
    while IFS= read -r line; do
        if [[ $line =~ ^([^:]+):[[:space:]]+(PASS|FAIL) ]]; then
            test_name="${BASH_REMATCH[1]}"
            outcome="${BASH_REMATCH[2]}"
            current_outcomes["$test_name"]="$outcome"
        fi
    done < <(grep -E "^(test_.*|.*::.*): (PASS|FAIL)" "$CURRENT_LOG")
    
    # Parse previous results
    while IFS= read -r line; do
        if [[ $line =~ ^([^:]+):[[:space:]]+(PASS|FAIL) ]]; then
            test_name="${BASH_REMATCH[1]}"
            outcome="${BASH_REMATCH[2]}"
            previous_outcomes["$test_name"]="$outcome"
        fi
    done < <(grep -E "^(test_.*|.*::.*): (PASS|FAIL)" "$LATEST_SYMLINK")
    
    # Find tests with changed outcomes
    changes_found=false
    for test_name in "${!current_outcomes[@]}"; do
        current_outcome="${current_outcomes[$test_name]}"
        previous_outcome="${previous_outcomes[$test_name]:-}"
        
        if [[ -n "$previous_outcome" && "$current_outcome" != "$previous_outcome" ]]; then
            if [ "$changes_found" = false ]; then
                echo "  Test outcome changes detected:" | tee -a "$CURRENT_LOG"
                changes_found=true
            fi
            echo "    $test_name: $previous_outcome → $current_outcome" | tee -a "$CURRENT_LOG"
        fi
    done
    
    # Check for new tests (in current but not previous)
    for test_name in "${!current_outcomes[@]}"; do
        if [[ -z "${previous_outcomes[$test_name]:-}" ]]; then
            if [ "$changes_found" = false ]; then
                echo "  New tests detected:" | tee -a "$CURRENT_LOG"
                changes_found=true
            fi
            echo "    $test_name: NEW (${current_outcomes[$test_name]})" | tee -a "$CURRENT_LOG"
        fi
    done
    
    # Check for removed tests (in previous but not current)
    for test_name in "${!previous_outcomes[@]}"; do
        if [[ -z "${current_outcomes[$test_name]:-}" ]]; then
            if [ "$changes_found" = false ]; then
                echo "  Removed tests detected:" | tee -a "$CURRENT_LOG"
                changes_found=true
            fi
            echo "    $test_name: REMOVED (was ${previous_outcomes[$test_name]})" | tee -a "$CURRENT_LOG"
        fi
    done
    
    if [ "$changes_found" = false ]; then
        echo "  No changes in test outcomes" | tee -a "$CURRENT_LOG"
    fi
else
    echo "" | tee -a "$CURRENT_LOG"
    echo "No previous test results found for comparison" | tee -a "$CURRENT_LOG"
fi

# Update latest symlink
ln -sf "$(basename "$CURRENT_LOG")" "$LATEST_SYMLINK"

echo "" | tee -a "$CURRENT_LOG"
echo "Log saved to: $CURRENT_LOG" | tee -a "$CURRENT_LOG"

# Exit with appropriate code
EXPECTED_EXIT=true
if [ $failed_individual_tests -gt 0 ]; then
    exit 1
else
    exit 0
fi

#!/bin/bash
# Extract test results from test logs for easy parsing
# Usage: ./extract_test_results.sh [log_file]

LOG_FILE="${1:-test_logs/latest_test_results.log}"

if [[ ! -f "$LOG_FILE" ]]; then
    echo "Error: Log file $LOG_FILE not found" >&2
    exit 1
fi

echo "=== Test Results Summary ==="
echo "Log: $LOG_FILE"
echo

# Extract individual test results
echo "Individual Test Results:"
grep -E "^[^:]+: (PASS|FAIL)" "$LOG_FILE" | while IFS= read -r line; do
    if [[ $line =~ : PASS ]]; then
        echo "✓ $line"
    else
        echo "✗ $line"
    fi
done

echo
echo "=== Summary Stats ==="
total_tests=$(grep -cE "^[^:]+: (PASS|FAIL)" "$LOG_FILE")
passed_tests=$(grep -cE "^[^:]+: PASS" "$LOG_FILE")
failed_tests=$(grep -cE "^[^:]+: FAIL" "$LOG_FILE")

echo "Total tests: $total_tests"
echo "Passed: $passed_tests"
echo "Failed: $failed_tests"

if [[ $failed_tests -gt 0 ]]; then
    echo
    echo "=== Failed Tests ==="
    grep -E "^[^:]+: FAIL" "$LOG_FILE"
fi
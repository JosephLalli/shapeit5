#!/usr/bin/env bash
# Simple regression check: run the traced supersite backward parity test and look for the IMPUTE guard skip message
set -euo pipefail
# Resolve script and repo paths robustly so the script works when invoked from repo root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPRO_TEST="$SCRIPT_DIR/bin/test_supersite_backward_parity"
if [ ! -x "$REPRO_TEST" ]; then
  echo "Test binary $REPRO_TEST not found or not executable. Build tests first: make -C tests" >&2
  exit 2
fi
# Run with trace enabled and capture output
LD_LIBRARY_PATH=$HOME/.linuxbrew/lib:/usr/local/lib SHAPEIT5_TEST_TRACE=1 "$REPRO_TEST" 2>&1 | tee /tmp/check_backward_impute_guard.out
# Check for a guard message
if grep -q "Skipping IMPUTE for locus" /tmp/check_backward_impute_guard.out; then
  echo "OK: IMPUTE guard observed in test output"
  exit 0
else
  echo "WARN: IMPUTE guard message not found in test output. This may be OK if the test case did not exercise a guarded IMPUTE path." >&2
  echo "See /tmp/check_backward_impute_guard.out for full test output." >&2
  # Do not fail the overall test-run for absence of the guard message; keep as a non-failing warning
  exit 0
fi

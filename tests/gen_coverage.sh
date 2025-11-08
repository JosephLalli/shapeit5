#!/usr/bin/env bash
# tests/gen_coverage.sh
# Makefile-based coverage for SHAPEIT5.
# Produces tests/coverage/lcov.info readable by VS Code Coverage Gutters.

set -euo pipefail

# --- Config ----
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MAKE_JOBS="${MAKE_JOBS:-24}"
COV_DIR="${REPO_ROOT}/tests/coverage"
# Globs to exclude from the report
EXCLUDE_GLOBS=(
  "/usr/*"
  "*/.git/*"
  "*/third_party/*"
  "*/external/*"
  "*/test/*"
  "*/tests/*"
  "*/build/*"
)

mkdir -p "${COV_DIR}"

echo "=== SHAPEIT5 Coverage Generation ==="
echo "Repository root: ${REPO_ROOT}"
echo "Coverage output: ${COV_DIR}"
echo "Make jobs: ${MAKE_JOBS}"
echo ""

detect_cc() {
  local bin="${CXX:-c++}"
  if "$bin" --version 2>/dev/null | grep -qi clang; then 
    echo "clang"
  else 
    echo "gcc"
  fi
}

# Clean previous coverage artifacts so the run is deterministic
preclean() {
  echo ">>> Cleaning previous coverage artifacts..."
  cd "${REPO_ROOT}"
  # GCC artifacts
  find . -name "*.gcno" -o -name "*.gcda" -delete 2>/dev/null || true
  # Clang artifacts
  find . -name "*.profraw" -delete 2>/dev/null || true
  rm -f "${COV_DIR}/lcov.info" "${COV_DIR}/coverage.profdata" "${COV_DIR}/lcov.raw.info" "${COV_DIR}/lcov.tmp.info"
  echo ">>> Preclean complete."
}

build_with_make() {
  local cc="$1"
  echo ">>> Building with ${cc} and coverage instrumentation..."
  cd "${REPO_ROOT}"
  
  # Clean and build
  echo ""
  echo ">>> Running: make clean"
  echo "========================================"
  if ! make clean 2>&1; then
    echo "========================================"
    echo "WARNING: make clean failed, continuing anyway..." >&2
  fi
  echo "========================================"
  echo ""
  
  echo ">>> Building main binaries with -j${MAKE_JOBS}..."
  echo "========================================"
  # Inject flags using the variable names used by this repo (CXXFLAG/LDFLAG)
  if ! make coverage -j"${MAKE_JOBS}" 2>&1; then
    echo "========================================"
    echo ""
    echo "ERROR: Main binary build failed!" >&2
    return 1
  fi
  echo "========================================"
  echo ""
  
  # Quick compile-time instrumentation sanity check
  local gcno_count
  gcno_count=$(find . -name "*.gcno" | wc -l | tr -d ' \t\n')
  echo ">>> Compile-time coverage artifacts (.gcno) found: ${gcno_count}"
  if [[ "${gcno_count:-0}" -eq 0 ]]; then
    echo "WARNING: No .gcno files found after build. Coverage flags may not have been applied." >&2
    echo "         Ensure the makefiles honor CXXFLAG/LDFLAG or pass them on the command line." >&2
  fi
  
  echo ">>> Build complete."
}

run_tests() {
  echo ">>> Running tests..."
  cd "${REPO_ROOT}"
  
  # Set library path per SHAPEIT5 requirements
  export LD_LIBRARY_PATH="${HOME}/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH:-}"
  
  # For clang, direct profile files to coverage directory
  if [[ "$(detect_cc)" == "clang" ]]; then
    export LLVM_PROFILE_FILE="${COV_DIR}/default-%p-%m.profraw"
  fi
  
  # Run the test script
  local test_runner="${REPO_ROOT}/tests/run_tests.sh"
  if [[ ! -x "$test_runner" ]]; then
    echo "ERROR: Test runner not found or not executable: $test_runner" >&2
    return 1
  fi
  
  echo ""
  echo ">>> Executing: $test_runner"
  echo "========================================"
  
  # Run tests - capture exit code but continue even on failure
  local test_exit_code=0
  bash "$test_runner" 2>&1 || test_exit_code=$?
  
  echo "========================================"
  echo ""
  
  if [[ "$test_exit_code" -ne 0 ]]; then
    echo "WARNING: Some tests failed (exit code: ${test_exit_code}), but continuing with coverage generation..." >&2
  fi
  
  echo ">>> Tests complete."
}

emit_lcov_gcc() {
  echo ">>> Generating LCOV report (GCC)..."
  
  if ! command -v lcov >/dev/null; then
    echo "ERROR: lcov not found. Install with: sudo apt-get install lcov" >&2
    exit 1
  fi
  
  local raw="${COV_DIR}/lcov.raw.info"
  local out="${COV_DIR}/lcov.info"
  rm -f "$raw" "$out"

  cd "${REPO_ROOT}"
  
  # Check for .gcda files
  local gcda_count=$(find . -name "*.gcda" | wc -l)
  echo ">>> Found ${gcda_count} .gcda files"
  
  if [[ "$gcda_count" -eq 0 ]]; then
    echo "ERROR: No .gcda files found. Coverage instrumentation may have failed." >&2
    echo ">>> Checking for .gcno files (compile-time coverage data):" >&2
    find . -name "*.gcno" | head -5 >&2
    exit 1
  fi

  # Capture from current tree (object + gcda scattered in-source)
  echo ">>> Running lcov --capture..."
  if ! lcov --directory . --capture --output-file "$raw" --ignore-errors inconsistent,mismatch --rc lcov_branch_coverage=1 2>&1 | tail -20; then
    echo "ERROR: lcov capture failed!" >&2
    exit 1
  fi
  
  if [[ ! -f "$raw" ]] || [[ ! -s "$raw" ]]; then
    echo "ERROR: lcov.raw.info was not created or is empty" >&2
    exit 1
  fi
  
  echo ">>> Raw LCOV file size: $(wc -l < "$raw") lines"

  # Apply exclusions
  echo ">>> Applying exclusions..."
  local filtered="$raw"
  for g in "${EXCLUDE_GLOBS[@]}"; do
    local tmp="${COV_DIR}/lcov.tmp.info"
    lcov --remove "$filtered" "$g" --output-file "$tmp" --ignore-errors inconsistent,mismatch,unused --rc lcov_branch_coverage=1 #2>&1 | tail -5
    filtered="$tmp"
  done
  echo ">>> LCOV tmp file created"
  mv "$filtered" "$out"
  rm -f "$raw" "${COV_DIR}/lcov.tmp.info"
  
  local line_count=$(grep -c '^DA:' "$out" || echo 0)
  echo ">>> LCOV written: $out"
  echo ">>> Lines covered: $line_count"
}

emit_lcov_clang() {
  echo ">>> Generating LCOV report (Clang)..."
  command -v llvm-profdata >/dev/null || { echo "ERROR: llvm-profdata not found" >&2; exit 1; }
  command -v llvm-cov >/dev/null || { echo "ERROR: llvm-cov not found" >&2; exit 1; }

  local profdata="${COV_DIR}/coverage.profdata"
  local out="${COV_DIR}/lcov.info"
  rm -f "$profdata" "$out"

  # Merge all .profraw emitted during tests
  mapfile -t raws < <(find "${COV_DIR}" -name "*.profraw" -type f)
  if (( ${#raws[@]} == 0 )); then
    echo "ERROR: No .profraw files found under ${COV_DIR}." >&2
    echo "Check that LLVM_PROFILE_FILE was set correctly during test execution." >&2
    exit 1
  fi
  
  echo ">>> Found ${#raws[@]} .profraw files, merging..."
  llvm-profdata merge -sparse "${raws[@]}" -o "$profdata"

  # Gather candidate binaries/libraries for coverage attribution
  mapfile -t bins < <(find "${REPO_ROOT}/phase_common/bin" "${REPO_ROOT}/tests/bin" \
      -type f -executable 2>/dev/null || true)

  if (( ${#bins[@]} == 0 )); then
    echo "ERROR: No executables found to attribute coverage." >&2
    echo "Expected binaries in phase_common/bin and tests/bin" >&2
    exit 1
  fi
  
  echo ">>> Found ${#bins[@]} binaries for coverage attribution"

  # Build ignore regex from globs (join with |, translate * to .*)
  join_regex() {
    local r=()
    for g in "${EXCLUDE_GLOBS[@]}"; do
      r+=( "$(printf '%s' "$g" | sed -e 's/[].[^$\\/*]/\\&/g' -e 's#/#/#g' -e 's/\*/.*/g')" )
    done
    (IFS='|'; echo "${r[*]}")
  }

  llvm-cov export \
    --format=lcov \
    --ignore-filename-regex="$(join_regex)" \
    --instr-profile "$profdata" \
    "${bins[@]}" > "$out"

  if [[ ! -s "$out" ]]; then
    echo "ERROR: Generated LCOV file is empty." >&2
    exit 1
  fi

  echo ">>> LCOV written: $out"
  echo ">>> Lines covered: $(grep -c '^DA:' "$out" || echo 0)"
}

generate_html_report() {
    # ---- Generate HTML coverage report (optional but highly recommended) ----
  # Requires: genhtml (part of the lcov package)
  HTML_DIR="${COV_DIR}/html"
  rm -rf "${HTML_DIR}"

  echo ">>> Generating HTML coverage report with branch coverage..."
  genhtml "${COV_DIR}/lcov.info" \
    --output-directory "${HTML_DIR}" \
    --branch-coverage \
    --demangle-cpp \
    --title "SHAPEIT5 Test Coverage" \
    --legend \
    --show-details \
    --function-coverage \
    --quiet \
    --rc geninfo_unexecuted_blocks=1 \
    --ignore-errors inconsistent,mismatch,deprecated,usage

  if [[ -f "${HTML_DIR}/index.html" ]]; then
    echo "HTML coverage report generated:"
    echo "  file://${HTML_DIR}/index.html"
  else
    echo "WARNING: genhtml did not produce an index.html file."
  fi

}

main() {
#  preclean
  local cc; cc="$(detect_cc)"
#  build_with_make "$cc"
#  run_tests
  if [[ "$cc" == "clang" ]]; then
    emit_lcov_clang
  else
    emit_lcov_gcc
  fi
  generate_html_report
  echo "Done. Open ${COV_DIR}/lcov.info with Coverage Gutters."
}

main "$@"

#!/usr/bin/env bash
# tools/gen_coverage.sh
# Makefile-based coverage for C++ projects (no CMake).
# Produces coverage/lcov.info readable by VS Code Coverage Gutters.

set -euo pipefail

# --- Config ----
TEST_RUNNER="${TEST_RUNNER:-tests/run_tests.sh}"
MAKE_JOBS="${MAKE_JOBS:-24}"
COV_DIR="${COV_DIR:-coverage}"
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

detect_cc() {
  local bin="${CXX:-c++}"
  if "$bin" --version 2>/dev/null | grep -qi clang; then echo "clang"; else echo "gcc"; fi
}

# Clean previous coverage artifacts so the run is deterministic
preclean() {
  # GCC artifacts
  find . -name "*.gcno" -o -name "*.gcda" -delete 2>/dev/null || true
  # Clang artifacts
  find . -name "*.profraw" -delete 2>/dev/null || true
  rm -f "${COV_DIR}/lcov.info" "${COV_DIR}/coverage.profdata"
}

build_with_make() {
  local cc="$1"
  if [[ "$cc" == "clang" ]]; then
    # Clang: instrumentation + mapping
    export CXXFLAGS="${CXXFLAGS:-} -O0 -g -fprofile-instr-generate -fcoverage-mapping"
    export CFLAGS="${CFLAGS:-} -O0 -g -fprofile-instr-generate -fcoverage-mapping"
    export LDFLAGS="${LDFLAGS:-} -fprofile-instr-generate"
  else
    # GCC: gcov
    export CXXFLAGS="${CXXFLAGS:-} -O0 -g --coverage"
    export CFLAGS="${CFLAGS:-} -O0 -g --coverage"
    export LDFLAGS="${LDFLAGS:-} --coverage"
  fi

  # Usual SHAPEIT5 build pattern
  make clean || true
  make -j "${MAKE_JOBS}"
}

run_tests() {
  # For clang, direct profile files to a known place
  if [[ "$(detect_cc)" == "clang" ]]; then
    export LLVM_PROFILE_FILE="${COV_DIR}/default-%p-%m.profraw"
  fi
  bash "${TEST_RUNNER}"
}

emit_lcov_gcc() {
  command -v lcov >/dev/null || { echo "lcov not found"; exit 1; }
  local raw="${COV_DIR}/lcov.raw.info"
  local out="${COV_DIR}/lcov.info"
  rm -f "$raw" "$out"

  # Capture from current tree (object + gcda scattered in-source)
  lcov --directory . --capture --output-file "$raw"

  # Apply exclusions
  for g in "${EXCLUDE_GLOBS[@]}"; do
    lcov --remove "$raw" "$g" --output-file "$raw"
  done

  mv "$raw" "$out"
  echo "LCOV written: $out"
}

emit_lcov_clang() {
  command -v llvm-profdata >/dev/null || { echo "llvm-profdata not found"; exit 1; }
  command -v llvm-cov >/dev/null || { echo "llvm-cov not found"; exit 1; }

  local profdata="${COV_DIR}/coverage.profdata"
  local out="${COV_DIR}/lcov.info"
  rm -f "$profdata" "$out"

  # Merge all .profraw emitted during tests
  mapfile -t raws < <(find "${COV_DIR}" -name "*.profraw" -type f)
  if (( ${#raws[@]} == 0 )); then
    echo "No .profraw files found under ${COV_DIR}."
    exit 1
  fi
  llvm-profdata merge -sparse "${raws[@]}" -o "$profdata"

  # Gather candidate binaries/libraries for coverage attribution
  # (adjust or extend these paths to match your repo layout)
  mapfile -t bins < <(find . -type f -perm -111 \
    -not -path "*/.git/*" -not -path "*/static_bins/*" 2>/dev/null)

  if (( ${#bins[@]} == 0 )); then
    echo "No executables found to attribute coverage."
    exit 1
  fi

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

  echo "LCOV written: $out"
}

main() {
  preclean
  local cc; cc="$(detect_cc)"
  build_with_make "$cc"
  run_tests
  if [[ "$cc" == "clang" ]]; then
    emit_lcov_clang
  else
    emit_lcov_gcc
  fi
  echo "Done. Open coverage/lcov.info with Coverage Gutters."
}

main "$@"

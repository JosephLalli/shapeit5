#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

usage() {
  cat <<EOF
Usage: $0 [--update] [--list]

Runs a curated, deterministic integration test suite on PAR2.

Options:
  --update   Regenerate expected MD5/VCF.gz baselines (UPDATE_EXPECTED=1)
  --list     Print the scripts included in the suite and exit

Environment overrides:
  TEST_REGION, TEST_SCAFFOLD_REGION, TEST_SEED (default seed 15052011)
EOF
}

UPDATE=0
if [[ ${1:-} == "--help" || ${1:-} == "-h" ]]; then
  usage; exit 0
fi
if [[ ${1:-} == "--list" ]]; then
  LIST_ONLY=1
else
  LIST_ONLY=0
fi
if [[ ${1:-} == "--update" ]]; then
  UPDATE=1
fi

# Minimal curated suite (fast, deterministic)
SUITE=(
  "${SCRIPT_DIR}/phase.wgs.unrelated.sh"
  "${SCRIPT_DIR}/phase.oneallele.wgs.unrelated.sh"
  "${SCRIPT_DIR}/phase.oneallele.wgs.unrelated.micro.sh"
  "${SCRIPT_DIR}/phase.oneallele.wgs.unrelated.micro_donor.sh"
)

if [[ $LIST_ONLY -eq 1 ]]; then
  printf '%s\n' "${SUITE[@]}"
  exit 0
fi

export TEST_SEED="${TEST_SEED:-15052011}"
export UPDATE_EXPECTED=$UPDATE

pass=0
fail=0

for t in "${SUITE[@]}"; do
  echo "--- Running: $(basename "$t") ---"
  if bash "$t"; then
    echo "✓ $(basename "$t") passed"
    pass=$((pass+1))
  else
    echo "✗ $(basename "$t") failed" >&2
    fail=$((fail+1))
    exit 1
  fi
done

echo "Summary: passed=$pass failed=$fail"
exit $(( fail == 0 ? 0 : 1 ))


#!/bin/bash
# Regression harness for the supersite finalization crash.
# Runs phase_common on a chr22 WGS slice with supersites enabled and captures
# stdout/stderr logs so the std::vector assertion can be reproduced locally.

set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

TEST_NAME="supersite_finalization_regression"
THREADS=${THREADS:-8}
REGION=${REGION:-chr22:19000000-20000000}
ITER_SCHEDULE=${ITER_SCHEDULE:-3b,1p,1b,1p,3m}
INPUT_BCF=${INPUT_BCF:-"$TEST_DIR/wgs/1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.unphased.native_maps.biallelic.18000000-25000000.bcf"}
MAP_FILE=${MAP_FILE:-"$TEST_DIR/info/chr22.gmap.gz"}
PHASE_BIN=${PHASE_BIN:-"$TEST_DIR/../phase_common/bin/phase_common"}
LOG_STEM=${LOG_STEM:-$TEST_NAME}
OUT_PREFIX="$TEST_DIR/tmp/$LOG_STEM"
STDOUT_LOG="$OUT_PREFIX.stdout.log"
STDERR_LOG="$OUT_PREFIX.stderr.log"
OUTPUT_BCF="$OUT_PREFIX.phased.bcf"

mkdir -p "$TEST_DIR/tmp"

if [[ ! -x "$PHASE_BIN" ]]; then
  echo "[$TEST_NAME] ERROR: phase_common binary not found at $PHASE_BIN" >&2
  exit 1
fi

if [[ ! -f "$INPUT_BCF" ]]; then
  echo "[$TEST_NAME] ERROR: input BCF not found at $INPUT_BCF" >&2
  exit 1
fi

if [[ ! -f "$MAP_FILE" ]]; then
  echo "[$TEST_NAME] ERROR: genetic map not found at $MAP_FILE" >&2
  exit 1
fi

rm -f "$STDOUT_LOG" "$STDERR_LOG" "$OUTPUT_BCF" "$OUTPUT_BCF.csi"

CMD=(
  "$PHASE_BIN"
  --input "$INPUT_BCF"
  --filter-maf 0.001
  --region "$REGION"
  --map "$MAP_FILE"
  --output "$OUTPUT_BCF"
  --enable-supersites
  --mcmc-iteration "$ITER_SCHEDULE"
  --thread "$THREADS"
)

echo "[$TEST_NAME] Starting phase_common run" | tee "$STDOUT_LOG"
echo "[$TEST_NAME] Command: ${CMD[*]}" | tee -a "$STDOUT_LOG"
echo "[$TEST_NAME] Logs: $STDOUT_LOG / $STDERR_LOG" | tee -a "$STDOUT_LOG"

start_time=$(date +%s.%N)
set +e
"${CMD[@]}" >>"$STDOUT_LOG" 2>>"$STDERR_LOG"
exit_code=$?
set -e
duration=$(echo "$(date +%s.%N) - $start_time" | bc -l 2>/dev/null || echo "0.000")

if [[ $exit_code -ne 0 ]]; then
  reason="phase_common exited with status $exit_code (see $STDERR_LOG)."
  if grep -q "__n < this->size()" "$STDERR_LOG"; then
    reason+=" Detected std::vector out-of-bounds assertion."
  fi
  test_fail "$TEST_NAME" "$duration" "$reason"
  echo "[$TEST_NAME] Tail of stderr:" >&2
  tail -n 25 "$STDERR_LOG" >&2 || true
  exit 1
fi

test_pass "$TEST_NAME" "$duration"
exit 0

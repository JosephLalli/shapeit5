#!/bin/bash
# Remove -e flag to prevent crashes on command failures
set -uo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

# Test configuration
TEST_NAME="phase.array.unrelated"
START_TIME=$(date +%s.%N)

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

output_bcf="$tmp_dir/target.phased.bcf"
region="${TEST_REGION:-1:5000000-6000000}"

# Run phase_common with error handling
log_file="$tmp_dir/phase_common.log"
if ! ../phase_common/bin/phase_common \
  --input array/target.unrelated.bcf \
  --region "$region" \
  --map info/chr1.gmap.gz \
  --output "$output_bcf" \
  --thread 8 >"$log_file" 2>&1; then
  
  duration=$(echo "$(date +%s.%N) - $START_TIME" | bc -l 2>/dev/null || echo "0.000")
  test_fail "$TEST_NAME" "$duration" "phase_common execution failed"
  echo "Phase common log:" >&2
  cat "$log_file" >&2
  exit 1
fi

# Validate output
assert_same_variants "$output_bcf" "$SCRIPT_DIR/expected/phase.array.unrelated.vcf" "$TEST_NAME" "$START_TIME"

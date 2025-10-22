#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

tmp_dir="$TEST_DIR/tmp"
mkdir -p "$tmp_dir"

scaffold_region="${TEST_SCAFFOLD_REGION:-chrX:153929053-154248138}"
comparison_region="${TEST_REGION:-chrX:153929053-154248138}"

FAM=info/1kgp_t2t.par2.ped
scaffold_bcf="$tmp_dir/target.scaffold.wgs.family.bcf"
output_bcf="$tmp_dir/target.phased.wgs.family.bcf"

# Phase common variants first
../phase_common/bin/phase_common \
  --input wgs/target.family.1kgp_t2t.par2.bcf \
  --filter-maf 0.001 \
  --pedigree "$FAM" \
  --region "$scaffold_region" \
  --map info/par2.gmap.gz \
  --seed 15052011 \
  --output "$scaffold_bcf" \
  --thread 1

# Phase rare variants in one go (no chunking)
log_file="$tmp_dir/phase_rare.log"
if ../phase_rare/bin/phase_rare \
    --input wgs/target.family.1kgp_t2t.par2.bcf \
    --scaffold "$scaffold_bcf" \
    --pedigree "$FAM" \
    --map info/par2.gmap.gz \
    --input-region "$scaffold_region" \
    --scaffold-region "$scaffold_region" \
    --output "$output_bcf" \
    --seed 15052011 \
    --thread 1 >"$log_file" 2>&1; then
  echo "Rare variant phasing completed successfully"
else
  if grep -q "No variants to be phased" "$log_file"; then
    echo "No rare variants found; using scaffold output only" >&2
    cp "$scaffold_bcf" "$output_bcf"
    cp "${scaffold_bcf}.csi" "${output_bcf}.csi"
  else
    echo "Rare variant phasing failed:" >&2
    cat "$log_file" >&2
    exit 1
  fi
fi

filtered_bcf="$tmp_dir/target.phased.wgs.family.filtered.bcf"
SSH_AUTH_SOCK= bcftools view -Ob -o "$filtered_bcf" -r "$comparison_region" "$output_bcf"

assert_same_md5 "$filtered_bcf" "phase.wgs.family"

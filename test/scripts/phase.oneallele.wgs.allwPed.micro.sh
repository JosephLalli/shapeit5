#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"


#tmp_dir="$TEST_DIR/tmp"
inbcf=$2
tmp_dir=$1
mkdir -p "$tmp_dir"

pedigree=$TEST_DIR/1kgp.ped
#scaffold_region="${TEST_SCAFFOLD_REGION:-chrX:153929053-154248138}"
#comparison_region="${TEST_REGION:-chrX:153929053-154248138}"
scaffold_region=$3
comparison_region=$scaffold_region

scaffold_bcf="$tmp_dir/target.scaffold.oneallele.micro.wgs.unrelated.bcf"
scaffold_stats="$tmp_dir/target.scaffold.oneallele.micro.wgs.unrelated.stats"
scaffold_debug="$tmp_dir/target.scaffold.oneallele.micro.wgs.unrelated.debug"
output_bcf="$tmp_dir/$(basename $inbcf).phased.bcf"

#e common variants with MICRO mode
../phase_common/bin/phase_common \
  --input $inbcf \
  --filter-maf 0.001 \
  --region "$scaffold_region" \
  --map info/par2.gmap.gz \
  --pedigree $pedigree \
  --enforce-oneallele \
  --oneallele-mode micro-donor \
  --oneallele-debug "$scaffold_debug" \
  --oneallele-stats "$scaffold_stats" \
  --seed 15052011 \
  --output "$scaffold_bcf"

# Phase rare variants (sparse-micro enforcement)
log_file="$tmp_dir/phase_rare.micro.log"
if ../phase_rare/bin/phase_rare \
    --input $inbcf \
    --scaffold "$scaffold_bcf" \
    --map info/par2.gmap.gz \
    --input-region "$scaffold_region" \
    --scaffold-region "$scaffold_region" \
    --pedigree $pedigree \
    --enforce-oneallele-rare \
    --oneallele-rare-mode sparse-micro \
    --output "$output_bcf" \
    --seed 15052011 >"$log_file" 2>&1; then
  echo "Rare variant phasing (MICRO) completed successfully"
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

filtered_bcf="$tmp_dir/target.phased.oneallele.micro.wgs.unrelated.filtered.bcf"
SSH_AUTH_SOCK= bcftools view -Ob -o "$filtered_bcf" -r "$comparison_region" "$output_bcf"

# Validate one-allele constraint and MD5s
assert_no_oneallele_violations "$scaffold_bcf"
assert_no_oneallele_violations "$filtered_bcf"
assert_same_md5 "$scaffold_bcf" "phase.wgs.unrelated.micro.common"
assert_same_md5 "$filtered_bcf" "phase.wgs.unrelated.micro"


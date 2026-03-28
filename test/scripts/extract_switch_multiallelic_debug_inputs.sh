#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

multi_ref_bcf="${TEST_DIR}/../1KGP.CHM13v2.0.chr22.18-25mb.snp_indel.phasing_qual_pass.fully_annotated.bcf"

region="${1:-chr22:19000000-19100000}"
out_dir="${2:-${TEST_DIR}/tmp/switch_multiallelic_debug}"
mkdir -p "$out_dir"

safe_region="${region//:/_}"
safe_region="${safe_region//-/_}"
prefix="${out_dir}/multi_ref.${safe_region}"

multiallelic_unphased="${prefix}.multiallelic.unphased.bcf"
biallelic_phased="${prefix}.biallelic.phased.bcf"
biallelic_unphased="${prefix}.biallelic.unphased.bcf"

bcftools view -r "$region" -Ou "$multi_ref_bcf" \
  | bcftools +setGT -Ou -- -t a -n u \
  | bcftools view -Ob -o "$multiallelic_unphased"
bcftools index -f "$multiallelic_unphased"

bcftools view -r "$region" -Ou "$multi_ref_bcf" \
  | bcftools norm -m -any -Ou \
  | bcftools view -Ob -o "$biallelic_phased"
bcftools index -f "$biallelic_phased"

bcftools view -r "$region" -Ou "$multi_ref_bcf" \
  | bcftools norm -m -any -Ou \
  | bcftools +setGT -Ou -- -t a -n u \
  | bcftools view -Ob -o "$biallelic_unphased"
bcftools index -f "$biallelic_unphased"

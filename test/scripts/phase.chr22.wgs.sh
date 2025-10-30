#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

source "$SCRIPT_DIR/lib/test_utils.sh"

tmp_dir=$TEST_DIR/tmp
mkdir -p tmp

scaffold_region=chr22:18000000-25000000
comparison_region=chr22:19000000-24000000

in_bcf='wgs/1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.unphased.native_maps.biallelic.18000000-25000000.bcf'
scaffold_bcf_prefix="$tmp_dir/chr22.1KGP.18-25mb.phase_common"
if [[ 'x' == 'y' ]]; then
if [[ ! -s $scaffold_bcf_prefix.og.bcf ]]; then
/usr/bin/time ./SHAPEIT5_phase_common_static_v1.1.1 \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.og.bcf \
  --thread 8
fi

/usr/bin/time ../phase_common/bin/phase_common \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.main_algo.bcf \
  --thread 8
fi
/usr/bin/time ../phase_common/bin/phase_common \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.supersites.bcf \
  --enable-supersites

#assert_same_variants "$filtered_bcf" "$SCRIPT_DIR/expected/phase.wgs.unrelated.vcf"

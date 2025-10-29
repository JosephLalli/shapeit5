#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

infile=$TEST_DIR/whole_genome_tmp/1KGP.CHM13v2.0.PAR1.snp_indel.phasing_qual_pass.biallelic.bcf
phased_truth=$TEST_DIR/whole_genome_tmp/1KGP.CHM13v2.0.PAR1.snp_indel.phasing_qual_pass.biallelic.bcf.phased.bcf
source "$SCRIPT_DIR/lib/test_utils.sh"
threads=32

tmp_dir=./tmp
mkdir -p $tmp_dir

scaffold_region="chrX:0-2300000"
comparison_region="chrX:0-2300000"

scaffold_bcf="$tmp_dir/target.scaffold"

if [[ 'false' == true ]] ; then
/usr/bin/time ../phase_common/bin/phase_common \
  --input $infile \
  --filter-maf 0.001 \
  --region "$scaffold_region" \
  --map info/PAR1.gmap.gz \
  --output "${scaffold_bcf}.new.bcf" \
  --thread $threads

/usr/bin/time ../../SHAPEIT5_phase_common_static_v1.1.1 \
  --input $infile \
  --filter-maf 0.001 \
  --region "$scaffold_region" \
  --map info/PAR1.gmap.gz \
  --output "${scaffold_bcf}.og.bcf" \
  --thread $threads

/usr/bin/time ../phase_common/bin/phase_common \
  --input $infile \
  --filter-maf 0.001 \
  --region "$scaffold_region" \
  --epsilon-indel=0.01 \
  --map info/PAR1.gmap.gz \
  --output "${scaffold_bcf}.indel_0.01.bcf" \
  --thread $threads
fi

../switch/bin/switch \
  --validation $phased_truth \
  --estimation ${scaffold_bcf}.new.bcf \
  --region "$scaffold_region" \
  --output ${scaffold_bcf}.new &

../switch/bin/switch \
  --validation $phased_truth \
  --estimation ${scaffold_bcf}.og.bcf \
  --region "$scaffold_region" \
  --output ${scaffold_bcf}.og &

../switch/bin/switch \
  --validation $phased_truth \
  --estimation ${scaffold_bcf}.indel_0.01.bcf \
  --region "$scaffold_region" \
  --output ${scaffold_bcf}.indel_0.01 &

#OUT="$tmp_dir/target.phased.unrelated.bcf"
#log_file="$tmp_dir/target.phased.unrelated.log"
#../phase_rare/bin/phase_rare \
#      --input $infile \
#      --scaffold "$scaffold_bcf" \
#      --map info/PAR1.gmap.gz \
#      --input-region "$scaffold_region" \
#      --scaffold-region "$scaffold_region" \
#      --output "$OUT" \
#      --thread $threads >"$log_file"


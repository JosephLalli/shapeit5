#!/bin/bash
set -euo pipefail

SCRIPT_PATH=$(realpath "${BASH_SOURCE[0]}")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
TEST_DIR=$(realpath "${SCRIPT_DIR}/..")
cd "$TEST_DIR"

seed=$1
run_suffix=$2
run_suffix=".${run_suffix}"

source "$SCRIPT_DIR/lib/test_utils.sh"

tmp_dir=$TEST_DIR/tmp
mkdir -p tmp

#scaffold_region=chr22:18000000-25000000
#comparison_region=chr22:18000000-25000000
scaffold_region=chr22:19000000-20000000
comparison_region=chr22:19000000-20000000

#scaffold_region=chr22:19000000-19100000
#comparison_region=chr22:19000000-19100000

#in_bcf="${TEST_DIR}/wgs/1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.unphased.native_maps.biallelic.18000000-25000000.renormed.bcf"
in_bcf="${TEST_DIR}/../1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.filtered.original_gt.bcf"
#in_bcf="${TEST_DIR}/../1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.filtered.sorted.bcf"
ref_bcf="${TEST_DIR}/wgs/chr22_t2t_reference_pangenome.filtered_variants.18000000-25000000.biallelic.filtered.bcf"
#in_bcf="${TEST_DIR}/wgs/snp_and_multi.unphased.vcf.gz"
#ref_bcf="${TEST_DIR}/wgs/snp_and_multi.phased.vcf.gz"

threads=32

scaffold_bcf_prefix="$tmp_dir/chr22.1KGP.18-25mb.phase_common"
# if [[ 'x' == 'y' ]]; then
# if [[ ! -s $scaffold_bcf_prefix.og.bcf ]]; then
/usr/bin/time ./SHAPEIT5_phase_common_static_v1.1.1 \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.og.small.bcf \
  --thread $threads &
#fi

/usr/bin/time ${TEST_DIR}/../phase_common/bin/phase_common \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.main_algo.small.${seed}${run_suffix}.bcf \
  --seed $seed \
  --thread $threads &
#fi

${TEST_DIR}/../phase_common/bin/phase_common \
  --input $in_bcf \
  --filter-maf 0.001 \
  --region $scaffold_region \
  --map info/chr22.gmap.gz \
  --output $scaffold_bcf_prefix.supersites.split_emissions.small.${seed}${run_suffix}.bcf \
  --enable-supersites \
  --seed $seed \
  --thread $threads &

#fi
wait

#if [[ 'a' == 'b' ]]; then
${TEST_DIR}/../switch/bin/switch \
  --validation "$ref_bcf" \
  --estimation "$scaffold_bcf_prefix.og.small.bcf" \
  --region $comparison_region \
  --output "$scaffold_bcf_prefix.og" \
  --log "$scaffold_bcf_prefix.og.log" >/dev/null
#fi

../switch/bin/switch \
  --validation "$ref_bcf" \
  --estimation "$scaffold_bcf_prefix.main_algo.small.${seed}${run_suffix}.bcf" \
  --region $comparison_region \
  --output "$scaffold_bcf_prefix.main_algo.${seed}${run_suffix}" \
  --log "$scaffold_bcf_prefix.main_algo.${seed}${run_suffix}.log" >/dev/null

../switch/bin/switch \
  --validation "$ref_bcf" \
  --estimation "$scaffold_bcf_prefix.supersites.split_emissions.small.${seed}${run_suffix}.bcf" \
  --region $comparison_region \
  --output "$scaffold_bcf_prefix.supersites.split_emissions.${seed}${run_suffix}" \
  --log "$scaffold_bcf_prefix.supersites.split_emissions.${seed}${run_suffix}.log" >/dev/null

# Integration testing: Extract and compare switch error rates
echo "=== Integration Test Results ==="

# Extract switch error rates from log files
extract_switch_error() {
    local log_file="$1"
    if [[ -f "$log_file" ]]; then
        grep "Overall switch error rate" "$log_file" | cut -f 9 -d ' '
    else
        echo "N/A"
    fi
}

og_error=$(extract_switch_error "$scaffold_bcf_prefix.og.log")
main_algo_error=$(extract_switch_error "$scaffold_bcf_prefix.main_algo.${seed}${run_suffix}.log")
supersites_error=$(extract_switch_error "$scaffold_bcf_prefix.supersites.split_emissions.${seed}${run_suffix}.log")

echo "Original algorithm switch error rate: $og_error"
echo "Main algorithm switch error rate: $main_algo_error"
echo "Supersites algorithm switch error rate: $supersites_error"

# Test logic: Fail if main_algo is worse than og
if [[ "$og_error" != "N/A" && "$main_algo_error" != "N/A" ]]; then
    if (( $(echo "$main_algo_error > $og_error" | bc -l) )); then
        echo "ERROR: Main algorithm switch error rate ($main_algo_error) is higher than original ($og_error)"
        exit 1
    fi
fi

# Compare supersites to main_algo
if [[ "$main_algo_error" != "N/A" && "$supersites_error" != "N/A" ]]; then
    improvement=$(echo "scale=6; $main_algo_error - $supersites_error" | bc -l)
    if (( $(echo "$improvement > 0" | bc -l) )); then
        percent_improvement=$(echo "scale=2; ($improvement / $main_algo_error) * 100" | bc -l)
        echo "SUCCESS: Supersites is ${percent_improvement}% better than main algorithm (improvement: $improvement)"
    elif (( $(echo "$improvement < 0" | bc -l) )); then
        percent_worse=$(echo "scale=2; (${improvement#-} / $main_algo_error) * 100" | bc -l)
        echo "WARNING: Supersites is ${percent_worse}% worse than main algorithm (degradation: ${improvement#-})"
    else
        echo "INFO: Supersites and main algorithm have identical switch error rates"
    fi
fi

echo "=== Test completed successfully ==="

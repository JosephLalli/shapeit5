#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 chr22:19300000-19400000" >&2
  exit 1
fi

REGION="$1"
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/tmp_multi/bench_leftshift}"
MAP_FILE="${MAP_FILE:-$ROOT_DIR/test/info/chr22.gmap.gz}"
REF_FASTA="${REF_FASTA:-/mnt/ssd/lalli/phasing_T2T/chm13v2.0.fa}"
THREADS_PHASE="${THREADS_PHASE:-16}"
THREADS_SWITCH="${THREADS_SWITCH:-8}"
SEED="${SEED:-42}"

export LD_LIBRARY_PATH="$HOME/.linuxbrew/lib/perl5/5.42/x86_64-linux-thread-multi/CORE:$HOME/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH:-}"

PHASE_COMMON="$ROOT_DIR/phase_common/bin/phase_common"
SWITCH_BIN="$ROOT_DIR/switch/bin/switch"

require_file() {
  if [[ ! -s "$1" ]]; then
    echo "missing required file: $1" >&2
    exit 1
  fi
}

require_exec() {
  if [[ ! -x "$1" ]]; then
    echo "missing executable: $1" >&2
    exit 1
  fi
}

require_exec "$PHASE_COMMON"
require_exec "$SWITCH_BIN"
command -v bcftools >/dev/null 2>&1 || { echo "missing bcftools in PATH" >&2; exit 1; }

MULTI_UNPHASED="$ROOT_DIR/tmp_multi/1KGP.CHM13v2.0.chr22.refiltered_from_scratch.multi.unphased.bcf"
BIAL_UNPHASED="$ROOT_DIR/tmp_multi/1KGP.CHM13v2.0.chr22.refiltered_from_scratch.biallelic.unphased.bcf"
BIAL_LEFTSHIFT_UNPHASED="$ROOT_DIR/1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.bcf"

TRUTH_MULTI="$ROOT_DIR/tmp_multi/hprc-v2.0-mc-chm13.multi.filtered.phased.bcf"
TRUTH_BIAL="$ROOT_DIR/tmp_multi/hprc-v2.0-mc-chm13.multi.filtered.biallelic.phased.bcf"
TRUTH_BIAL_LEFTSHIFT="$ROOT_DIR/tmp_multi/hprc-v2.0-mc-chm13.multi.filtered.biallelic.phased.leftshifted.bcf"

require_file "$MULTI_UNPHASED"
require_file "$BIAL_UNPHASED"
require_file "$BIAL_LEFTSHIFT_UNPHASED"
require_file "$TRUTH_MULTI"
require_file "$TRUTH_BIAL"
require_file "$TRUTH_BIAL_LEFTSHIFT"
require_file "$MAP_FILE"
require_file "$REF_FASTA"

mkdir -p "$OUT_DIR"

PHASE_MULTI="$OUT_DIR/phase.1kgp.multi.noshift.phased.bcf"
PHASE_BIAL="$OUT_DIR/phase.1kgp.bial.noshift.phased.bcf"
PHASE_BIAL_LEFTSHIFT="$OUT_DIR/phase.1kgp.bial.leftshift.phased.bcf"
SPLIT_MULTI="$OUT_DIR/phase.1kgp.multi.noshift.phased.split.biallelic.bcf"
SPLIT_MULTI_LEFTSHIFT="$OUT_DIR/phase.1kgp.multi.noshift.phased.split.biallelic.leftshifted.bcf"

SW_MULTI_PREFIX="$OUT_DIR/eval.1kgp.multi.noshift.switch"
SW_SPLIT_PREFIX="$OUT_DIR/eval.1kgp.multi.noshift.split_bial.switch"
SW_SPLIT_LEFTSHIFT_PREFIX="$OUT_DIR/eval.1kgp.multi.noshift.split_bial.leftshift.switch"
SW_BIAL_PREFIX="$OUT_DIR/eval.1kgp.bial.noshift.switch"
SW_BIAL_LEFTSHIFT_PREFIX="$OUT_DIR/eval.1kgp.bial.leftshift.switch"
if [[ 'x' == 'y' ]]; then
"$PHASE_COMMON" \
  --input "$MULTI_UNPHASED" \
  --map "$MAP_FILE" \
  --region "$REGION" \
  --output "$PHASE_MULTI" \
  --output-format bcf \
  --thread "$THREADS_PHASE" \
  --seed "$SEED" \
  --enable-supersites \
  --log "${PHASE_MULTI%.bcf}.log" && \
bcftools index -f "$PHASE_MULTI" &

"$PHASE_COMMON" \
  --input "$BIAL_UNPHASED" \
  --map "$MAP_FILE" \
  --region "$REGION" \
  --output "$PHASE_BIAL" \
  --output-format bcf \
  --thread "$THREADS_PHASE" \
  --seed "$SEED" \
  --log "${PHASE_BIAL%.bcf}.log" && \
bcftools index -f "$PHASE_BIAL" &

"$PHASE_COMMON" \
  --input "$BIAL_LEFTSHIFT_UNPHASED" \
  --map "$MAP_FILE" \
  --region "$REGION" \
  --output "$PHASE_BIAL_LEFTSHIFT" \
  --output-format bcf \
  --thread "$THREADS_PHASE" \
  --seed "$SEED" \
  --log "${PHASE_BIAL_LEFTSHIFT%.bcf}.log" && \
bcftools index -f "$PHASE_BIAL_LEFTSHIFT" &

wait

bcftools norm -m -any -O b -W -o "$SPLIT_MULTI" "$PHASE_MULTI" && bcftools norm -f "$REF_FASTA" -O b -W -o "$SPLIT_MULTI_LEFTSHIFT" "$SPLIT_MULTI"

wait

"$SWITCH_BIN" \
  --validation "$TRUTH_MULTI" \
  --estimation "$PHASE_MULTI" \
  --region "$REGION" \
  --output "$SW_MULTI_PREFIX" \
  --log "${SW_MULTI_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIAL" \
  --estimation "$SPLIT_MULTI" \
  --region "$REGION" \
  --output "$SW_SPLIT_PREFIX" \
  --log "${SW_SPLIT_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIAL_LEFTSHIFT" \
  --estimation "$SPLIT_MULTI_LEFTSHIFT" \
  --region "$REGION" \
  --output "$SW_SPLIT_LEFTSHIFT_PREFIX" \
  --log "${SW_SPLIT_LEFTSHIFT_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIAL" \
  --estimation "$PHASE_BIAL" \
  --region "$REGION" \
  --output "$SW_BIAL_PREFIX" \
  --log "${SW_BIAL_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIAL_LEFTSHIFT" \
  --estimation "$PHASE_BIAL_LEFTSHIFT" \
  --region "$REGION" \
  --output "$SW_BIAL_LEFTSHIFT_PREFIX" \
  --log "${SW_BIAL_LEFTSHIFT_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

wait
fi

python3 "$ROOT_DIR/scripts/bench_leftshift_report.py" --out-dir "$OUT_DIR"

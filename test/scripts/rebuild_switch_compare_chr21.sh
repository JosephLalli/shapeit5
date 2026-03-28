#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)

WORK_DIR="${1:-${ROOT_DIR}/test/tmp/hprc_compare_chr21_11300000_11400000}"
SRC_VCF="${2:-${ROOT_DIR}/../hprc-v2.0-mc-chm13.pgin.vcf.gz}"
REGION="${3:-chr21:11300000-11400000}"

HPRC_SAMPLES="${HPRC_SAMPLES:-${WORK_DIR}/samples.hprc.txt}"
PANEL_SAMPLES="${PANEL_SAMPLES:-${WORK_DIR}/samples.panel.txt}"

BCFTOOLS_LD="/mnt/ssd/lalli/.linuxbrew/lib/perl5/5.42/x86_64-linux-thread-multi/CORE:${LD_LIBRARY_PATH:-}"
SWITCH_LD="/mnt/ssd/lalli/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH:-}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
SWITCH_BIN="${SWITCH_BIN:-${ROOT_DIR}/switch/bin/switch}"
COMPARE_PY="${COMPARE_PY:-${ROOT_DIR}/test/scripts/lib/compare_switch_bial_multiallelic.py}"

overlap_file="${WORK_DIR}/samples.overlap.txt"

hprc_bial_out="${WORK_DIR}/hprc.chr21_11300000_11400000.biallelic.le50.gt_ac_an_maf.bcf"
hprc_multi_out="${WORK_DIR}/hprc.chr21_11300000_11400000.multiallelic.le50.gt_ac_an_maf.bcf"
panel_bial_out="${WORK_DIR}/panel.chr21_11300000_11400000.biallelic.le50.gt_ac_an_maf.bcf"
panel_multi_out="${WORK_DIR}/panel.chr21_11300000_11400000.multiallelic.le50.gt_ac_an_maf.bcf"

bcftools_cmd() {
  LD_LIBRARY_PATH="$BCFTOOLS_LD" "$BCFTOOLS_BIN" "$@"
}

switch_cmd() {
  LD_LIBRARY_PATH="$SWITCH_LD" "$SWITCH_BIN" "$@"
}

if [[ ! -f "$SRC_VCF" ]]; then
  echo "ERROR: missing source VCF: $SRC_VCF" >&2
  exit 1
fi

if [[ ! -f "$SRC_VCF.tbi" && ! -f "$SRC_VCF.csi" ]]; then
  echo "ERROR: missing index for source VCF: ${SRC_VCF}.tbi/.csi" >&2
  exit 1
fi

if [[ ! -f "$HPRC_SAMPLES" ]]; then
  echo "ERROR: missing HPRC sample list: $HPRC_SAMPLES" >&2
  exit 1
fi

if [[ ! -f "$PANEL_SAMPLES" ]]; then
  echo "ERROR: missing panel sample list: $PANEL_SAMPLES" >&2
  exit 1
fi

if [[ ! -f "$overlap_file" ]]; then
  echo "ERROR: missing overlap file: $overlap_file" >&2
  exit 1
fi

build_raw() {
  local samples=$1
  local out_bial=$2
  local out_multi=$3

  bcftools_cmd view -r "$REGION" -S "$samples" --force-samples -Ou "$SRC_VCF" \
    | bcftools_cmd norm -m -any -Ou \
    | bcftools_cmd view -i 'strlen(REF)<=50 && strlen(ALT)<=50' -Ob -o "$out_bial"

  bcftools_cmd view -r "$REGION" -S "$samples" --force-samples -Ou "$SRC_VCF" \
    | bcftools_cmd norm -m -any -Ou \
    | bcftools_cmd view -i 'strlen(REF)<=50 && strlen(ALT)<=50' -Ou \
    | bcftools_cmd norm -m +any -Ob -o "$out_multi"
}

clean_bcf() {
  local in_bcf=$1
  local out_bcf=$2
  local tag_label=$3
  local tmp_fill="${WORK_DIR}/tmp.${tag_label}.fill.bcf"
  local tmp_keep="${WORK_DIR}/tmp.${tag_label}.keep.bcf"

  bcftools_cmd +fill-tags "$in_bcf" -Ob -o "$tmp_fill" -- -t AC,AN,MAF
  bcftools_cmd annotate -x '^INFO/AC,^INFO/AN,^INFO/MAF,^FORMAT/GT' "$tmp_fill" -Ob -o "$tmp_keep"
  mv -f "$tmp_keep" "$out_bcf"
  bcftools_cmd index -f "$out_bcf"
  rm -f "$tmp_fill"

  headfile="${WORK_DIR}/tmp.${tag_label}.head.txt"
  bcftools_cmd view -H "$out_bcf" -o "$headfile"
  echo "### ${tag_label}"
  head -n 1 "$headfile" | cut -f 1-12
  fmt=$(head -n 1 "$headfile" | cut -f 9)
  rm -f "$headfile"

  if [[ "$fmt" != *GT* ]]; then
    echo "ERROR: FORMAT column missing GT for $out_bcf" >&2
    exit 1
  fi
}

tmp_hprc_bial="${WORK_DIR}/tmp.hprc.raw.biallelic.bcf"
tmp_hprc_multi="${WORK_DIR}/tmp.hprc.raw.multiallelic.bcf"
tmp_panel_bial="${WORK_DIR}/tmp.panel.raw.biallelic.bcf"
tmp_panel_multi="${WORK_DIR}/tmp.panel.raw.multiallelic.bcf"

build_raw "$HPRC_SAMPLES" "$tmp_hprc_bial" "$tmp_hprc_multi"
build_raw "$PANEL_SAMPLES" "$tmp_panel_bial" "$tmp_panel_multi"

clean_bcf "$tmp_hprc_bial" "$hprc_bial_out" "hprc.biallelic"
clean_bcf "$tmp_hprc_multi" "$hprc_multi_out" "hprc.multiallelic"
clean_bcf "$tmp_panel_bial" "$panel_bial_out" "panel.biallelic"
clean_bcf "$tmp_panel_multi" "$panel_multi_out" "panel.multiallelic"

rm -f "$tmp_hprc_bial" "$tmp_hprc_multi" "$tmp_panel_bial" "$tmp_panel_multi"

val_bial="$hprc_bial_out"
est_bial="$panel_bial_out"
val_multi="$hprc_multi_out"
est_multi="$panel_multi_out"

tmp_val_bial="${WORK_DIR}/tmp.switch.validation.biallelic.bcf"
tmp_est_bial="${WORK_DIR}/tmp.switch.estimation.biallelic.bcf"
tmp_val_multi="${WORK_DIR}/tmp.switch.validation.multiallelic.bcf"
tmp_est_multi="${WORK_DIR}/tmp.switch.estimation.multiallelic.bcf"

bcftools_cmd view -S "$overlap_file" -Ob -o "$tmp_val_bial" "$val_bial"
bcftools_cmd view -S "$overlap_file" -Ob -o "$tmp_est_bial" "$est_bial"
bcftools_cmd view -S "$overlap_file" -Ob -o "$tmp_val_multi" "$val_multi"
bcftools_cmd view -S "$overlap_file" -Ob -o "$tmp_est_multi" "$est_multi"
bcftools_cmd index -f "$tmp_val_bial"
bcftools_cmd index -f "$tmp_est_bial"
bcftools_cmd index -f "$tmp_val_multi"
bcftools_cmd index -f "$tmp_est_multi"

bial_prefix="${WORK_DIR}/switch.panel_vs_hprc.biallelic"
multi_prefix="${WORK_DIR}/switch.panel_vs_hprc.multiallelic"

switch_cmd \
  --validation "$tmp_val_bial" \
  --estimation "$tmp_est_bial" \
  --region "$REGION" \
  --output "$bial_prefix" \
  --log "${bial_prefix}.log"

switch_cmd \
  --validation "$tmp_val_multi" \
  --estimation "$tmp_est_multi" \
  --region "$REGION" \
  --output "$multi_prefix" \
  --log "${multi_prefix}.log"

compare_out="$(python3 "$COMPARE_PY" \
  --multi-variant "${multi_prefix}.variant.switch.txt.gz" \
  --bial-variant "${bial_prefix}.variant.switch.txt.gz" \
  --multi-sample "${multi_prefix}.sample.switch.txt.gz" \
  --bial-sample "${bial_prefix}.sample.switch.txt.gz" \
  --out-prefix "${WORK_DIR}/switch.compare")"

printf '%s\n' "$compare_out"

if echo "$compare_out" | rg -q "positions with differing \\(errors,checked\\): 0" \
  && echo "$compare_out" | rg -q "samples with differing \\(errors,checked\\): 0"; then
  exit 0
fi

bial_site_log="${bial_prefix}.site.log"
multi_site_log="${multi_prefix}.site.log"

switch_cmd \
  --validation "$tmp_val_bial" \
  --estimation "$tmp_est_bial" \
  --region "$REGION" \
  --output "$bial_prefix" \
  --log "${bial_prefix}.site.log.run.log" \
  --site-log "$bial_site_log"

switch_cmd \
  --validation "$tmp_val_multi" \
  --estimation "$tmp_est_multi" \
  --region "$REGION" \
  --output "$multi_prefix" \
  --log "${multi_prefix}.site.log.run.log" \
  --site-log "$multi_site_log"

python3 "$COMPARE_PY" \
  --multi-variant "${multi_prefix}.variant.switch.txt.gz" \
  --bial-variant "${bial_prefix}.variant.switch.txt.gz" \
  --multi-sample "${multi_prefix}.sample.switch.txt.gz" \
  --bial-sample "${bial_prefix}.sample.switch.txt.gz" \
  --multi-site-log "$multi_site_log" \
  --bial-site-log "$bial_site_log" \
  --out-prefix "${WORK_DIR}/switch.compare.site"

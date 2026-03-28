#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT_DIR=$(cd "${SCRIPT_DIR}/../.." && pwd)

REGION="${2:-chr21:19000000-19100000}"
REGION_TAG="${REGION//:/_}"
REGION_TAG="${REGION_TAG//-/_}"

WORK_DIR="${1:-${ROOT_DIR}/test/tmp/hprc_compare_${REGION_TAG}}"
MAP="${3:-${ROOT_DIR}/resources/maps/b38/chr21.b38.gmap.gz}"
THREADS="${THREADS:-1}"

BCFTOOLS_LD="/mnt/ssd/lalli/.linuxbrew/lib/perl5/5.42/x86_64-linux-thread-multi/CORE:${LD_LIBRARY_PATH:-}"
PHASE_LD="/mnt/ssd/lalli/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH:-}"
SWITCH_LD="/mnt/ssd/lalli/.linuxbrew/lib:/usr/local/lib:${LD_LIBRARY_PATH:-}"

BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
PHASE_BIN="${PHASE_BIN:-${ROOT_DIR}/phase_common/bin/phase_common}"
SWITCH_BIN="${SWITCH_BIN:-${ROOT_DIR}/switch/bin/switch}"

NOT_PARENTS="${NOT_PARENTS:-${HOME}/phasing_T2T/sample_subsets/not_parents.txt}"
HPRC_SRC_VCF="${HPRC_SRC_VCF:-${ROOT_DIR}/../hprc-v2.0-mc-chm13.pgin.vcf.gz}"
PANEL_BIAL_SRC="${PANEL_BIAL_SRC:-${ROOT_DIR}/../phased_T2T_panel_T2T_scaled_newimpute_092424/1KGP.CHM13v2.0.chr21.recalibrated.snp_indel.pass.phased.native_maps.biallelic.2504.bcf}"
PANEL_MULTI_SRC="${PANEL_MULTI_SRC:-${ROOT_DIR}/../phased_T2T_panel_T2T_scaled_newimpute_092424/1KGP.CHM13v2.0.chr21.recalibrated.snp_indel.pass.phased.native_maps.2504.vcf.gz}"

DEFAULT_SAMPLE_DIR="${ROOT_DIR}/test/tmp/hprc_compare_chr21_11300000_11400000"
HPRC_SAMPLES="${HPRC_SAMPLES:-${WORK_DIR}/samples.hprc.txt}"

HPRC_PREFIX="${WORK_DIR}/hprc.${REGION_TAG}"
PANEL_PREFIX="${WORK_DIR}/panel.${REGION_TAG}"
HPRC_BIAL="${HPRC_PREFIX}.biallelic.le50.gt_ac_an_maf.bcf"
HPRC_MULTI="${HPRC_PREFIX}.multiallelic.le50.gt_ac_an_maf.bcf"

bcftools_cmd() {
  LD_LIBRARY_PATH="$BCFTOOLS_LD" "$BCFTOOLS_BIN" "$@"
}

phase_cmd() {
  LD_LIBRARY_PATH="$PHASE_LD" "$PHASE_BIN" "$@"
}

switch_cmd() {
  LD_LIBRARY_PATH="$SWITCH_LD" "$SWITCH_BIN" "$@"
}

mkdir -p "$WORK_DIR"

if [[ ! -f "$HPRC_SAMPLES" && -f "${DEFAULT_SAMPLE_DIR}/samples.hprc.txt" ]]; then
  HPRC_SAMPLES="${DEFAULT_SAMPLE_DIR}/samples.hprc.txt"
fi

if [[ ! -f "${WORK_DIR}/samples.hprc.txt" ]]; then
  cp "$HPRC_SAMPLES" "${WORK_DIR}/samples.hprc.txt"
  HPRC_SAMPLES="${WORK_DIR}/samples.hprc.txt"
fi

for path in "$MAP" "$NOT_PARENTS" "$HPRC_SAMPLES" "$HPRC_SRC_VCF" "$PANEL_BIAL_SRC" "$PANEL_MULTI_SRC"; do
  if [[ ! -f "$path" ]]; then
    echo "ERROR: missing file: $path" >&2
    exit 1
  fi
done

if [[ ! -f "${HPRC_SRC_VCF}.tbi" && ! -f "${HPRC_SRC_VCF}.csi" ]]; then
  echo "ERROR: missing index for HPRC VCF: ${HPRC_SRC_VCF}.tbi/.csi" >&2
  exit 1
fi

if [[ ! -f "${PANEL_BIAL_SRC}.tbi" && ! -f "${PANEL_BIAL_SRC}.csi" ]]; then
  echo "ERROR: missing index for panel biallelic: ${PANEL_BIAL_SRC}.tbi/.csi" >&2
  exit 1
fi

if [[ ! -f "${PANEL_MULTI_SRC}.tbi" && ! -f "${PANEL_MULTI_SRC}.csi" ]]; then
  echo "ERROR: missing index for panel multiallelic: ${PANEL_MULTI_SRC}.tbi/.csi" >&2
  exit 1
fi

panel_subset_samples="${WORK_DIR}/samples.panel.not_parents.no_hprc.txt"
tmp_src_samples="${WORK_DIR}/tmp.src.samples.txt"
tmp_not_parents="${WORK_DIR}/tmp.not_parents.no_hprc.txt"

bcftools_cmd query -l "$PANEL_BIAL_SRC" | sort -u > "$tmp_src_samples"
awk 'NR==FNR{h[$1]=1;next}!h[$1]' "$HPRC_SAMPLES" "$NOT_PARENTS" | sort -u > "$tmp_not_parents"
comm -12 "$tmp_src_samples" "$tmp_not_parents" > "$panel_subset_samples"

if [[ ! -s "$panel_subset_samples" ]]; then
  echo "ERROR: panel subset sample list is empty: $panel_subset_samples" >&2
  exit 1
fi

clean_bial_from_src() {
  local samples=$1
  local src=$2
  local out_bial=$3
  local label=$4

  local tmp_bial="${WORK_DIR}/tmp.${label}.biallelic.bcf"
  local tmp_multi="${WORK_DIR}/tmp.${label}.multiallelic.bcf"
  local tmp_fill="${WORK_DIR}/tmp.${label}.fill.bcf"
  local tmp_keep="${WORK_DIR}/tmp.${label}.keep.bcf"

  bcftools_cmd view -r "$REGION" -S "$samples" --force-samples -Ou "$src" \
    | bcftools_cmd norm -m -any -Ou \
    | bcftools_cmd view -i 'strlen(REF)<=50 && strlen(ALT)<=50' -Ob -o "$tmp_bial"

  bcftools_cmd +fill-tags "$tmp_bial" -Ob -o "$tmp_fill" -- -t AC,AN,MAF
  bcftools_cmd annotate -x '^INFO/AC,^INFO/AN,^INFO/MAF,^FORMAT/GT' "$tmp_fill" -Ob -o "$tmp_keep"
  mv -f "$tmp_keep" "$out_bial"
  bcftools_cmd index -f "$out_bial"
  rm -f "$tmp_fill" "$tmp_bial"
}

clean_multi_from_src() {
  local samples=$1
  local src=$2
  local out_multi=$3
  local label=$4

  local tmp_multi="${WORK_DIR}/tmp.${label}.multiallelic.bcf"
  local tmp_fill="${WORK_DIR}/tmp.${label}.fill.bcf"
  local tmp_keep="${WORK_DIR}/tmp.${label}.keep.bcf"

  bcftools_cmd view -r "$REGION" -S "$samples" --force-samples -Ou "$src" \
    | bcftools_cmd norm -m -any -Ou \
    | bcftools_cmd view -i 'strlen(REF)<=50 && strlen(ALT)<=50' -Ou \
    | bcftools_cmd norm -m +any -Ob -o "$tmp_multi"

  bcftools_cmd +fill-tags "$tmp_multi" -Ob -o "$tmp_fill" -- -t AC,AN,MAF
  bcftools_cmd annotate -x '^INFO/AC,^INFO/AN,^INFO/MAF,^FORMAT/GT' "$tmp_fill" -Ob -o "$tmp_keep"
  mv -f "$tmp_keep" "$out_multi"
  bcftools_cmd index -f "$out_multi"
  rm -f "$tmp_fill" "$tmp_multi"
}

panel_bial_subset="${PANEL_PREFIX}.biallelic.le50.gt_ac_an_maf.not_parents.no_hprc.bcf"
panel_multi_subset="${PANEL_PREFIX}.multiallelic.le50.gt_ac_an_maf.not_parents.no_hprc.bcf"

clean_bial_from_src "$panel_subset_samples" "$PANEL_BIAL_SRC" "$panel_bial_subset" "panel"
clean_multi_from_src "$panel_subset_samples" "$PANEL_MULTI_SRC" "$panel_multi_subset" "panel"
clean_bial_from_src "$HPRC_SAMPLES" "$HPRC_SRC_VCF" "$HPRC_BIAL" "hprc"
clean_multi_from_src "$HPRC_SAMPLES" "$HPRC_SRC_VCF" "$HPRC_MULTI" "hprc"

hprc_bial_unphased="${HPRC_PREFIX}.biallelic.le50.gt_ac_an_maf.unphased.bcf"
hprc_multi_unphased="${HPRC_PREFIX}.multiallelic.le50.gt_ac_an_maf.unphased.bcf"

bcftools_cmd +setGT "$HPRC_BIAL" -Ob -o "$hprc_bial_unphased" -- -t a -n u
bcftools_cmd index -f "$hprc_bial_unphased"
bcftools_cmd +setGT "$HPRC_MULTI" -Ob -o "$hprc_multi_unphased" -- -t a -n u
bcftools_cmd index -f "$hprc_multi_unphased"

hprc_bial_rephased="${HPRC_PREFIX}.biallelic.rephased.bcf"
hprc_bial_rephased_supersites="${HPRC_PREFIX}.biallelic.rephased.supersites.bcf"
hprc_multi_rephased_supersites="${HPRC_PREFIX}.multiallelic.rephased.supersites.bcf"

phase_cmd \
  --input "$hprc_bial_unphased" \
  --reference "$panel_bial_subset" \
  --map "$MAP" \
  --region "$REGION" \
  --output "$hprc_bial_rephased" \
  --log "${hprc_bial_rephased}.log" \
  --thread "$THREADS"

phase_cmd \
  --input "$hprc_bial_unphased" \
  --reference "$panel_bial_subset" \
  --map "$MAP" \
  --region "$REGION" \
  --output "$hprc_bial_rephased_supersites" \
  --log "${hprc_bial_rephased_supersites}.log" \
  --thread "$THREADS" \
  --enable-supersites

phase_cmd \
  --input "$hprc_multi_unphased" \
  --reference "$panel_multi_subset" \
  --map "$MAP" \
  --region "$REGION" \
  --output "$hprc_multi_rephased_supersites" \
  --log "${hprc_multi_rephased_supersites}.log" \
  --thread "$THREADS" \
  --enable-supersites

hprc_multi_rephased_split="${HPRC_PREFIX}.multiallelic.rephased.supersites.split.bcf"
bcftools_cmd norm -m -any -Ob -o "$hprc_multi_rephased_split" "$hprc_multi_rephased_supersites"
bcftools_cmd index -f "$hprc_multi_rephased_split"

switch_cmd \
  --validation "$HPRC_BIAL" \
  --estimation "$hprc_bial_rephased" \
  --region "$REGION" \
  --output "${WORK_DIR}/switch.${REGION_TAG}.rephased.biallelic" \
  --log "${WORK_DIR}/switch.${REGION_TAG}.rephased.biallelic.log"

switch_cmd \
  --validation "$HPRC_BIAL" \
  --estimation "$hprc_bial_rephased_supersites" \
  --region "$REGION" \
  --output "${WORK_DIR}/switch.${REGION_TAG}.rephased.biallelic.supersites" \
  --log "${WORK_DIR}/switch.${REGION_TAG}.rephased.biallelic.supersites.log"

switch_cmd \
  --validation "$HPRC_MULTI" \
  --estimation "$hprc_multi_rephased_supersites" \
  --region "$REGION" \
  --output "${WORK_DIR}/switch.${REGION_TAG}.rephased.multiallelic.supersites" \
  --log "${WORK_DIR}/switch.${REGION_TAG}.rephased.multiallelic.supersites.log"

switch_cmd \
  --validation "$HPRC_BIAL" \
  --estimation "$hprc_multi_rephased_split" \
  --region "$REGION" \
  --output "${WORK_DIR}/switch.${REGION_TAG}.rephased.multiallelic.supersites.split_vs_bial" \
  --log "${WORK_DIR}/switch.${REGION_TAG}.rephased.multiallelic.supersites.split_vs_bial.log"

rm -f "$tmp_src_samples" "$tmp_not_parents"

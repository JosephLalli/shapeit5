#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <validation.vcf.gz> <estimation.vcf.gz> <frequency.vcf.gz> <region> <out_prefix>"
  echo "Optional env: SWITCH_BIN (default: switch/bin/switch)"
  exit 1
fi

VALIDATION="$1"
ESTIMATION="$2"
FREQUENCY="$3"
REGION="$4"
OUT_PREFIX="$5"

SWITCH_BIN="${SWITCH_BIN:-switch/bin/switch}"
SITE_LOG="${OUT_PREFIX}.variant.switch.txt.gz"
BINNED_OUT="${OUT_PREFIX}.variant.switch.binned.csv"

if [[ ! -x "${SWITCH_BIN}" ]]; then
  echo "Missing switch binary at ${SWITCH_BIN}"
  exit 1
fi

"${SWITCH_BIN}" \
  --validation "${VALIDATION}" \
  --estimation "${ESTIMATION}" \
  --frequency "${FREQUENCY}" \
  --region "${REGION}" \
  --output "${OUT_PREFIX}"

python3 tests/scripts/switch_multialt_bins.py \
  --switch-variant "${SITE_LOG}" \
  --frequency-vcf "${FREQUENCY}" \
  --out "${BINNED_OUT}"

echo "Wrote ${BINNED_OUT}"

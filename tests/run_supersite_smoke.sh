#!/usr/bin/env bash
set -euo pipefail
DIR=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$DIR/.." && pwd)

BIN="$ROOT/phase_common/bin/phase_common"
DATA="$DIR/data"
OUTDIR="$DIR/out"
mkdir -p "$OUTDIR"

# Ensure dynamic libs are found for boost/hts from linuxbrew default
export LD_LIBRARY_PATH="${HOME}/.linuxbrew/lib:${LD_LIBRARY_PATH:-}"

log() { echo "[test] $*"; }

# Ensure VCF is bgzip-compressed as required by htslib
if [ ! -f "$DATA/supersite.vcf.gz" ]; then
  log "bgzip compressing supersite.vcf"
  bgzip -c "$DATA/supersite.vcf" > "$DATA/supersite.vcf.gz"
  tabix -p vcf "$DATA/supersite.vcf.gz" || true
fi

if [ ! -f "$DATA/supersite_ref.vcf.gz" ]; then
  log "bgzip compressing supersite_ref.vcf"
  bgzip -c "$DATA/supersite_ref.vcf" > "$DATA/supersite_ref.vcf.gz"
  tabix -p vcf "$DATA/supersite_ref.vcf.gz" || true
fi

log "Running supersite smoke test with --enable-supersites"
"$BIN" \
  --input "$DATA/supersite.vcf.gz" \
  --reference "$DATA/supersite_ref.vcf.gz" \
  --region 1:1-200 \
  --output "$OUTDIR/supersite_enable.graph" \
  --output-format graph \
  --mcmc-iterations 1b,1m \
  --thread 1 \
  --enable-supersites \
  --progress || { echo "supersite run failed"; exit 1; }

log "Running supersite smoke test without --enable-supersites"
"$BIN" \
  --input "$DATA/supersite.vcf.gz" \
  --reference "$DATA/supersite_ref.vcf.gz" \
  --region 1:1-200 \
  --output "$OUTDIR/supersite_disable.graph" \
  --output-format graph \
  --mcmc-iterations 1b,1m \
  --thread 1 \
  --progress || { echo "supersite run (disabled) failed"; exit 1; }

log "OK"

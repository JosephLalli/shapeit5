#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "usage: $0 chr22:19300000-19400000 [--keep-phased]" >&2
  exit 1
fi

REGION="$1"
KEEP_PHASED=0
if [[ $# -eq 2 ]]; then
  if [[ "$2" == "--keep-phased" ]]; then
    KEEP_PHASED=1
  else
    echo "unknown option: $2" >&2
    exit 1
  fi
fi
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${OUT_DIR:-$ROOT_DIR/tmp_multi/bench_leftshift}"
MAP_FILE="${MAP_FILE:-$ROOT_DIR/test/info/chr22.gmap.gz}"
REF_FASTA="${REF_FASTA:-/mnt/ssd/lalli/phasing_T2T/chm13v2.0.fa}"
THREADS_PHASE="${THREADS_PHASE:-64}"
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

MULTIALLELIC_UNPHASED="$ROOT_DIR/tmp_multi/1KGP.CHM13v2.0.chr22.refiltered_from_scratch.multi.unphased.bcf"
TRUTH_MULTIALLELIC="$ROOT_DIR/tmp_multi/hprc-v2.0-mc-chm13.multi.filtered.phased.bcf"
TRUTH_BIALLELIC="$ROOT_DIR/tmp_multi/hprc-v2.0-mc-chm13.multi.filtered.biallelic.phased.bcf"

require_file "$MULTIALLELIC_UNPHASED"
require_file "$TRUTH_MULTIALLELIC"
require_file "$TRUTH_BIALLELIC"
require_file "$MAP_FILE"
require_file "$REF_FASTA"

mkdir -p "$OUT_DIR"

PHASE_MULTIALLELIC="$OUT_DIR/phase.1kgp.multi.noshift.phased.bcf"
MULTIALLELIC_PHASED_SPLIT_BIALLELIC="$OUT_DIR/phase.1kgp.multi.noshift.phased.split.biallelic.bcf"
MULTIALLELIC_UNPHASED_SPLIT_BIALLELIC="$OUT_DIR/1kgp.multi.unphased.split.biallelic.bcf"
PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT="$OUT_DIR/phase.1kgp.multi.unphased.split.biallelic.phased.bcf"

SWITCH_MULTIALLELIC_SPLIT_PREFIX="$OUT_DIR/eval.1kgp.multi.noshift.split_bial.switch"
SWITCH_BIALLELIC_FROM_MULTIALLELIC_SPLIT_PREFIX="$OUT_DIR/eval.1kgp.multi.unphased.split_bial.phased.switch"
SWITCH_MULTIALLELIC_PREFIX="$OUT_DIR/eval.1kgp.multi.noshift.switch"

VT_BIALLELIC="$OUT_DIR/eval.1kgp.multi.unphased.split_bial.phased.switch.variant.typing.txt.gz"
VT_MULTIALLELIC_SPLIT="$OUT_DIR/eval.1kgp.multi.noshift.split_bial.switch.variant.typing.txt.gz"

if [[ "$KEEP_PHASED" -eq 1 && -s "$PHASE_MULTIALLELIC" ]]; then
  echo "Reusing existing phased multiallelic output: $PHASE_MULTIALLELIC"
else
  echo "Phasing multiallelic input"
  "$PHASE_COMMON" \
    --input "$MULTIALLELIC_UNPHASED" \
    --map "$MAP_FILE" \
    --region "$REGION" \
    --output "$PHASE_MULTIALLELIC" \
    --output-format bcf \
    --thread "$THREADS_PHASE" \
    --seed "$SEED" \
    --enable-supersites \
    --filter-maf 0.001 \
    --log "${PHASE_MULTIALLELIC%.bcf}.log" && \
  bcftools index -f "$PHASE_MULTIALLELIC" &
fi

echo "Prep: split unphased multiallelic input and phase biallelic"
bcftools norm -m -any -O b -o "$MULTIALLELIC_UNPHASED_SPLIT_BIALLELIC" "$MULTIALLELIC_UNPHASED"
bcftools index -f "$MULTIALLELIC_UNPHASED_SPLIT_BIALLELIC"

if [[ "$KEEP_PHASED" -eq 1 && -s "$PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT" ]]; then
  echo "Reusing existing phased split-biallelic output: $PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT"
else
  echo "Phasing split-biallelic input"
  "$PHASE_COMMON" \
    --input "$MULTIALLELIC_UNPHASED_SPLIT_BIALLELIC" \
    --map "$MAP_FILE" \
    --region "$REGION" \
    --output "$PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT" \
    --output-format bcf \
    --thread "$THREADS_PHASE" \
    --seed "$SEED" \
    --filter-maf 0.001 \
    --log "${PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT%.bcf}.log" && \
  bcftools index -f "$PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT" &
fi
wait

bcftools norm -m -any -O b -o "$MULTIALLELIC_PHASED_SPLIT_BIALLELIC" "$PHASE_MULTIALLELIC"
bcftools index -f "$MULTIALLELIC_PHASED_SPLIT_BIALLELIC"

"$SWITCH_BIN" \
  --validation "$TRUTH_MULTIALLELIC" \
  --estimation "$PHASE_MULTIALLELIC" \
  --region "$REGION" \
  --output "$SWITCH_MULTIALLELIC_PREFIX" \
  --log "${SWITCH_MULTIALLELIC_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIALLELIC" \
  --estimation "$MULTIALLELIC_PHASED_SPLIT_BIALLELIC" \
  --region "$REGION" \
  --output "$SWITCH_MULTIALLELIC_SPLIT_PREFIX" \
  --log "${SWITCH_MULTIALLELIC_SPLIT_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &

"$SWITCH_BIN" \
  --validation "$TRUTH_BIALLELIC" \
  --estimation "$PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT" \
  --region "$REGION" \
  --output "$SWITCH_BIALLELIC_FROM_MULTIALLELIC_SPLIT_PREFIX" \
  --log "${SWITCH_BIALLELIC_FROM_MULTIALLELIC_SPLIT_PREFIX}.log" \
  --thread "$THREADS_SWITCH" &
#fi
wait

echo "Step 1: summarize switch results and typing deltas"
python3 "$ROOT_DIR/scripts/bench_leftshift_report.py" \
  --out-dir "$OUT_DIR" \
  --vt-biallelic "$VT_BIALLELIC" \
  --vt-multiallelic-split "$VT_MULTIALLELIC_SPLIT"

gt_compare() {
  local a="$1"
  local b="$2"
  local prefix="$3"
  GT_A="$a" GT_B="$b" REGION="$REGION" OUT_PREFIX="$prefix" python - <<'PY'
import os
import re
import subprocess
from pathlib import Path

gt_a = os.environ["GT_A"]
gt_b = os.environ["GT_B"]
region = os.environ["REGION"]
out_prefix = Path(os.environ["OUT_PREFIX"])

def sample_list(path):
    out = subprocess.check_output(["bcftools", "query", "-l", path], text=True)
    return [s for s in out.splitlines() if s]

samples_a = sample_list(gt_a)
samples_b = sample_list(gt_b)
if samples_a != samples_b:
    raise SystemExit("sample order mismatch between estimates")

nsamp = len(samples_a)
per_sample_compared = [0] * nsamp
per_sample_mismatch = [0] * nsamp

def parse_gt(gt):
    if gt in {".", "./.", ".|."}:
        return None
    alleles = re.split(r"[|/]", gt)
    if any(a == "." for a in alleles):
        return None
    if len(alleles) == 0:
        return None
    return tuple(sorted(alleles))

def iter_query(path):
    cmd = ["bcftools", "query", "-r", region, "-f", "%CHROM\t%POS\t%REF\t%ALT[\\t%GT]\\n", path]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    try:
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            key = (parts[0], int(parts[1]), parts[2], parts[3])
            yield key, parts[4:]
    finally:
        proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"bcftools query failed for {path}")

it_a = iter_query(gt_a)
it_b = iter_query(gt_b)

key_a, gts_a = next(it_a, (None, None))
key_b, gts_b = next(it_b, (None, None))

def key_lt(a, b):
    return a < b

variants_common = 0
variants_only_a = 0
variants_only_b = 0
compared = 0
mismatches = 0

while key_a is not None and key_b is not None:
    if key_a == key_b:
        variants_common += 1
        if len(gts_a) != nsamp or len(gts_b) != nsamp:
            raise SystemExit("genotype column count mismatch")
        for i in range(nsamp):
            ga = parse_gt(gts_a[i])
            gb = parse_gt(gts_b[i])
            if ga is None or gb is None:
                continue
            per_sample_compared[i] += 1
            compared += 1
            if ga != gb:
                per_sample_mismatch[i] += 1
                mismatches += 1
        key_a, gts_a = next(it_a, (None, None))
        key_b, gts_b = next(it_b, (None, None))
    elif key_lt(key_a, key_b):
        variants_only_a += 1
        key_a, gts_a = next(it_a, (None, None))
    else:
        variants_only_b += 1
        key_b, gts_b = next(it_b, (None, None))

variants_only_a += sum(1 for _ in it_a)
variants_only_b += sum(1 for _ in it_b)

summary_path = out_prefix.with_suffix(".summary.txt")
per_sample_path = out_prefix.with_suffix(".per_sample.tsv")

with summary_path.open("w") as out:
    out.write(f"variants_common\t{variants_common}\n")
    out.write(f"variants_only_a\t{variants_only_a}\n")
    out.write(f"variants_only_b\t{variants_only_b}\n")
    out.write(f"genotypes_compared\t{compared}\n")
    out.write(f"genotypes_mismatched\t{mismatches}\n")
    rate = (mismatches / compared) if compared else 0.0
    out.write(f"genotype_mismatch_rate\t{rate:.6f}\n")

with per_sample_path.open("w") as out:
    out.write("sample\tcompared\tmismatched\trate\n")
    for s, c, m in zip(samples_a, per_sample_compared, per_sample_mismatch):
        r = (m / c) if c else 0.0
        out.write(f"{s}\t{c}\t{m}\t{r:.6f}\n")

print(f"wrote {summary_path}")
print(f"wrote {per_sample_path}")
PY
}

#echo "Step 2: genotype mismatch between biallelic vs multiallelic_split estimates"
#gt_compare "$PHASE_BIALLELIC_FROM_MULTIALLELIC_SPLIT" "$MULTIALLELIC_PHASED_SPLIT_BIALLELIC" "$OUT_DIR/compare.gt_mismatch.biallelic_vs_multiallelic_split"

echo "Step 3: evaluate both pipelines against biallelic truth"
echo "  multiallelic->split: ${SWITCH_MULTIALLELIC_SPLIT_PREFIX}.log"
echo "  split->biallelic:    ${SWITCH_BIALLELIC_FROM_MULTIALLELIC_SPLIT_PREFIX}.log"

echo "done"

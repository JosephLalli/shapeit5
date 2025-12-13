#!/usr/bin/env bash
set -euo pipefail

# Triage per-variant SER deltas from a single .tsv file produced by the
# supersite comparison harness.
#
# Focus: variants where supersites increase the switch error rate vs "main".
# Uses column 15: d_SER_super_main.
#
# Usage:
#   scripts/triage_ser_deltas.sh path/to/*.ser_deltas.tsv [min_delta] [min_nb] [max_per_group]
#
# Defaults:
#   min_delta     = 5   (minimum absolute increase in SER, percentage points)
#   min_nb        = 5   (minimum nb_super, to avoid tiny denominators)
#   max_per_group = 20  (max variants shown per sup_status group)

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <ser_deltas.tsv> [min_delta] [min_nb] [max_per_group]" >&2
  exit 1
fi

FILE=$1
MIN_DELTA=${2:-5}
MIN_NB=${3:-5}
MAX_PER_GROUP=${4:-20}

if [[ ! -f "$FILE" ]]; then
  echo "Error: file not found: $FILE" >&2
  exit 1
fi

HEADER=$(head -n1 "$FILE")

echo "Triage of SER deltas file:"
echo "  file        : $FILE"
echo "  min_delta   : $MIN_DELTA (d_SER_super_main, percentage points)"
echo "  min_nb      : $MIN_NB (nb_super)"
echo "  max_per_grp : $MAX_PER_GROUP"
echo

echo "=== Summary: variants with d_SER_super_main > 0 ==="
awk '
  NR > 1 && $15 != "" && ($15 + 0) > 0 {
    total++;
    by_status[$13]++;
  }
  END {
    printf("total_positive\t%d\n", total + 0);
    for (s in by_status) {
      printf("%s\t%d\n", s, by_status[s]);
    }
  }
' "$FILE"

echo
echo "=== Positive variants by sup_status and var_count ==="
awk '
  NR > 1 && $15 != "" && ($15 + 0) > 0 {
    key = $13 "\tvar_count=" $14;
    by_key[key]++;
  }
  END {
    for (k in by_key) {
      printf("%s\t%d\n", k, by_key[k]);
    }
  }
' "$FILE" | sort

show_group() {
  local status=$1
  echo
  echo "=== Top variants (sup_status=${status}, d_SER_super_main >= ${MIN_DELTA}, nb_super >= ${MIN_NB}) ==="
  echo "$HEADER"
  awk -v status="$status" -v min_delta="$MIN_DELTA" -v min_nb="$MIN_NB" '
    NR > 1 && $13 == status && $15 != "" && ($15 + 0) >= min_delta && ($4 + 0) >= min_nb {
      print
    }
  ' "$FILE" | sort -k15,15nr -k4,4nr | head -n "$MAX_PER_GROUP"
}

show_group "anchor"
show_group "sibling"
show_group "non"


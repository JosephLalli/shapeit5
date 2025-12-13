#!/usr/bin/env bash
set -euo pipefail

SUPER="$1"
MAIN="$2"
OG="$3"
OUT="$4"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Normalize a switch file → id pos SER
normalize_switch() {
    local in="$1"
    local label="$2"
    local out="$3"

    zcat "$in" \
    | awk -v lbl="$label" 'BEGIN{FS=OFS=" "}
        {
            id=$1; pos=$2
            # Detect presence of SER (5th col) or compute from col3/col4
            if (NF >= 5) {
                ser=$5
            } else {
                # Fallback: compute SER = 100 * n_switch_err / n_biallelic
                nse=$3; nb=$4
                if (nb==0) ser="nan"; else ser=100*(nse/nb)
            }
            print id, pos, ser
        }' \
    > "$out"
}

normalize_switch "$SUPER" "supersites" "$tmpdir/super.tsv"
normalize_switch "$MAIN"  "main"       "$tmpdir/main.tsv"
normalize_switch "$OG"    "og"         "$tmpdir/og.tsv"

# Convert "id pos ser" → "id|pos ser"
awk 'BEGIN{FS=OFS=" "} {print $1"|" $2, $3}' "$tmpdir/super.tsv" > "$tmpdir/super.key"
awk 'BEGIN{FS=OFS=" "} {print $1"|" $2, $3}' "$tmpdir/main.tsv"  > "$tmpdir/main.key"
awk 'BEGIN{FS=OFS=" "} {print $1"|" $2, $3}' "$tmpdir/og.tsv"    > "$tmpdir/og.key"

# Join supersites + main
join -a1 -a2 -e "NA" -o 0,1.2,2.2 \
    -1 1 -2 1 \
    <(sort "$tmpdir/super.key") \
    <(sort "$tmpdir/main.key") \
    > "$tmpdir/join1"

# Join with og
join -a1 -a2 -e "NA" -o 0,1.2,2.2 \
    -1 1 -2 1 \
    <(sort "$tmpdir/join1") \
    <(sort "$tmpdir/og.key") \
    > "$tmpdir/join2"

# Output format:
# join2: "id|pos SER_super SER_main SER_og"

awk 'BEGIN{FS=OFS=" "}
    {
        split($1, kp, "|")
        id=kp[1]; pos=kp[2]
        s_super=$2; s_main=$3; s_og=$4

        # Cast "NA" to numeric NA
        d_sm="NA"; d_so="NA"; d_mo="NA"
        if(s_super!="NA" && s_main!="NA") d_sm=s_super - s_main
        if(s_super!="NA" && s_og!="NA")   d_so=s_super - s_og
        if(s_main!="NA"  && s_og!="NA")   d_mo=s_main  - s_og

        print id, pos, s_super, s_main, s_og, d_sm, d_so, d_mo
    }' "$tmpdir/join2" \
    > "$OUT"

echo "Wrote $OUT"

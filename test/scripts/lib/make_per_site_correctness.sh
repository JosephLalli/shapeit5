#!/usr/bin/env bash
set -euo pipefail

if [[ "$#" -ne 6 ]]; then
    echo "Usage: $0 TRUTH_VCF RUN_VCF REF_FASTA TRUTH_BED OUTDIR RUN_NAME" >&2
    echo "  TRUTH_VCF  : phased truth VCF/BCF (multi-sample or single-sample)" >&2
    echo "  RUN_VCF    : phased query VCF/BCF for this run (bial or supersite)" >&2
    echo "  REF_FASTA  : reference FASTA used by vcfdist" >&2
    echo "  TRUTH_BED  : BED file describing evaluation regions" >&2
    echo "  OUTDIR     : output directory" >&2
    echo "  RUN_NAME   : label for this run (e.g. bial, supersite)" >&2
    exit 1
fi

TRUTH_VCF="$1"
RUN_VCF="$2"
REF_FASTA="$3"
TRUTH_BED="$4"
OUTDIR="$5"
RUN_NAME="$6"

mkdir -p "${OUTDIR}"

# Output file for this run
OUT_TSV="${OUTDIR}/${RUN_NAME}.per_site.tsv"
echo -e "chrom\tpos\tref\talt\tsample\tis_correct" > "${OUT_TSV}"

# Make sure vcfdist and bcftools exist
command -v bcftools >/dev/null 2>&1 || { echo "ERROR: bcftools not in PATH" >&2; exit 1; }
command -v vcfdist  >/dev/null 2>&1 || { echo "ERROR: vcfdist not in PATH" >&2; exit 1; }

# List samples in the run VCF
SAMPLES=$(bcftools query -l "${RUN_VCF}")

for SAMPLE in ${SAMPLES}; do
    echo "Processing sample ${SAMPLE} for run ${RUN_NAME}..." >&2

    SAMPLE_DIR="${OUTDIR}/${RUN_NAME}/${SAMPLE}"
    mkdir -p "${SAMPLE_DIR}"

    QUERY_VCF="${SAMPLE_DIR}/query.${SAMPLE}.vcf.gz"
    TRUTH_SAMPLE_VCF="${SAMPLE_DIR}/truth.${SAMPLE}.vcf.gz"

    # Extract this sample from query and truth VCFs
    bcftools view -s "${SAMPLE}" -Oz -o "${QUERY_VCF}" "${RUN_VCF}"
    tabix -p vcf "${QUERY_VCF}"

    bcftools view -s "${SAMPLE}" -Oz -o "${TRUTH_SAMPLE_VCF}" "${TRUTH_VCF}"
    tabix -p vcf "${TRUTH_SAMPLE_VCF}"

    # Run vcfdist for this sample
    # Outputs (including query.tsv) will be placed directly in SAMPLE_DIR
    vcfdist \
        "${QUERY_VCF}" \
        "${TRUTH_SAMPLE_VCF}" \
        "${REF_FASTA}" \
        -b "${TRUTH_BED}" \
        -p "${SAMPLE_DIR}/" \
        -v 0 >/dev/null 2>&1

    QUERY_TSV="${SAMPLE_DIR}/query.tsv"
    if [[ ! -f "${QUERY_TSV}" ]]; then
        echo "ERROR: Expected ${QUERY_TSV} from vcfdist for sample ${SAMPLE}" >&2
        exit 1
    fi

    # Collapse query.tsv to per-site correctness:
    #   - group by (CONTIG, POS)
    #   - if any CREDIT > 0 at that position, mark is_correct=1, else 0
    #
    # query.tsv columns (from vcfdist docs) are:
    #   1: CONTIG
    #   2: POS (0-based)
    #   3: HAP
    #   4: REF
    #   5: ALT
    #   6: QUAL
    #   7: TYPE
    #   8: ERR_TYPE
    #   9: CREDIT
    #   ...
    #
    # We sort by chrom,pos and then do a streaming group-by in awk.
    #
    # NOTE: vcfdist uses 0-based POS; if you want 1-based, add +1 in awk.

    tail -n +2 "${QUERY_TSV}" \
    | sort -k1,1 -k2,2n -k4,4 -k5,5 \
    | awk -v sample="${SAMPLE}" '
            BEGIN {
                FS = OFS = "\t";
                prev_key = "";
                max_credit = 0;
                have_prev = 0;
            }
            {
                chrom = $1;
                pos   = $2;       # 0-based POS from vcfdist
                ref   = $4;
                alt   = $5;
                credit = $9 + 0;  # CREDIT column

                key = chrom "_" pos "_" ref "_" alt;

                # Emit previous group if key changed
                if (have_prev && key != prev_key) {
                    is_correct = (max_credit > 0.0 ? 1 : 0);

                    # variant_id = chrom:pos:ref:alt (1-based POS optional; choose consistent)
                    variant_id = prev_chrom "_" (prev_pos+1) "_" prev_ref "_" prev_alt;

                    print prev_chrom, prev_pos+1, prev_ref, prev_alt, sample, is_correct;

                    max_credit = 0;
                }

                # Update group accumulation
                if (credit > max_credit) {
                    max_credit = credit;
                }

                prev_key   = key;
                prev_chrom = chrom;
                prev_pos   = pos;
                prev_ref   = ref;
                prev_alt   = alt;
                have_prev  = 1;
            }
            END {
                if (have_prev) {
                    is_correct = (max_credit > 0.0 ? 1 : 0);
                    variant_id = prev_chrom "_" (prev_pos+1) "_" prev_ref "_" prev_alt;

                    print prev_chrom, prev_pos+1, prev_ref, prev_alt, sample, is_correct;
                }
            }
    ' >> "${OUT_TSV}"

done

echo "Done. Per-site correctness for run '${RUN_NAME}' written to:"
echo "  ${OUT_TSV}"

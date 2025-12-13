#!/usr/bin/env python3
"""
Enforce supersite missingness rule:
  - Group records by (contig, position) supersite.
  - For each sample, if any variant in the supersite is missing and all others
    are 0/0 (or also missing), rewrite every variant in the supersite for that
    sample to missing ("./.").
  - If a supersite mixes a missing call with any 0/1 or 1/1 for the same
    sample, raise an exception.
"""

import argparse
import sys
from typing import Iterable, List, Tuple

import pysam


def is_missing(gt: Tuple[int, ...] | None) -> bool:
    return gt is None or all(allele is None for allele in gt)


def has_alt(gt: Tuple[int, ...] | None) -> bool:
    return gt is not None and any(allele not in (None, 0) for allele in gt)


def is_ref(gt: Tuple[int, ...] | None) -> bool:
    return gt is not None and all(allele == 0 for allele in gt)


def detect_ploidy(genotypes: Iterable[Tuple[int, ...] | None]) -> int:
    for gt in genotypes:
        if gt is not None:
            return len(gt)
    return 2


def process_group(records: List[pysam.VariantRecord], out_vcf: pysam.VariantFile) -> None:
    if not records:
        return

    samples = records[0].samples.keys()

    for sample in samples:
        sample_gts = [rec.samples[sample].get("GT") for rec in records]
        missing_seen = any(is_missing(gt) for gt in sample_gts)
        alt_seen = any(has_alt(gt) for gt in sample_gts)

        if not missing_seen:
            continue

        if alt_seen:
            raise RuntimeError(
                f"Found missing + alt mix at supersite {records[0].chrom}:{records[0].pos} "
                f"for sample {sample} with GTs {sample_gts}"
            )

        if all(is_missing(gt) or is_ref(gt) for gt in sample_gts):
            ploidy = detect_ploidy(sample_gts)
            missing_gt = tuple(None for _ in range(ploidy))
            for rec in records:
                rec.samples[sample]["GT"] = missing_gt

    for rec in records:
        out_vcf.write(rec)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rewrite supersite genotypes so partial missingness becomes fully missing per sample."
    )
    parser.add_argument("input", help="Input VCF/BCF (must be position-sorted).")
    parser.add_argument(
        "output",
        help="Output VCF/BCF path. Use '-' for stdout.",
    )
    args = parser.parse_args()

    in_vcf = pysam.VariantFile(args.input)
    out_vcf = pysam.VariantFile(args.output, "w", header=in_vcf.header)

    current_group: List[pysam.VariantRecord] = []
    last_key: Tuple[str, int] | None = None

    for record in in_vcf:
        key = (record.chrom, record.pos)
        if last_key is None or key == last_key:
            current_group.append(record)
            last_key = key
            continue

        process_group(current_group, out_vcf)
        current_group = [record]
        last_key = key

    process_group(current_group, out_vcf)

    out_vcf.close()
    in_vcf.close()


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001 - surface context directly
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

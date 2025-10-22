#!/usr/bin/env python3
import argparse
from cyvcf2 import VCF

def hap_state(gt):
    """
    Return which haplotype(s) carry ALT for a phased biallelic GT.
    - {0}  -> ALT on left hap (e.g., 1|0)
    - {1}  -> ALT on right hap (e.g., 0|1)
    - {0,1}-> both haps ALT (e.g., 1|1)
    - None -> unphased / missing / not biallelic 0/1
    - set() -> no ALT (0|0)
    """
    if gt is None or len(gt) < 3:
        return None
    a, b, phased = gt[0], gt[1], bool(gt[2])
    if not phased:
        return None
    if a < 0 or b < 0:
        return None
    # Expect split biallelic => alleles in {0,1}
    if a not in (0,1) or b not in (0,1):
        return None
    if a == 0 and b == 0:
        return set()
    if a == 1 and b == 0:
        return {0}
    if a == 0 and b == 1:
        return {1}
    # 1|1
    return {0,1}

def process_group(records, sample_names, include_filtered=False, outfh=None):
    """
    records: list of cyvcf2.Variant at the same CHROM+POS
    For each sample, track which ALT(s) appear on each hap. Emit conflicts where
    a hap has >1 distinct ALT at this position.
    """
    if not records:
        return 0
    # Map: sample_index -> {0: set(ALTs), 1: set(ALTs)}
    conflicts = 0
    hap_alts = [ {0:set(), 1:set()} for _ in sample_names ]

    # Preload per-record ALT strings (first ALT only; split-biallelic)
    rec_alts = []
    for v in records:
        if (not include_filtered) and v.FILTER not in (None, 'PASS', '.'):
            rec_alts.append(None)  # treated as skipped
            continue
        alt = (v.ALT or [None])[0]
        rec_alts.append(alt)

    for r_idx, v in enumerate(records):
        alt = rec_alts[r_idx]
        if not alt:  # skip filtered/invalid
            continue
        gts = v.genotypes  # list of [a,b,phased,(GQ...)]
        for i, gt in enumerate(gts):
            hs = hap_state(gt)
            if hs is None or len(hs) == 0:
                continue
            if hs == {0,1}:
                # ALT on both haps => add to both
                for h in (0,1):
                    hap_alts[i][h].add(alt)
            else:
                h = next(iter(hs))
                hap_alts[i][h].add(alt)

    # After collecting, check each sample for conflicts (hap set size > 1)
    for i, per_hap in enumerate(hap_alts):
        for h in (0,1):
            if len(per_hap[h]) > 1:
                # Conflict: same hap has multiple different ALTs at this position
                sname = sample_names[i]
                chrom = records[0].CHROM
                pos = records[0].POS
                alts_sorted = sorted(per_hap[h])
                outfh.write(f"{chrom}\t{pos}\t{sname}\thap{h}\t" +
                            ",".join(alts_sorted) + "\n")
                conflicts += 1
    return conflicts

def iterate_positions(vcf):
    """
    Generator yielding lists of Variant objects grouped by CHROM+POS.
    """
    curr = []
    last_key = None
    for v in vcf:
        key = (v.CHROM, v.POS)
        if last_key is None:
            curr = [v]
            last_key = key
            continue
        if key == last_key:
            curr.append(v)
        else:
            yield curr
            curr = [v]
            last_key = key
    if curr:
        yield curr

def main():
    ap = argparse.ArgumentParser(description="Detect per-sample, per-position haplotype conflicts across split-biallelic records (two different ALTs on the same hap).")
    ap.add_argument("vcf", help="Input phased, split-biallelic VCF/BCF (.vcf.gz/.bcf), indexed")
    ap.add_argument("-o", "--out", default="-", help="Output TSV (default: stdout)")
    ap.add_argument("--include-filtered", action="store_true",
                    help="Include records with non-PASS FILTER")
    args = ap.parse_args()

    v = VCF(args.vcf)
    samples = v.samples

    outfh = open(args.out, "w") if args.out != "-" else __import__("sys").stdout
    outfh.write("#CHROM\tPOS\tSAMPLE\tHAP\tALTS_AT_SAME_POS\n")

    total_conflicts = 0
    for group in iterate_positions(v):
        if len(group) < 2:
            continue  # only meaningful if multiple records at same position
        total_conflicts += process_group(group, samples, args.include_filtered, outfh)

    if outfh is not __import__("sys").stdout:
        outfh.close()

    # Print a quick summary to stderr
    import sys
    print(f"[summary] conflicts reported: {total_conflicts}", file=sys.stderr)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Enhanced violation analysis script to identify transition mode cases.

Extends find_one_alt_per_hap_violation.py to categorize violation types:
- trans_case: One haplotype has exactly 2 ALTs, other haplotype has 0 ALTs
- mixed_case: Both haplotypes have ALTs (various combinations)
- complex_case: One haplotype has >2 ALTs

This helps understand why transition mode only fixes ~200 out of 8600+ violations.
"""

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

def classify_violation_type(hap0_alts, hap1_alts):
    """
    Classify the type of violation based on ALT distribution across haplotypes.
    
    Returns:
    - 'trans_case': One haplotype has exactly 2 ALTs, other has 0 ALTs (transition mode should fix)
    - 'mixed_case': Both haplotypes have ALTs 
    - 'complex_case': One haplotype has >2 ALTs
    - 'no_violation': No violation detected
    """
    h0_count = len(hap0_alts)
    h1_count = len(hap1_alts)
    
    # No violation
    if h0_count <= 1 and h1_count <= 1:
        return 'no_violation'
    
    # Transition cases: exactly 2 ALTs on one hap, 0 on the other
    if h0_count == 2 and h1_count == 0:
        return 'trans_case'
    if h1_count == 2 and h0_count == 0:
        return 'trans_case'
    
    # Complex cases: >2 ALTs on one haplotype
    if h0_count > 2 or h1_count > 2:
        return 'complex_case'
    
    # Mixed cases: both haplotypes have ALTs (any other combination)
    return 'mixed_case'

def process_group(records, sample_names, include_filtered=False, outfh=None, stats=None):
    """
    records: list of cyvcf2.Variant at the same CHROM+POS
    For each sample, track which ALT(s) appear on each hap. Emit conflicts where
    a hap has >1 distinct ALT at this position, with violation type classification.
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

    # After collecting, check each sample for conflicts and classify violation types
    for i, per_hap in enumerate(hap_alts):
        hap0_alts = per_hap[0]
        hap1_alts = per_hap[1]
        
        # Check for violations on either haplotype
        violation_found = False
        for h in (0,1):
            if len(per_hap[h]) > 1:
                violation_found = True
                break
        
        if violation_found:
            # Classify the violation type
            violation_type = classify_violation_type(hap0_alts, hap1_alts)
            
            # Update statistics
            if stats is not None:
                stats[violation_type] = stats.get(violation_type, 0) + 1
            
            # Output the violating haplotype(s)
            for h in (0,1):
                if len(per_hap[h]) > 1:
                    sname = sample_names[i]
                    chrom = records[0].CHROM
                    pos = records[0].POS
                    alts_sorted = sorted(per_hap[h])
                    other_hap_count = len(per_hap[1-h])
                    other_hap_alts = ",".join(sorted(per_hap[1-h])) if per_hap[1-h] else "none"
                    
                    outfh.write(f"{chrom}\t{pos}\t{sname}\thap{h}\t" +
                                ",".join(alts_sorted) + f"\t{violation_type}\t{other_hap_count}\t{other_hap_alts}\n")
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
    ap = argparse.ArgumentParser(description="Analyze violation patterns and identify transition mode cases: violations where one haplotype has exactly 2 ALTs and the other has 0 ALTs.")
    ap.add_argument("vcf", help="Input phased, split-biallelic VCF/BCF (.vcf.gz/.bcf), indexed")
    ap.add_argument("-o", "--out", default="-", help="Output TSV (default: stdout)")
    ap.add_argument("--include-filtered", action="store_true",
                    help="Include records with non-PASS FILTER")
    args = ap.parse_args()

    v = VCF(args.vcf)
    samples = v.samples

    outfh = open(args.out, "w") if args.out != "-" else __import__("sys").stdout
    outfh.write("#CHROM\tPOS\tSAMPLE\tHAP\tALTS_AT_SAME_POS\tVIOLATION_TYPE\tOTHER_HAP_ALT_COUNT\tOTHER_HAP_ALTS\n")

    total_conflicts = 0
    violation_stats = {}
    
    for group in iterate_positions(v):
        if len(group) < 2:
            continue  # only meaningful if multiple records at same position
        total_conflicts += process_group(group, samples, args.include_filtered, outfh, violation_stats)

    if outfh is not __import__("sys").stdout:
        outfh.close()

    # Print detailed summary to stderr
    import sys
    print(f"[summary] Total violations reported: {total_conflicts}", file=sys.stderr)
    print(f"[summary] Violation type breakdown:", file=sys.stderr)
    for vtype, count in sorted(violation_stats.items()):
        pct = (count / total_conflicts * 100) if total_conflicts > 0 else 0
        print(f"  {vtype}: {count} ({pct:.1f}%)", file=sys.stderr)
    
    trans_cases = violation_stats.get('trans_case', 0)
    print(f"\n[analysis] Transition mode should theoretically fix {trans_cases} violations", file=sys.stderr)
    print(f"[analysis] This represents {(trans_cases / total_conflicts * 100):.1f}% of all violations", file=sys.stderr)

if __name__ == "__main__":
    main()
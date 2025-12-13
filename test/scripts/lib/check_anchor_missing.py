#!/usr/bin/env python3
"""
Check if incompatible sites have missing anchors or non-missing anchors.
This helps diagnose whether the bug is in anchor_has_missing detection.
"""

import subprocess
import sys
from collections import defaultdict

def parse_gt(gt_str):
    if '|' in gt_str:
        parts = gt_str.split('|')
    elif '/' in gt_str:
        parts = gt_str.split('/')
    else:
        return None, None
    try:
        return int(parts[0]), int(parts[1])
    except (ValueError, IndexError):
        return None, None

def is_missing(gt_str):
    """Check if genotype is missing (./. or .|.)"""
    gt = gt_str.split(':')[0]
    return '.' in gt

def main():
    if len(sys.argv) < 2:
        print("Usage: python check_anchor_missing.py <phased.bcf> [region]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    region = sys.argv[2] if len(sys.argv) > 2 else None

    cmd = ['bcftools', 'view', '-H', bcf_file]
    if region:
        cmd.extend(['-r', region])

    sample_cmd = ['bcftools', 'query', '-l', bcf_file]
    samples = subprocess.check_output(sample_cmd, text=True).strip().split('\n')

    # Group variants by position
    position_data = defaultdict(list)
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        genotypes = fields[9:]
        position_data[pos].append({'ref': ref, 'alt': alt, 'genotypes': genotypes})

    # Find incompatible cases and check anchor status
    anchor_missing_count = 0
    anchor_not_missing_count = 0

    incompatible_anchor_missing = []
    incompatible_anchor_not_missing = []

    for pos in sorted(position_data.keys()):
        variants = position_data[pos]
        if len(variants) <= 1:
            continue

        # First variant is the anchor (lowest index in file = anchor by convention)
        anchor = variants[0]

        for sample_idx, sample_name in enumerate(samples):
            hap0_alts = []
            hap1_alts = []
            anchor_is_missing = is_missing(anchor['genotypes'][sample_idx])

            for var in variants:
                gt_field = var['genotypes'][sample_idx]
                gt_str = gt_field.split(':')[0]
                hap0, hap1 = parse_gt(gt_str)

                if hap0 is None:
                    continue

                if hap0 > 0:
                    hap0_alts.append(var['alt'])
                if hap1 > 0:
                    hap1_alts.append(var['alt'])

            # Check for incompatibility
            if len(hap0_alts) > 1 or len(hap1_alts) > 1:
                if anchor_is_missing:
                    anchor_missing_count += 1
                    if len(incompatible_anchor_missing) < 5:
                        incompatible_anchor_missing.append({
                            'pos': pos,
                            'sample': sample_name,
                            'anchor_gt': anchor['genotypes'][sample_idx].split(':')[0],
                            'hap0': hap0_alts,
                            'hap1': hap1_alts
                        })
                else:
                    anchor_not_missing_count += 1
                    if len(incompatible_anchor_not_missing) < 5:
                        incompatible_anchor_not_missing.append({
                            'pos': pos,
                            'sample': sample_name,
                            'anchor_gt': anchor['genotypes'][sample_idx].split(':')[0],
                            'hap0': hap0_alts,
                            'hap1': hap1_alts
                        })

    print("=" * 70)
    print(f"ANCHOR MISSING STATUS FOR INCOMPATIBLE HAPLOTYPES")
    print("=" * 70)
    print(f"\nTotal incompatible cases where ANCHOR is MISSING: {anchor_missing_count}")
    print(f"Total incompatible cases where ANCHOR is NOT missing: {anchor_not_missing_count}")

    if anchor_missing_count + anchor_not_missing_count > 0:
        pct_missing = anchor_missing_count / (anchor_missing_count + anchor_not_missing_count) * 100
        print(f"\nPercentage with anchor missing: {pct_missing:.1f}%")

    if incompatible_anchor_missing:
        print(f"\nExamples where ANCHOR IS MISSING:")
        for ex in incompatible_anchor_missing:
            print(f"  pos={ex['pos']} sample={ex['sample']} anchor_gt={ex['anchor_gt']}")
            print(f"    hap0: {ex['hap0']}, hap1: {ex['hap1']}")

    if incompatible_anchor_not_missing:
        print(f"\nExamples where ANCHOR IS NOT MISSING:")
        for ex in incompatible_anchor_not_missing:
            print(f"  pos={ex['pos']} sample={ex['sample']} anchor_gt={ex['anchor_gt']}")
            print(f"    hap0: {ex['hap0']}, hap1: {ex['hap1']}")

if __name__ == '__main__':
    main()

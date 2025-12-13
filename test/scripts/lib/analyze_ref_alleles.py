#!/usr/bin/env python3
"""
Analyze reference allele patterns at multiallelic positions.
Identify positions where variants have different reference alleles.
"""

import subprocess
import sys
from collections import defaultdict

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_ref_alleles.py <bcf_file> [region]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    region = sys.argv[2] if len(sys.argv) > 2 else None

    cmd = ['bcftools', 'view', '-H', bcf_file]
    if region:
        cmd.extend(['-r', region])

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Group by position
    position_data = defaultdict(list)
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4]
        position_data[pos].append({'ref': ref, 'alt': alt})

    # Analyze multiallelic positions
    same_ref_positions = 0
    diff_ref_positions = 0

    same_ref_variants = 0
    diff_ref_variants = 0

    diff_ref_examples = []

    for pos in sorted(position_data.keys()):
        variants = position_data[pos]
        if len(variants) <= 1:
            continue

        refs = set(v['ref'] for v in variants)

        if len(refs) == 1:
            same_ref_positions += 1
            same_ref_variants += len(variants)
        else:
            diff_ref_positions += 1
            diff_ref_variants += len(variants)
            if len(diff_ref_examples) < 20:
                diff_ref_examples.append({
                    'pos': pos,
                    'variants': variants,
                    'unique_refs': refs
                })

    total_multi = same_ref_positions + diff_ref_positions

    print("=" * 70)
    print("REFERENCE ALLELE ANALYSIS AT MULTIALLELIC POSITIONS")
    print("=" * 70)
    print(f"\nTotal multiallelic positions: {total_multi}")
    print(f"\nPositions with SAME reference allele: {same_ref_positions} ({same_ref_positions/total_multi*100:.1f}%)")
    print(f"  Total variants: {same_ref_variants}")
    print(f"\nPositions with DIFFERENT reference alleles: {diff_ref_positions} ({diff_ref_positions/total_multi*100:.1f}%)")
    print(f"  Total variants: {diff_ref_variants}")

    print(f"\n{'='*70}")
    print("EXAMPLES OF POSITIONS WITH DIFFERENT REFERENCE ALLELES:")
    print("=" * 70)

    for ex in diff_ref_examples:
        print(f"\nPosition {ex['pos']} - {len(ex['unique_refs'])} unique refs:")
        for v in ex['variants']:
            print(f"  {v['ref']} > {v['alt']}")

if __name__ == '__main__':
    main()

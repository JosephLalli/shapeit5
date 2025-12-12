#!/usr/bin/env python3
"""
Debug incompatible haplotype sites in supersite output.
Output detailed information about each incompatible case.
"""

import subprocess
import sys
from collections import defaultdict

def parse_gt(gt_str):
    """Parse genotype string like '0|1' into (hap0_allele, hap1_allele)"""
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

def main():
    if len(sys.argv) < 2:
        print("Usage: python debug_incompatible_sites.py <phased.bcf> [region] [max_examples]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    region = sys.argv[2] if len(sys.argv) > 2 else None
    max_examples = int(sys.argv[3]) if len(sys.argv) > 3 else 50

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
        chrom = fields[0]
        pos = int(fields[1])
        var_id = fields[2]
        ref = fields[3]
        alt = fields[4]
        genotypes = fields[9:]

        position_data[pos].append({
            'chrom': chrom,
            'ref': ref,
            'alt': alt,
            'var_id': var_id,
            'genotypes': genotypes
        })

    # Find and report incompatible cases
    incompatible_cases = []
    positions_by_variant_count = defaultdict(int)

    for pos in sorted(position_data.keys()):
        variants = position_data[pos]
        if len(variants) <= 1:
            continue

        positions_by_variant_count[len(variants)] += 1

        for sample_idx, sample_name in enumerate(samples):
            hap0_details = []  # (alt, gt_string)
            hap1_details = []

            for var in variants:
                gt_field = var['genotypes'][sample_idx]
                gt_str = gt_field.split(':')[0]
                hap0, hap1 = parse_gt(gt_str)

                if hap0 is None:
                    continue

                if hap0 > 0:
                    hap0_details.append((var['alt'], gt_str, var['ref']))
                if hap1 > 0:
                    hap1_details.append((var['alt'], gt_str, var['ref']))

            # Check for incompatibility
            if len(hap0_details) > 1 or len(hap1_details) > 1:
                incompatible_cases.append({
                    'pos': pos,
                    'sample': sample_name,
                    'sample_idx': sample_idx,
                    'hap0': hap0_details,
                    'hap1': hap1_details,
                    'num_variants': len(variants),
                    'variants': variants
                })

    # Summary statistics
    print("=" * 70)
    print(f"INCOMPATIBLE HAPLOTYPE ANALYSIS: {bcf_file}")
    print("=" * 70)
    print(f"\nTotal incompatible (sample, position) pairs: {len(incompatible_cases)}")

    # Count by number of variants at position
    incompat_by_var_count = defaultdict(int)
    for case in incompatible_cases:
        incompat_by_var_count[case['num_variants']] += 1

    print(f"\nMultiallelic positions by variant count:")
    for count in sorted(positions_by_variant_count.keys()):
        incompat = incompat_by_var_count.get(count, 0)
        print(f"  {count} variants: {positions_by_variant_count[count]} positions, {incompat} incompatible cases")

    # Count unique positions with incompatibility
    unique_positions = set(case['pos'] for case in incompatible_cases)
    print(f"\nUnique positions with incompatibility: {len(unique_positions)}")

    # Categorize by type of incompatibility
    both_haps_multi = 0  # Both haplotypes have multiple alts
    hap0_only_multi = 0  # Only hap0 has multiple alts
    hap1_only_multi = 0  # Only hap1 has multiple alts

    for case in incompatible_cases:
        h0_multi = len(case['hap0']) > 1
        h1_multi = len(case['hap1']) > 1
        if h0_multi and h1_multi:
            both_haps_multi += 1
        elif h0_multi:
            hap0_only_multi += 1
        else:
            hap1_only_multi += 1

    print(f"\nIncompatibility type breakdown:")
    print(f"  Both haplotypes have >1 ALT: {both_haps_multi}")
    print(f"  Only hap0 has >1 ALT: {hap0_only_multi}")
    print(f"  Only hap1 has >1 ALT: {hap1_only_multi}")

    # Show detailed examples
    print(f"\n{'='*70}")
    print(f"DETAILED EXAMPLES (first {max_examples}):")
    print("=" * 70)

    for i, case in enumerate(incompatible_cases[:max_examples]):
        print(f"\n--- Case {i+1}: pos={case['pos']}, sample={case['sample']} ---")
        print(f"Number of variants at this position: {case['num_variants']}")

        print("\nAll variants at this position:")
        for var in case['variants']:
            gt_field = var['genotypes'][case['sample_idx']]
            gt_str = gt_field.split(':')[0]
            print(f"  {var['ref']}>{var['alt']}  GT={gt_str}")

        print(f"\nHap0 ALTs ({len(case['hap0'])}): ", end="")
        for alt, gt, ref in case['hap0']:
            print(f"{ref}>{alt}({gt}) ", end="")
        print()

        print(f"Hap1 ALTs ({len(case['hap1'])}): ", end="")
        for alt, gt, ref in case['hap1']:
            print(f"{ref}>{alt}({gt}) ", end="")
        print()

    # Group by position to see which positions are most problematic
    pos_counts = defaultdict(int)
    for case in incompatible_cases:
        pos_counts[case['pos']] += 1

    print(f"\n{'='*70}")
    print("TOP 20 POSITIONS BY INCOMPATIBILITY COUNT:")
    print("=" * 70)
    for pos, count in sorted(pos_counts.items(), key=lambda x: -x[1])[:20]:
        num_vars = len(position_data[pos])
        print(f"  pos={pos}: {count} incompatible samples, {num_vars} variants")

if __name__ == '__main__':
    main()

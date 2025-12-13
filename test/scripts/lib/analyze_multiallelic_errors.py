#!/usr/bin/env python3
"""
Analyze switch errors at multiallelic positions to determine how many
are caused by haplotype incompatibility vs other phasing errors.

This script:
1. Identifies multiallelic positions with switch errors
2. Checks if those samples have incompatible haplotypes at those positions
3. Reports overlap between switch errors and incompatibility errors
"""

import subprocess
import sys
from collections import defaultdict
import gzip

def parse_gt(gt_str):
    """Parse genotype string like '0|1' or '1|0' into (hap0_allele, hap1_allele)"""
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

def get_incompatible_samples(bcf_file, region=None):
    """Get set of (position, sample) tuples with incompatible haplotypes"""
    cmd = ['bcftools', 'view', '-H', bcf_file]
    if region:
        cmd.extend(['-r', region])

    sample_cmd = ['bcftools', 'query', '-l', bcf_file]
    samples = subprocess.check_output(sample_cmd, text=True).strip().split('\n')

    position_data = defaultdict(list)
    result = subprocess.run(cmd, capture_output=True, text=True)

    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        pos = int(fields[1])
        alt = fields[4]
        genotypes = fields[9:]
        position_data[pos].append({'alt': alt, 'genotypes': genotypes})

    incompatible = set()
    for pos in position_data:
        variants = position_data[pos]
        if len(variants) <= 1:
            continue

        for sample_idx, sample_name in enumerate(samples):
            hap0_alts = []
            hap1_alts = []

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

            if len(hap0_alts) > 1 or len(hap1_alts) > 1:
                incompatible.add((pos, sample_name))

    return incompatible

def get_switch_errors_at_multiallelic(switch_file, multiallelic_positions):
    """Get switch errors that occur at multiallelic positions"""
    # The switch file format is: variant_id position errors hets rate
    # We need the sample-level switch file or derive from variant-level

    # Actually, we need the flipsAndSwitches file which has sample-level detail
    # Let's read the variant.switch file and identify positions with errors

    errors_at_multiallelic = {}

    opener = gzip.open if switch_file.endswith('.gz') else open
    with opener(switch_file, 'rt') as f:
        for line in f:
            fields = line.strip().split()
            if len(fields) < 3:
                continue
            variant_id = fields[0]
            pos = int(fields[1])
            errors = int(fields[2])

            if pos in multiallelic_positions and errors > 0:
                errors_at_multiallelic[pos] = errors_at_multiallelic.get(pos, 0) + errors

    return errors_at_multiallelic

def main():
    if len(sys.argv) < 3:
        print("Usage: python analyze_multiallelic_errors.py <phased.bcf> <variant.switch.txt.gz> [region]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    switch_file = sys.argv[2]
    region = sys.argv[3] if len(sys.argv) > 3 else None

    print("Getting incompatible samples...")
    incompatible = get_incompatible_samples(bcf_file, region)

    # Get multiallelic positions
    cmd = ['bcftools', 'view', '-H', bcf_file]
    if region:
        cmd.extend(['-r', region])
    result = subprocess.run(cmd, capture_output=True, text=True)

    pos_counts = defaultdict(int)
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        pos = int(line.split('\t')[1])
        pos_counts[pos] += 1

    multiallelic_positions = {pos for pos, count in pos_counts.items() if count > 1}

    print(f"Found {len(multiallelic_positions)} multiallelic positions")
    print(f"Found {len(incompatible)} incompatible (position, sample) pairs")

    # Get switch errors at multiallelic positions
    errors_at_multi = get_switch_errors_at_multiallelic(switch_file, multiallelic_positions)

    total_errors_at_multi = sum(errors_at_multi.values())
    positions_with_errors = len(errors_at_multi)

    print(f"\nSwitch errors at multiallelic positions: {total_errors_at_multi}")
    print(f"Multiallelic positions with errors: {positions_with_errors}")

    # Check overlap - positions that have both incompatibility AND switch errors
    incompatible_positions = {pos for pos, sample in incompatible}
    positions_with_both = incompatible_positions & set(errors_at_multi.keys())

    errors_at_incompatible_positions = sum(errors_at_multi.get(pos, 0) for pos in positions_with_both)

    print(f"\nPositions with BOTH incompatibility errors AND switch errors: {len(positions_with_both)}")
    print(f"Switch errors at those positions: {errors_at_incompatible_positions}")

    if total_errors_at_multi > 0:
        print(f"Percentage of multiallelic switch errors at incompatible positions: {errors_at_incompatible_positions/total_errors_at_multi*100:.1f}%")

if __name__ == '__main__':
    main()

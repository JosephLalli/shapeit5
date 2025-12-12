#!/usr/bin/env python3
"""
Count haplotype-incompatible errors at multiallelic positions.

A haplotype-incompatible error occurs when two different ALT alleles at the same
genomic position are phased to the same haplotype. This is biologically impossible
since each position can only have one allele per haplotype.

Usage:
    python count_incompatible_haplotypes.py <phased.bcf> [region]
"""

import subprocess
import sys
from collections import defaultdict

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

def main():
    if len(sys.argv) < 2:
        print("Usage: python count_incompatible_haplotypes.py <phased.bcf> [region]")
        sys.exit(1)

    bcf_file = sys.argv[1]
    region = sys.argv[2] if len(sys.argv) > 2 else None

    # Build bcftools command
    cmd = ['bcftools', 'view', '-H', bcf_file]
    if region:
        cmd.extend(['-r', region])

    # Get sample names
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
        ref = fields[3]
        alt = fields[4]

        # Get genotypes (field 9 onwards)
        genotypes = fields[9:]

        position_data[pos].append({
            'ref': ref,
            'alt': alt,
            'genotypes': genotypes
        })

    # Analyze multiallelic positions
    total_multiallelic_positions = 0
    total_incompatible_samples = 0
    total_het_samples_at_multiallelic = 0
    incompatible_details = []

    for pos in sorted(position_data.keys()):
        variants = position_data[pos]
        if len(variants) <= 1:
            continue  # Not multiallelic

        total_multiallelic_positions += 1

        # For each sample, check if multiple ALT alleles are on the same haplotype
        for sample_idx, sample_name in enumerate(samples):
            hap0_alts = []  # ALT alleles on haplotype 0
            hap1_alts = []  # ALT alleles on haplotype 1
            is_het_at_any = False

            for var in variants:
                gt_field = var['genotypes'][sample_idx]
                gt_str = gt_field.split(':')[0]  # Get just GT part
                hap0, hap1 = parse_gt(gt_str)

                if hap0 is None:
                    continue

                # Track if this sample is het at any variant here
                if hap0 != hap1:
                    is_het_at_any = True

                # If hap0 has ALT (allele > 0), record it
                if hap0 > 0:
                    hap0_alts.append(var['alt'])
                if hap1 > 0:
                    hap1_alts.append(var['alt'])

            if is_het_at_any:
                total_het_samples_at_multiallelic += 1

            # Check for incompatibility: >1 ALT allele on same haplotype
            if len(hap0_alts) > 1 or len(hap1_alts) > 1:
                total_incompatible_samples += 1
                incompatible_details.append({
                    'pos': pos,
                    'sample': sample_name,
                    'hap0_alts': hap0_alts,
                    'hap1_alts': hap1_alts
                })

    # Print results
    print("=" * 60)
    print(f"Haplotype Incompatibility Analysis: {bcf_file}")
    print("=" * 60)
    print(f"Total multiallelic positions: {total_multiallelic_positions}")
    print(f"Total het sample-positions at multiallelic sites: {total_het_samples_at_multiallelic}")
    print(f"Total incompatible sample-positions: {total_incompatible_samples}")
    if total_het_samples_at_multiallelic > 0:
        rate = total_incompatible_samples / total_het_samples_at_multiallelic * 100
        print(f"Incompatibility rate: {rate:.2f}%")
    print()

    # Show some examples
    if incompatible_details:
        print("First 10 incompatible cases:")
        for detail in incompatible_details[:10]:
            print(f"  pos={detail['pos']} sample={detail['sample']}")
            print(f"    hap0 ALTs: {detail['hap0_alts']}")
            print(f"    hap1 ALTs: {detail['hap1_alts']}")
    else:
        print("No incompatible haplotypes found!")

    return total_incompatible_samples

if __name__ == '__main__':
    main()

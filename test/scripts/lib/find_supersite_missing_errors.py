#!/usr/bin/env python3
"""
Find variants that are:
1. At a supersite position (2+ variants at same bp with MAF > 0.1%)
2. Missing in input for a truth sample
3. Correct in biallelic phased output
4. Incorrect in supersite phased output
"""

import subprocess
import sys
from collections import defaultdict

# File paths
INPUT_BCF = "1KGP.CHM13v2.0.chr22.snp_indel.phasing_qual_pass.biallelic.filtered.bcf"
TRUTH_BCF = "test/wgs/chr22_t2t_reference_pangenome.filtered_variants.18000000-25000000.biallelic.filtered.bcf"
BIAL_BCF = "test/tmp/chr22.1KGP.18-25mb.phase_common.main_algo.small.42.full_unphased.bcf"
SUPER_BCF = "test/tmp/chr22.1KGP.18-25mb.phase_common.supersites.split_emissions.small.42.full_unphased.bcf"
REGION = "chr22:19000000-20000000"

def run_bcftools(cmd):
    """Run bcftools command and return output lines."""
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip().split('\n') if result.stdout.strip() else []

def normalize_gt(gt):
    """Normalize genotype by sorting alleles (ignore phase)."""
    if not gt or gt in ('./.', '.|.', '.'):
        return None
    alleles = sorted(gt.replace('|', '/').split('/'))
    return '/'.join(alleles)

def get_genotypes_for_samples(bcf, region, samples):
    """Get all genotypes at positions in region for given samples.
    Returns dict: (pos, ref, alt, sample) -> gt
    """
    sample_str = ','.join(samples)
    cmd = f"bcftools query -f '%POS\\t%REF\\t%ALT[\\t%GT]\\n' -r {region} -s {sample_str} {bcf} 2>/dev/null"
    lines = run_bcftools(cmd)

    result = {}
    for line in lines:
        if not line:
            continue
        parts = line.split('\t')
        pos = parts[0]
        ref = parts[1]
        alt = parts[2]
        gts = parts[3:]

        for i, gt in enumerate(gts):
            sample = samples[i]
            result[(pos, ref.upper(), alt.upper(), sample)] = gt

    return result

def main():
    print("Step 1: Getting overlap samples...", file=sys.stderr)
    # Get samples in both input and truth
    input_samples = set(run_bcftools(f"bcftools query -l {INPUT_BCF}"))
    truth_samples = set(run_bcftools(f"bcftools query -l {TRUTH_BCF}"))
    overlap_samples = sorted(input_samples & truth_samples)
    print(f"  Found {len(overlap_samples)} overlapping samples", file=sys.stderr)

    print("Step 2: Getting genotypes from all files...", file=sys.stderr)

    print("  Loading input genotypes...", file=sys.stderr)
    input_gts = get_genotypes_for_samples(INPUT_BCF, REGION, overlap_samples)
    print(f"    {len(input_gts)} entries", file=sys.stderr)

    print("  Loading truth genotypes...", file=sys.stderr)
    truth_gts = get_genotypes_for_samples(TRUTH_BCF, REGION, overlap_samples)
    print(f"    {len(truth_gts)} entries", file=sys.stderr)

    print("  Loading biallelic phased genotypes...", file=sys.stderr)
    bial_gts = get_genotypes_for_samples(BIAL_BCF, REGION, overlap_samples)
    print(f"    {len(bial_gts)} entries", file=sys.stderr)

    print("  Loading supersite phased genotypes...", file=sys.stderr)
    super_gts = get_genotypes_for_samples(SUPER_BCF, REGION, overlap_samples)
    print(f"    {len(super_gts)} entries", file=sys.stderr)

    print("Step 3: Finding supersite positions with MAF > 0.1%...", file=sys.stderr)
    # Count variants per position in MAF-filtered input
    cmd = f"bcftools view -q 0.001:minor {INPUT_BCF} -r {REGION} 2>/dev/null | bcftools query -f '%POS\\n'"
    pos_counts = defaultdict(int)
    for line in run_bcftools(cmd):
        if line:
            pos_counts[line] += 1

    supersite_positions = {pos for pos, count in pos_counts.items() if count >= 2}
    print(f"  Found {len(supersite_positions)} supersite positions", file=sys.stderr)

    print("Step 4: Finding candidates...", file=sys.stderr)

    # Print header
    print("POS\tREF\tALT\tSAMPLE\tINPUT_GT\tTRUTH_GT\tBIAL_GT\tSUPER_GT\tBIAL_CORRECT\tSUPER_CORRECT")

    candidates = []
    for key, input_gt in input_gts.items():
        pos, ref, alt, sample = key

        # Check if at supersite position
        if pos not in supersite_positions:
            continue

        # Check if input is missing
        if input_gt not in ('./.', '.|.', '.'):
            continue

        # Get truth GT (case-insensitive match)
        truth_gt = None
        for (t_pos, t_ref, t_alt, t_sample), t_gt in truth_gts.items():
            if t_pos == pos and t_sample == sample and t_ref.upper() == ref and t_alt.upper() == alt:
                truth_gt = t_gt
                break

        # Skip if truth is missing
        if not truth_gt or normalize_gt(truth_gt) is None:
            continue

        # Get phased GTs
        bial_gt = bial_gts.get(key, '')
        super_gt = super_gts.get(key, '')

        if not bial_gt or not super_gt:
            continue

        # Compare
        truth_norm = normalize_gt(truth_gt)
        bial_norm = normalize_gt(bial_gt)
        super_norm = normalize_gt(super_gt)

        bial_correct = 'Y' if bial_norm == truth_norm else 'N'
        super_correct = 'Y' if super_norm == truth_norm else 'N'

        print(f"{pos}\t{ref}\t{alt}\t{sample}\t{input_gt}\t{truth_gt}\t{bial_gt}\t{super_gt}\t{bial_correct}\t{super_correct}")

        # Track candidates where biallelic is correct but supersite is wrong
        if bial_correct == 'Y' and super_correct == 'N':
            candidates.append((pos, ref, alt, sample, input_gt, truth_gt, bial_gt, super_gt))

    print(f"\nFound {len(candidates)} candidates (biallelic correct, supersite incorrect)", file=sys.stderr)

if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""
Comprehensive multiallelic violation analysis script.

Enhanced version that collects detailed metrics about multiallelic violations:
- ALT allele counts per haplotype (total and common-only)
- Number of variants at each genomic position
- Missing genotype counts
- Rare variant identification and counting
- Violation type classification

This provides comprehensive data to understand the nature of multiallelic violations
and why different enforcement modes succeed or fail.
"""

import argparse
from cyvcf2 import VCF
from collections import defaultdict

def count_alts_in_genotype(gt):
    """
    Count the number of ALT alleles in a genotype (0, 1, or 2).
    Returns -1 for missing genotypes.
    """
    if gt is None or len(gt) < 2:
        return -1
    a, b = gt[0], gt[1]
    if a < 0 or b < 0:
        return -1  # missing
    # Count ALT alleles (non-zero)
    return (1 if a > 0 else 0) + (1 if b > 0 else 0)

def hap_state(gt):
    """
    Return which haplotype(s) carry ALT for a phased biallelic GT.
    - {0}  -> ALT on left hap (e.g., 1|0)
    - {1}  -> ALT on right hap (e.g., 0|1)
    - {0,1}-> both haps ALT (e.g., 1|1)
    - None -> unphased / missing / not biallelic 0/1
    - set() -> no ALT (0|0)
    - 'missing' -> missing genotype (./. or .|.)
    """
    if gt is None or len(gt) < 3:
        return None
    a, b, phased = gt[0], gt[1], bool(gt[2])
    if not phased:
        return None
    if a < 0 or b < 0:
        return 'missing'
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

def is_rare_variant(ac, an, rare_threshold=0.0001):
    """
    Determine if a variant is rare based on allele frequency.
    Rare defined as AF < 0.01% or AF > 99.99%
    """
    if an == 0:
        return True  # No data, consider rare
    af = ac / an
    return af < rare_threshold or af > (1 - rare_threshold)

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

def load_unphased_data(unphased_vcf_path, target_positions):
    """
    Load unphased genotype data for specific positions.
    Returns dict: {(chrom, pos): {variant_id: [genotypes]}}
    """
    if not unphased_vcf_path:
        return {}
    
    unphased_data = defaultdict(dict)
    unphased_vcf = VCF(unphased_vcf_path)
    
    for variant in unphased_vcf:
        pos_key = (variant.CHROM, variant.POS)
        if pos_key in target_positions:
            variant_id = f"{variant.CHROM}_{variant.POS}_{variant.REF}_{variant.ALT[0] if variant.ALT else 'N'}"
            unphased_data[pos_key][variant_id] = variant.genotypes
    
    unphased_vcf.close()
    return unphased_data

def process_group(records, sample_names, include_filtered=False, outfh=None, stats=None, unphased_data=None):
    """
    Enhanced processing of multiallelic position groups.
    
    records: list of cyvcf2.Variant at the same CHROM+POS
    Collects comprehensive metrics about violations including:
    - ALT counts per haplotype (total and common-only)
    - Number of variants at position
    - Missing genotype counts
    - Rare variant involvement
    """
    if not records:
        return 0
    
    # Map: sample_index -> {0: set(ALTs), 1: set(ALTs)}
    conflicts = 0
    hap_alts = [ {0:set(), 1:set()} for _ in sample_names ]
    hap_common_alts = [ {0:set(), 1:set()} for _ in sample_names ]
    missing_counts = [0] * len(sample_names)
    
    # For imputation analysis
    pos_key = (records[0].CHROM, records[0].POS)
    original_missing_counts = [0] * len(sample_names)
    imputed_alts_counts = [0] * len(sample_names)

    # Preload per-record ALT strings and rare variant info
    rec_data = []
    variants_at_pos = 0
    
    for v in records:
        if (not include_filtered) and v.FILTER not in (None, 'PASS', '.'):
            rec_data.append(None)  # treated as skipped
            continue
            
        alt = (v.ALT or [None])[0]
        if not alt:
            rec_data.append(None)
            continue
            
        variants_at_pos += 1
        
        # Extract AC and AN from INFO field
        ac = v.INFO.get('AC', 0)
        an = v.INFO.get('AN', 0)
        
        # Handle AC as list (for multiallelic sites split into biallelic)
        if isinstance(ac, (list, tuple)):
            ac = ac[0] if len(ac) > 0 else 0
        if isinstance(an, (list, tuple)):
            an = an[0] if len(an) > 0 else 0
            
        is_rare = is_rare_variant(ac, an)
        
        rec_data.append({
            'alt': alt,
            'is_rare': is_rare,
            'ac': ac,
            'an': an,
            'variant_id': f"{v.CHROM}_{v.POS}_{v.REF}_{alt}"
        })

    # Process genotypes for each record
    for r_idx, v in enumerate(records):
        rec_info = rec_data[r_idx]
        if not rec_info:  # skip filtered/invalid
            continue
            
        alt = rec_info['alt']
        is_rare = rec_info['is_rare']
        
        gts = v.genotypes  # list of [a,b,phased,(GQ...)]
        
        # Get unphased data for comparison if available
        unphased_gts = None
        if unphased_data and pos_key in unphased_data:
            variant_id = rec_info['variant_id']
            if variant_id in unphased_data[pos_key]:
                unphased_gts = unphased_data[pos_key][variant_id]
        
        for i, gt in enumerate(gts):
            hs = hap_state(gt)
            
            # Compare with unphased data if available
            if unphased_gts and i < len(unphased_gts):
                unphased_gt = unphased_gts[i]
                unphased_alt_count = count_alts_in_genotype(unphased_gt)
                phased_alt_count = count_alts_in_genotype(gt)
                
                if unphased_alt_count == -1:
                    # Was missing in unphased
                    original_missing_counts[i] += 1
                    if phased_alt_count > 0:
                        # ALTs were imputed
                        imputed_alts_counts[i] += phased_alt_count
                elif phased_alt_count > unphased_alt_count:
                    # More ALTs in phased than unphased (imputation increased)
                    imputed_alts_counts[i] += (phased_alt_count - unphased_alt_count)
            
            if hs == 'missing':
                missing_counts[i] += 1
                continue
            elif hs is None or len(hs) == 0:
                continue
                
            if hs == {0,1}:
                # ALT on both haps => add to both
                for h in (0,1):
                    hap_alts[i][h].add(alt)
                    if not is_rare:
                        hap_common_alts[i][h].add(alt)
            else:
                h = next(iter(hs))
                hap_alts[i][h].add(alt)
                if not is_rare:
                    hap_common_alts[i][h].add(alt)

    # After collecting, check each sample for conflicts and generate detailed metrics
    for i, per_hap in enumerate(hap_alts):
        hap0_alts = per_hap[0]
        hap1_alts = per_hap[1]
        hap0_common_alts = hap_common_alts[i][0]
        hap1_common_alts = hap_common_alts[i][1]
        
        # Check for violations on either haplotype
        violation_found = False
        for h in (0,1):
            if len(per_hap[h]) > 1:
                violation_found = True
                break
        
        if violation_found:
            # Count rare ALTs involved in violation
            all_alts_in_violation = hap0_alts | hap1_alts
            rare_alt_count = 0
            
            for r_idx, rec_info in enumerate(rec_data):
                if rec_info and rec_info['alt'] in all_alts_in_violation and rec_info['is_rare']:
                    rare_alt_count += 1
            
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
                    
                    # Metrics
                    hap0_alt_count = len(hap0_alts)
                    hap1_alt_count = len(hap1_alts)
                    hap0_common_alt_count = len(hap0_common_alts)
                    hap1_common_alt_count = len(hap1_common_alts)
                    missing_count = missing_counts[i]
                    original_missing = original_missing_counts[i]
                    imputed_alts = imputed_alts_counts[i]
                    
                    outfh.write(f"{chrom}\t{pos}\t{sname}\thap{h}\t" +
                                ",".join(alts_sorted) + 
                                f"\t{violation_type}\t{hap0_alt_count}\t{hap1_alt_count}\t" +
                                f"{hap0_common_alt_count}\t{hap1_common_alt_count}\t" +
                                f"{variants_at_pos}\t{missing_count}\t{rare_alt_count}\t" +
                                f"{original_missing}\t{imputed_alts}\n")
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
    ap = argparse.ArgumentParser(description="Comprehensive analysis of multiallelic violations with detailed metrics including ALT counts, missing data, rare variant involvement, and imputation analysis.")
    ap.add_argument("vcf", help="Input phased, split-biallelic VCF/BCF (.vcf.gz/.bcf), indexed")
    ap.add_argument("-o", "--out", default="-", help="Output TSV (default: stdout)")
    ap.add_argument("--unphased", help="Optional unphased VCF/BCF for imputation analysis")
    ap.add_argument("--include-filtered", action="store_true",
                    help="Include records with non-PASS FILTER")
    ap.add_argument("--rare-threshold", type=float, default=0.0001,
                    help="Allele frequency threshold for rare variants (default: 0.0001 = 0.01%%)")
    args = ap.parse_args()

    v = VCF(args.vcf)
    samples = v.samples

    # First pass: collect all multiallelic positions for unphased data loading
    multiallelic_positions = set()
    if args.unphased:
        for group in iterate_positions(v):
            if len(group) >= 2:
                multiallelic_positions.add((group[0].CHROM, group[0].POS))
        v.close()
        v = VCF(args.vcf)  # Reopen for second pass
    
    # Load unphased data if provided
    unphased_data = load_unphased_data(args.unphased, multiallelic_positions) if args.unphased else None

    outfh = open(args.out, "w") if args.out != "-" else __import__("sys").stdout
    outfh.write("#CHROM\tPOS\tSAMPLE\tHAP\tALTS_AT_SAME_POS\tVIOLATION_TYPE\t" +
                "HAP0_ALT_COUNT\tHAP1_ALT_COUNT\tHAP0_COMMON_ALT_COUNT\tHAP1_COMMON_ALT_COUNT\t" +
                "VARIANTS_AT_POS\tMISSING_COUNT\tNUMBER_RARE_ALTS\t" +
                "ORIGINAL_MISSING\tIMPUTED_ALTS\n")

    total_conflicts = 0
    violation_stats = {}
    
    for group in iterate_positions(v):
        if len(group) < 2:
            continue  # only meaningful if multiple records at same position
        total_conflicts += process_group(group, samples, args.include_filtered, outfh, violation_stats, unphased_data)

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
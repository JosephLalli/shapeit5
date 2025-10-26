#!/usr/bin/env python3
"""
Create a simple test VCF with split multiallelic violations to test SHAPEIT5 one-allele enforcement.

This creates a minimal VCF with:
- A few samples
- Split multiallelic sites at the same position
- Genotypes that will violate the one-allele constraint when phased independently
"""

import sys

def create_multiallelic_test_vcf():
    """Create test VCF with split multiallelic sites that will cause violations"""
    
    # Create enough samples for SHAPEIT5 (minimum 50)
    n_samples = 60
    sample_names = [f"SAMPLE{i+1}" for i in range(n_samples)]
    
    header = f"""##fileformat=VCFv4.2
##contig=<ID=1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=1,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	{chr(9).join(sample_names)}"""

    # Create genotypes for 60 samples
    # First few samples have violations, rest are mostly REF or random
    def make_genotypes(g1, g2, g3, g4):
        """Make genotypes: first 4 samples get specified genotypes, rest get REF or random"""
        import random
        random.seed(42)  # Deterministic
        gts = [g1, g2, g3, g4]
        
        # Add more samples with mostly REF genotypes and some random variation
        for i in range(4, n_samples):
            if i % 5 == 0:  # Every 5th sample gets ALT
                gts.append(random.choice(["1|0", "0|1"]))
            else:
                gts.append("0|0")
        return "\t".join(gts)
    
    # Create split multiallelic sites at position 1000
    # These genotypes will be unphased - SHAPEIT5 will phase them independently
    # which should create violations when split multiallelic sites get phased separately
    
    def unphase_genotype(gt):
        """Convert phased genotype to unphased: 0|1 or 1|0 -> 0/1, others unchanged"""
        if gt in ['0|1', '1|0']:
            return '0/1'
        elif gt in ['0|0', '1|1']:
            return gt.replace('|', '/')
        else:
            return gt  # Missing data unchanged
    
    def generate_pl_values(gt):
        """Generate realistic PL values for a genotype"""
        import random
        if gt == '0/0':
            # Homozygous reference: high confidence in 0/0, low in others
            return f"{random.randint(0,5)},{random.randint(30,60)},{random.randint(60,120)}"
        elif gt == '0/1':
            # Heterozygous: high confidence in 0/1, lower in homozygous states
            return f"{random.randint(20,50)},{random.randint(0,5)},{random.randint(20,50)}"
        elif gt == '1/1':
            # Homozygous alt: high confidence in 1/1, low in others
            return f"{random.randint(60,120)},{random.randint(30,60)},{random.randint(0,5)}"
        else:
            # Missing or unusual genotype
            return "0,0,0"
    
    def make_unphased_genotypes(g1, g2, g3, g4):
        """Make unphased genotypes with PL values from phased ones"""
        import random
        random.seed(42)  # Deterministic
        gts_with_pl = []
        
        # Process first 4 samples
        for gt in [unphase_genotype(g1), unphase_genotype(g2), unphase_genotype(g3), unphase_genotype(g4)]:
            pl = generate_pl_values(gt)
            gts_with_pl.append(f"{gt}:{pl}")
        
        # Add more samples with mostly REF genotypes and some random variation  
        for i in range(4, n_samples):
            if i % 5 == 0:  # Every 5th sample gets ALT
                gt = '0/1'  # Unphased heterozygote
            else:
                gt = '0/0'
            pl = generate_pl_values(gt)
            gts_with_pl.append(f"{gt}:{pl}")
        return "\t".join(gts_with_pl)
    
    # Create unphased genotypes - SHAPEIT5 will phase these independently
    # When it phases the split multiallelic sites independently, it should create violations
    variants = [
        f"1\t900\t.\tG\tA\t60\tPASS\tAC=12;AN={2*n_samples}\tGT:PL\t{make_unphased_genotypes('1|0', '0|1', '1|0', '0|1')}",
        f"1\t1000\t.\tA\tG\t60\tPASS\tAC=15;AN={2*n_samples}\tGT:PL\t{make_unphased_genotypes('1|0', '0|1', '0|1', '0|0')}",  # Potential violations after phasing
        f"1\t1000\t.\tA\tT\t60\tPASS\tAC=18;AN={2*n_samples}\tGT:PL\t{make_unphased_genotypes('0|1', '0|1', '0|1', '0|0')}",  # Potential violations after phasing
        f"1\t1000\t.\tA\tC\t60\tPASS\tAC=10;AN={2*n_samples}\tGT:PL\t{make_unphased_genotypes('0|0', '0|0', '0|1', '0|0')}",  # Potential violations after phasing
        f"1\t1100\t.\tT\tC\t60\tPASS\tAC=20;AN={2*n_samples}\tGT:PL\t{make_unphased_genotypes('0|1', '1|0', '0|1', '1|0')}"
    ]
    
    print(header)
    for variant in variants:
        print(variant)

if __name__ == "__main__":
    create_multiallelic_test_vcf()
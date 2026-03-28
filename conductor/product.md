# Initial Concept

SHAPEIT5 estimates haplotypes in large datasets, with a special focus on rare variants.

## Product Guide

### Target Users
- Bioinformatics researchers and geneticists.

### Core Features & Capabilities
- **State-of-the-Art Phasing:** Provides industry-leading accuracy for both common variants (SNP arrays) and rare variants (scaffold-based phasing).
- **Native Multiallelic Support:** Extends the core phasing algorithms to natively handle multiallelic variants, removing the need for "split" records and the legacy "anchor-sibling" workaround.
- **High Performance:** Aims to maintain high efficiency (within 10-20% of the biallelic baseline) while processing large-scale datasets.
- **Supersite Refinement:** Builds upon previous "supersite" logic to implement true native reading, writing, and internal representation of multiallelic sites.

### Project Goals
- **Zero Regression:** Paramount importance is placed on maintaining identical phasing accuracy for biallelic sites compared to the previous stable version.
- **Native Data Structures:** Transition from 4-bit packed genotypes (max 15 alleles) to 8-bit (1-byte) genotypes to support up to 255 alternative alleles per site.
- **Architectural Overhaul:**
    - Update Input/Output (XCF/BCF) to read/write multiallelic records as single units.
    - Refactor the HMM to emit probabilities for >2 alleles natively.
    - Remove the "anchor-sibling" dichotomy in favor of a unified site representation.

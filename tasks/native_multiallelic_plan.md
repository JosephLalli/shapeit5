# Native Multiallelic Supersites Plan

## Goal
- Keep the existing biallelic vs supersite split.
- Build supersites directly from multiallelic BCF records (no anchor+siblings).
- Use 8-bit allele codes: 0=REF, 1..n_alts (up to 255).
- Missing rule: any missing allele -> supersite missing (both haps).
- MAF/MAC: ref vs nonref; ALT count is the sum across all ALT alleles.
- Binary haplotype formats are unsupported when multiallelic records are present.

## Data Model
- `SuperSite`:
  - `n_alts` / `n_classes` as `uint16_t` (1..256 classes).
  - `panel_offset`/`panel_span_bytes` for 1 byte per hap code.
  - `rare_code_mask` as 256-bit bitset (e.g., `std::array<uint64_t,4>`).
- Panel codes:
  - `packed_allele_codes` is a flat `std::vector<uint8_t>`; layout `[supersite][hap]`.
  - No 4-bit packing.
- Genotype supersites:
  - `ss_observed_gts` / `ss_phased_gts` store allele codes per supersite.
  - `ss_missing_mask` (2 bits per supersite) replaces missing sentinel.
  - Anchor `Variants` MIS/HET/HOM flags derived from codes + missing mask.
- Haplotype set:
  - `H_supersite_codes` stores per-hap codes for all supersites.

## Input Rules
- Multiallelic records must be read from VCF/BCF.
- If any multiallelic record is present, binary haplotype formats for reference/scaffold are rejected.

## Implementation Steps
1) **Reader changes**
   - `xcf_reader`: accept `n_allele >= 2`, expose `n_allele` + full ALT list; sum AC across all ALTs for AF.
   - `genotype_reader_scaning`: stop dropping multiallelic; track presence; enforce binary-haplotype restriction; update SNP filter to require all ALTs be SNPs.
   - `genotype_reader_reading`: decode GT allele indices into per-hap supersite codes + missing mask.

2) **Supersite construction**
   - One supersite per multiallelic record.
   - Build `packed_allele_codes` from `H_supersite_codes` (1 byte per hap).
   - Compute `rare_code_mask` from code frequencies.

3) **Genotype integration**
   - Update `snapshotSupersiteObservedGts`, `setSupersitePhasedGt`, `projectSupersites` for direct allele codes.
   - Remove sibling logic in segment building and ambiguous indexing.

4) **PBWT + HMM**
   - `class_code()` uses 8-bit supersite codes.
   - Emissions and imputation use missing masks + class codes; biallelic path unchanged.

5) **Output + tests**
   - Write multiallelic GT/HP using stored allele codes.
   - Add multiallelic BCF test coverage.

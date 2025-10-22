# SHAPEIT5 Testing Architecture — 2025-10-22

## Headerless VCF.gz + MD5 Validation System

The SHAPEIT5 test suite uses a deterministic validation system designed to eliminate non-deterministic elements that cause false test failures.

### Problem Solved
Traditional VCF comparison suffered from:
- Timestamp headers (`##fileDate=21/10/2025 - 23:05:17`)
- Tool version metadata (`##bcftools_viewVersion=1.22+htslib-1.22.1`)
- Temporary file paths (`##bcftools_viewCommand=view /tmp/tmp.xyz/...`)
- Varying contig formats between tool versions (`##contig=<ID=1>` vs `##contig=<ID=1,length=10000000>`)

### Solution: Headerless + MD5
1. **Strip all headers**: Use `bcftools view -H` to extract only variant data lines
2. **Compress deterministically**: Store as `.vcf.gz` files
3. **MD5 validation**: Compare MD5 checksums instead of line-by-line diff
4. **Ground truth files**: Generated using published SHAPEIT5 binaries for consistency

### File Structure
```
test/scripts/expected/
├── phase.array.family.vcf.gz     # Headerless compressed variants
├── phase.array.family.md5        # MD5 checksum for validation
├── phase.array.reference.vcf.gz  
├── phase.array.reference.md5
└── ... (8 test scenarios total)
```

### Validation Functions (`test/scripts/lib/test_utils.sh`)
- `extract_variants()`: Converts BCF to headerless VCF.gz using `bcftools view -H`
- `assert_same_md5()`: Compares MD5 checksums for validation
- `generate_headerless_vcf_gz()`: Creates ground truth files from published binaries

### Benefits
- **100% deterministic**: Immune to timestamp/metadata variations
- **Faster validation**: MD5 comparison vs multi-megabyte diff operations
- **Smaller files**: ~70% size reduction through compression
- **Version-agnostic**: Works across bcftools/SHAPEIT5 versions
- **Clear semantics**: MD5 mismatch = functional difference in phasing results

### Test Coverage
- Array-based phasing (5 scenarios): reference, unrelated, family, scaffold, haploid
- WGS-based phasing (3 scenarios): unrelated, family, haploid
- Single-threaded execution for reproducibility (`--thread 1`) and fixed seed (`--seed 15052011`)
- Region-restricted datasets for speed; expected files are headerless `.vcf.gz` with `.md5`

### Usage Guidelines for AI Agents
When modifying SHAPEIT5 tests:
1. **Always use MD5 validation functions** from `test_utils.sh`
2. **Generate expected files using published SHAPEIT5 binaries** for consistency
3. **Verify all tests pass** before committing changes
4. **Never compare raw VCF files with headers** due to non-determinism
5. **Use single-threaded execution** (`--thread 1`) for reproducible results

### Implementation Status
- ✅ Single-threaded execution confirmed
- ✅ MD5-based validation in place for all scenarios
- ✅ All 8 test scenarios pass functional validation

## WGS PAR2 Fixtures (1KGP T2T)

The WGS tests use CHM13 T2T 1KGP PAR2 inputs on chrX, with known multiallelic violations to exercise one-allele enforcement.

### Region and Map
- Region: `chrX:153929053-154248138` (PAR2)
- Map: `test/info/par2.gmap.gz`

### Inputs (BCF)
- Family: `test/wgs/target.family.1kgp_t2t.par2.bcf` (≈1803 samples)
- Unrelated: `test/wgs/target.unrelated.1kgp_t2t.par2.bcf` (≈1399 samples)
- Haploid cohort: `test/wgs/target.haploid.1kgp_t2t.par2.bcf` (≈1599 male samples)

These BCFs carry only `FORMAT:GT` (no `FORMAT:PL`). This is sufficient for phasing workflows used in tests.

### Sample Lists and Pedigree
- Family union from PED (Son, Father, Mother; drop NA/-1): `test/wgs/samples.family.1kgp_t2t.par2.txt`
- All samples in master BCF: `test/wgs/samples.all.1kgp_t2t.par2.txt`
- Unrelated complement: `test/wgs/samples.unrelated.1kgp_t2t.par2.txt`
- Haploid sample list (first column only): `test/wgs/samples.haploid.1kgp_t2t.par2.txt`
- Pedigree TSV used by scripts: `test/info/1kgp_t2t.par2.ped`

### Haploid Handling
- The haploid BCF remains diploid in representation (GT); haploid samples are declared via `--haploids info/1kgp_t2t.par2.haploid.txt` so heterozygous GTs are reset to missing within the pipeline.
- We do not require `bcftools +fixploidy` in the final setup; if ever needed, `test/wgs/par2.ploidy.txt` and `test/wgs/samples.haploid.sex.txt` illustrate ploidy control.

### Scripts
- `test/scripts/phase.wgs.family.sh`: uses family BCF + PED, region and map above
- `test/scripts/phase.wgs.unrelated.sh`: uses unrelated BCF
- `test/scripts/phase.wgs.haploid.sh`: uses haploid BCF + haploid list

All three use `--map test/info/par2.gmap.gz`, `--region chrX:153929053-154248138`, `--seed 15052011`, and `--thread 1`.

### Expectations
- Ground-truth outputs live in `test/scripts/expected/` as headerless `.vcf.gz` with corresponding `.md5` files:
  - `phase.wgs.family.{vcf.gz,md5}`
  - `phase.wgs.unrelated.{vcf.gz,md5}`
  - `phase.wgs.haploid.{vcf.gz,md5}`

### Violation Spot-Checks
- The repository includes `test/find_one_alt_per_hap_violation.py` to confirm split multiallelic positions with same-haplotype ALT duplications are present in the inputs. Running it on the three BCFs reports multiple violations in each subset, ensuring the constraint logic is exercised.

## Array Fixtures

Array tests remain unchanged functionally but now also pin `--seed 15052011` and `--thread 1`. Expected outputs are maintained as headerless `.vcf.gz` + `.md5` pairs in `test/scripts/expected/`.

## Determinism & Reproducibility

- Always run tests with `--thread 1` and a fixed `--seed` to avoid haplotype-orientation drift.
- Validate with `assert_same_md5` (not raw VCF diffs) to avoid spurious header/metadata mismatches.
- Expected files should be regenerated using known-good binaries and then committed.

This architecture ensures reliable regression testing and eliminates false failures from metadata differences, providing a solid foundation for SHAPEIT5 development.

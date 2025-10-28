Project: SHAPEIT5 — phase_common super-site integration

Objectives
- Implement optional super-site phasing for multi-allelic loci (e.g., STRs) behind `--enable-supersites`.
- Maintain identical default behavior when the option is not set.

Plan
1) Data Model
- Add `phaser` members for super-site indexing and membership:
  - `std::vector<int> locus_to_super_idx` — maps variant locus to super-site index or -1.
  - `std::vector<int> super_site_var_index` — flattened storage of per–super-site variant indices.
- Extend `SuperSite` to reference its member variants:
  - `uint32_t var_start;` and `uint16_t var_count;` (span in `super_site_var_index`).
- Keep `packed_allele_codes` (conditioning panel codes, 2 haplotypes per byte with 4-bit codes).

2) Super-site Builder
- Implement and wire `buildSuperSites(variant_map&, conditioning_set&, ...)` to fill:
  - `super_sites`, `locus_to_super_idx`, `packed_allele_codes`, `super_site_var_index`.
- Group split biallelic records by (chr,bp). For >15 ALTs, split deterministically into multiple super-sites.
- Encode panel haplotype codes using `H.Hhap` and pack two 4-bit codes per byte.

3) Accessors
- Add helper to derive a sample’s super-site allele code on-demand from the constituent variants:
  - `getSampleSuperSiteAlleleCode(const genotype*, const SuperSite&, const std::vector<int>&, int hap)`.
- Continue using existing unpackers (`unpackSuperSiteCode*`) and AVX2 emission precompute.

4) HMM Integration
- Pass super-site state to `haplotype_segment_single/double` via optional pointers:
  - `super_sites`, `locus_to_super_idx`, `packed_allele_codes.data()`, `super_site_var_index`.
- At each locus inside forward/backward:
  - `int ss_idx = (*locus_to_super_idx)[curr_abs_locus];` gate super-site path when `ss_idx >= 0`, else run current biallelic path.
- Provide both float/double paths:
  - Start with scalar or lightly vectorized float path for correctness, keep AVX2 path in double; optimize float with AVX2 later.

5) Driver Lifecycle
- Build supersites once per run (first iteration): call builder when `--enable-supersites && current_iteration == 0`.
- When disabled (default), pass `nullptr` to segments and keep existing behavior.

6) Validation
- Unit tests:
  - Builder: synthetic groups → verify `locus_to_super_idx`, `super_site_var_index`, packing correctness.
  - Accessor: sample-code derivation matches expected from input variants.
  - Emissions: AVX2 vs scalar precompute agree within tolerance.
- Smoke tests:
  - Biallelic-only input with and without `--enable-supersites` yields matching outputs (stochastic variation aside).

Notes
- P0 fixed: a single canonical `aligned_vector32` alias is used from `utils/otools.h`.
- Keep changes isolated to phase_common; default runs remain unaffected unless `--enable-supersites` is set.


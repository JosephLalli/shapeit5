# SHAPEIT5 — Multiallelic Site Constraint Enforcement (Handover)

This document explains, at a practical engineering level, what SHAPEIT5 does, why multiallelic sites are tricky, the three remedies we implemented (in rising complexity), how they’re wired into the codebase, and the main challenges we hit. It is written for a strong C++ developer without a genetics background.

## 1) What SHAPEIT5 currently does

- Phases diploid genotypes into haplotypes using:
  - PBWT (Positional Burrows–Wheeler Transform) to select a small panel of donor haplotypes per window.
  - A Li–Stephens HMM to sample each sample’s haplotypes across windows (burn‑in, pruning, main iterations).
- Two pipelines:
  - Common variants: `phase_common` (array/WGS)
  - Rare variants on a scaffold: `phase_rare`
- Variant representation is biallelic (REF vs ALT). Real multiallelic sites are split into several biallelic rows with identical (chr, pos).

## 2) Why multiallelic handling is a problem

**The biological constraint**: At any genomic position, a haplotype can carry at most one allele (REF or one specific ALT). 

**The input format**: SHAPEIT5 takes split biallelic sites as input. A multiallelic site like:
```
chr1  1000  .  A  G,T,C    0|2   1|3  2|2    ./.
```
gets split into separate biallelic rows:
```
chr1  1000  .  A  G  0|0    1|0    0|0    ./.
chr1  1000  .  A  T  0|1    0|0    1|1    ./.
chr1  1000  .  A  C  0|0    0|1    0|0    ./.
```

**The problem**: SHAPEIT5's HMM treats these rows independently during phasing and imputation. This can produce biologically impossible results where a single haplotype has multiple ALT alleles at the same position:
```
chr1  1000  .  A  G  0|0    0|1    0|0    1|0
chr1  1000  .  A  T  0|1    0|0    1|1    1|0  
chr1  1000  .  A  C  0|0    0|1    0|0    1|0
```
Here, sample 2 has both G and C on haplotype 1, and sample 4 has G, T, and C all on haplotype 0 — both biologically impossible. This misleads downstream tools (imputation, burden tests, STR analysis).

## 3) Three solutions (in order of increasing complexity)

All three run after an HMM sampling step — they are small, local consistency fixes that do not alter the core HMM.

- Transition‑only post‑hoc flipping (simplest)
  - Use nearest phased heterozygous anchors p (left) and n (right) within the same phase set. For a violation (e.g., two ALT rows on hap0), evaluate the two ways to flip one row to hap1 using Li–Stephens stay/switch probabilities on [p→site] and [site→n]; choose the higher likelihood.
  - Pros: O(1) per violation, minimally invasive.
  - Limits: Does not consider per‑site emissions; does not use the actual donor identities; with >2 rows or missing anchors, may not fully resolve.

- Micro re‑decode (donor‑agnostic)
  - Enumerate all valid compound assignments for the K split rows (≤ (K+1)^2, pruned to ≤1 ALT per haplotype). Score each assignment with the run’s transition model (and a simple emission proxy or GL/PL if present). Apply the ML assignment and rewrite bits.
  - Pros: Handles arbitrary K rows; can fix both allele identity and haplotype assignment together.
  - Limits: Still donor‑agnostic; emissions in `phase_common` are approximate without GL/PL.

- Micro re‑decode with frozen donors (micro‑donor, most integrated)
  - Use the donor haplotypes (PBWT Kstates) selected for the sample/window. For each candidate, compute a small donor‑weighted chain score across p→site→n to reward candidates supported by donors.
  - Pros: Better ties resolution and consistency with the model actually used in the pass.
  - Limits: Current code implements a donor‑weighted chain (emissions + stay/switch) rather than a full donor‑conditioned DP; requires window/donor plumbing.

## 4) Current status in `phase_common` and `phase_rare`

### phase_common

- Feature flags:
  - `--enforce-oneallele` (enable) and `--oneallele-stats <file>`
  - `--oneallele-mode {transition|micro|micro-donor}` (select algorithm)

- Where it runs:
  - Per sample, inside worker loop after each sampling step: `phase_common/src/phaser/phaser_algorithm.cpp`
  - Final “belt‑and‑suspenders” pass after `G.solve()`: `phase_common/src/phaser/phaser_finalise.cpp`

- Key modules and functions:
  - `phase_common/src/modules/multiallelic_position_map.{h,cpp}`: build groups of split rows by (chr, pos).
  - `phase_common/src/modules/transition_scorer.{h,cpp}`: stay/switch scoring helpers.
  - `phase_common/src/modules/oneallele_enforcer.{h,cpp}`:
    - `set_enabled()`, `set_mode()`, `stats()` (API)
    - `enforce_sample(...)`: per‑sample entry point in worker; dispatches to the selected mode per position group.
    - `enforce(...)`: batch fallback (used at finalize step)
    - `enforce_group_transition(...)`: transition‑only fix at a position
    - `enforce_group_micro(...)`: micro/micro‑donor enumeration across arbitrary K rows
    - `evaluate_candidate_micro(...)`: scores a single compound assignment
    - `compute_chain_score_with_donors(...)`: donor‑weighted chain score
    - Helpers: `find_left_neighbor(...)`, `find_right_neighbor(...)`, `flip_alt_to_other_hap(...)`, `rebalance_excess_alts(...)`

### phase_rare

- Feature flags:
  - `--enforce-oneallele-rare`, `--oneallele-rare-stats`
- Where it runs:
  - After merging rare data back into the scaffold (step 3.5): `phase_rare/src/phaser/phaser_algorithm.cpp`
- Approach:
  - Rare pipeline stores per‑row phase probabilities (PP). It resolves violations by flipping the lowest‑confidence ALT rows on the offending haplotype until each hap has ≤1 ALT. This is conceptually similar to a simplified micro decision but tailored to the sparse data structure.

## 5) Challenges encountered

- Access to donors/windows outside the worker: we solved this by running enforcement per sample inside the worker loop.
- Anchor detection: respecting phase‑set/segment boundaries is required; anchors outside the segment are ignored.
- Emissions in `phase_common`: without GL/PL we rely on a simple error model for micro scoring; accuracy improves if GL/PL are present.
- Determinism/ties: we use deterministic tie‑breaking to keep runs reproducible.
- Performance: transition‑only is essentially free; micro/micro‑donor remain cheap with small K and m (typically ≤8).

## Code map (quick reference)

- Common pipeline
  - `phase_common/src/phaser/phaser_algorithm.cpp`: per‑sample enforcement call sites
  - `phase_common/src/phaser/phaser_finalise.cpp`: final pass and stats output
  - `phase_common/src/phaser/phaser_parameters.cpp`: CLI flags and mode selection
  - `phase_common/src/modules/*`: implementation of grouping, scoring and enforcement

- Rare pipeline
  - `phase_rare/src/phaser/phaser_algorithm.cpp`: rare enforcement (step 3.5)
  - `phase_rare/src/phaser/phaser_parameters.cpp`: rare flags

## Running and validating

- Typical command (common):
  - `phase_common/bin/phase_common --input X.bcf --map GMAP.gz --region 1 --output out.bcf --enforce-oneallele --oneallele-mode micro --seed 42`
- Checker: `tools/check_oneallele out.bcf` reports per‑haplotype violation counts per position.
- See `test/scripts` for integration exercises, including a multiallelic‑focused script that runs all three modes.

## Testing System for Multiallelic Enforcement

**Critical testing approach**: To validate the multiallelic enforcement algorithms:

1. **Create unphased input data** with potential multiallelic violations:
   - Use split biallelic sites at the same genomic position (e.g., chr1:1000 A→G, chr1:1000 A→T, chr1:1000 A→C)  
   - Input genotypes must be unphased (0/1, not 0|1) so SHAPEIT5 phases them independently
   - Include samples with multiple heterozygous calls at the same position

2. **Let SHAPEIT5 phase independently** without enforcement:
   - Run `phase_common --seed 42` (no --threads) to ensure deterministic results
   - Independent phasing of split multiallelic sites creates violations where single haplotypes get multiple ALT alleles at same position

3. **Test enforcement algorithms** on the same input:
   - Run with `--enforce-oneallele --oneallele-mode {transition|micro|micro-donor} --seed 42`
   - Compare outputs to verify violations are detected and corrected
   - Check final statistics: `violations=N / flips=M` in console output

**Deterministic execution requirements**:
- Always use `--seed 42` (or any fixed seed) for reproducible results
- Never specify `--threads` flag as it causes non-deterministic behavior even with `--threads 1`
- Use identical input data and parameters for before/after comparisons

This testing methodology ensures that enforcement algorithms are validated against real violations created by the core SHAPEIT5 phasing process, rather than artificial pre-phased data.


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

## Handover: Where We Left Off and Next Steps

This section captures the precise development state in the codebase and the next engineering steps to complete multiallelic one-allele enforcement as described above.

### Current Implementation Snapshot (October 2025)

- Common pipeline (phase_common)
  - ✅ **Fully Implemented and Working**
    - Mode selection: `--oneallele-mode {transition|micro|micro-donor}` flag parsing with validation in `phase_common/src/phaser/phaser_parameters.cpp`
    - **Transition mode**: Complete transition-only enforcement with left/right heterozygous anchors per phase segment using `transition_scorer.{h,cpp}`
      - Statistics reporting fixed (cumulative counting bug resolved)
      - Handles violations where all ALT alleles are concentrated on one haplotype
      - Li-Stephens scoring with anchor detection working correctly
      - Produces biologically valid results (both haplotypes ≤1 ALT)
    - Final enforcement pass: Belt-and-suspenders cleanup after `G.solve()` in finalize step
    - Epoch reporting: Per-iteration stats with timing: `"Multiallelic correction (mode) [violations=N / flipped=N] (X.XXs)"`
    - Stats infrastructure: `OneAlleleStats` for cumulative stats, `OneAlleleEpochStats` for per-epoch tracking (bugs fixed)
    - Core execution: Enforcement runs after each MCMC iteration in `phase_common/src/phaser/phaser_algorithm.cpp` with proper timing measurement
    - CLI integration: All flags (`--enforce-oneallele`, `--oneallele-mode`, `--oneallele-stats`) parsed and validated
  
  - ✅ **Mostly Implemented (Functional)**
    - **Micro mode**: Donor-agnostic enumeration algorithm in `oneallele_enforcer.{h,cpp}`:
      - `MicroCandidate` structure for assignment enumeration  
      - `enumerate_micro_candidates()` generates all valid assignments (≤ (K+1)² combinations, ≤1 ALT per hap)
      - `evaluate_candidate_micro()` scores with Li-Stephens transition model + simple emission model
      - `apply_micro_candidate()` applies maximum likelihood assignment
      - **Status**: Implementation complete, needs comprehensive testing

  - ⚠️ **Partially Implemented (Incomplete)**
    - **Micro-donor mode**: Intended to use PBWT donor-weighted scoring:
      - `enforce_sample()` method for per-sample enforcement with donor context access ✅
      - `enforce_group_micro()` with donor context parameter ✅
      - **Missing**: `evaluate_candidate_micro_donor()` currently stubs out to basic micro mode (line 905-907)
      - **Missing**: `compute_chain_score_with_donors()` is incomplete stub
      - **Current behavior**: Falls back to micro mode (donor-agnostic scoring)
      - Integration with worker threads for donor Kstates access during sampling ✅

- Rare pipeline (phase_rare)
  - Implemented
    - Enforcement after merging rare variants onto scaffold (step 3.5) with PP-based flipping
    - Flags `--enforce-oneallele-rare` and `--oneallele-rare-stats` with reporting
    - Enhanced violation reporting: tracks extreme violations (>2 ALT alleles at same position in same sample)
    - Stats include positions_checked, sample_violations_found, extreme_violations_found, flips_applied
  - Current approach sufficient
    - PP-based approach works well for sparse rare variant data structure

- Testing and validation 
  - Unit tests: `test/unit/test_oneallele_enforcer.cpp` covers transition scoring and violation resolution
  - Integration tests: WGS PAR2 and array scripts in `test/scripts/` with deterministic MD5 validation
  - Verification tools: `tools/check_oneallele` for output validation  
  - Test coverage: All three modes (transition/micro/micro-donor) can be tested
  - Build verification: Compiles successfully with `-j 24`, all existing tests pass
  - Test script fixes: Fixed line continuation in `test/scripts/phase.oneallele.array.family.sh`

### Next Steps (Remaining Work - Priority Ordered)

**IMMEDIATE PRIORITIES (High Impact)**

1) **Complete micro-donor mode implementation** - Core Feature Completion
   - **Goal**: Implement full PBWT donor-weighted scoring for micro-donor mode
   - **Tasks**: 
     - Implement `evaluate_candidate_micro_donor()` with actual donor weighting (currently stub at line 905-907)
     - Complete `compute_chain_score_with_donors()` to use PBWT Kstates for scoring
     - Integrate donor haplotype information into Li-Stephens emission probability calculations
     - Test donor-weighted vs donor-agnostic scoring differences
   - **Files**: `phase_common/src/modules/oneallele_enforcer.cpp` (lines 899-907, 913+)
   - **Complexity**: Medium-High (requires understanding PBWT integration)
   - **Current Status**: Infrastructure exists, core algorithm missing

2) **Comprehensive testing and validation** - Quality Assurance
   - **Goal**: Validate all three enforcement modes work correctly and identify remaining edge cases
   - **Tasks**: 
     - Test micro mode comprehensively (implementation appears complete but untested)
     - Create test datasets with complex multiallelic violations (>2 ALTs, mixed patterns)
     - Compare mode effectiveness: transition vs micro vs micro-donor on same violations
     - Performance testing with large datasets and varying violation densities
     - Validate final statistics accuracy across all modes
   - **Files**: `test/scripts/`, custom test data generation, validation scripts
   - **Complexity**: Medium (requires careful test design and analysis)

**MAINTENANCE TASKS (Low Risk)**

3) **Enhanced debugging and statistics** - User Experience
   - **Goal**: Improve debugging capabilities and statistics reporting
   - **Tasks**:
     - Fix debug target infrastructure (currently not producing output)
     - Add per-mode violation type counters to distinguish what each mode handles
     - Enhanced finalization reporting with mode-specific statistics breakdown
     - Optionally expose `--oneallele-min-cm` for transition scoring distance control
   - **Files**: `oneallele_enforcer.{h,cpp}`, `phaser_parameters.cpp`, `phaser_finalise.cpp`
   - **Complexity**: Low-Medium (mostly bookkeeping and debugging infrastructure)

**OPTIONAL ENHANCEMENTS (Future Work)**

4) **Performance optimization** - Algorithmic Efficiency
   - **Goal**: Optimize enforcement algorithms for large datasets with many violations
   - **Tasks**: 
     - Profile micro enumeration performance with large K values
     - Implement early termination heuristics for candidate evaluation
     - Memory optimization for large multiallelic groups
   - **Files**: `oneallele_enforcer.{h,cpp}`
   - **Priority**: Low (current performance adequate for typical use cases)

## Summary for Next Developer

**Current Status**: The multiallelic constraint enforcement feature is **partially implemented** with transition mode fully working and micro/micro-donor modes needing completion.

**What Works Now**:
- **Mode selection**: `--oneallele-mode {transition|micro|micro-donor}` with validation ✅
- **Transition mode**: Fast, lightweight enforcement for simple cases (fully functional, tested, statistics validated) ✅
- **Micro mode**: Sophisticated enumeration algorithm handling arbitrary K variants (implemented, needs testing) ⚠️
- **Micro-donor mode**: Advanced donor-weighted scoring using PBWT context (partially implemented, core algorithm missing) ❌
- **Final enforcement pass**: Belt-and-suspenders cleanup before output (implemented) ✅
- **Epoch reporting**: Real-time stats per iteration showing violations and timing ✅
- **Integration**: Seamlessly integrated into existing MCMC loop with proper timing ✅
- **Testing**: All existing tests pass, builds successfully ✅
- **Rare pipeline**: Enhanced with extreme violation reporting (>2 ALT alleles tracking) ✅

**Immediate Next Steps**: 
1. **Complete micro-donor mode implementation** (priority #1)
2. **Test micro mode comprehensively** to validate enumeration algorithm
3. **Compare all three modes** on same violation datasets

**Estimated Completion**: 
- **Micro-donor mode completion**: 2-3 days (implement donor scoring)
- **Comprehensive testing**: 3-5 days (validate all modes)
- **Total remaining work**: ~1 week for full feature completion

**Key Files for Next Developer**:
- `phase_common/src/modules/oneallele_enforcer.{h,cpp}` - Core enforcement algorithms
- `phase_common/src/phaser/phaser_algorithm.cpp` - Integration point for enforcement calls  
- `phase_common/src/phaser/phaser_parameters.cpp` - CLI flag parsing and validation
- `test/scripts/phase.*.sh` - Test scripts for validation

**Specific Implementation Guidance for Micro-Donor Mode**:

The missing pieces are in `oneallele_enforcer.cpp`:

1. **Lines 905-907**: `evaluate_candidate_micro_donor()` currently stubs out:
   ```cpp
   // TODO: Implement full donor-weighted scoring
   // For now, fall back to donor-agnostic scoring
   return evaluate_candidate_micro(candidate, g, V, variant_to_segment, position_group_indices);
   ```
   **Action**: Replace with donor-weighted Li-Stephens emission scoring using the `Kstates` parameter

2. **Line 913+**: `compute_chain_score_with_donors()` is incomplete
   **Action**: Implement PBWT donor haplotype integration for transition scoring across multiallelic sites

**Reference Implementation Patterns**:
- Look at existing PBWT usage in `phase_common/src/objects/genotype/genotype_*` files
- Study Li-Stephens emission calculations in `haplotype_segment_*` classes
- The `Kstates` parameter contains donor haplotype indices - use these to weight scoring

This implementation provides significant value in current state with transition mode fully functional for the most common violation patterns.

## Design Notes for Future Enhancement

### Transition Scoring Improvements

**User-Configurable Effective Population Size (Ne)**:
The current `transScore` function in `phase_common/src/modules/transition_scorer.cpp` hardcodes `Ne = 10000.0`. This should be made user-configurable via command-line parameter to allow fine-tuning for different populations and datasets.

**Zero Distance Handling for Same-Position Variants**:
For multiallelic sites where multiple biallelic variants exist at the same genomic position, the transition score between variants should be zero (no recombination opportunity). The current implementation uses `std::max(e.distance_cM, 1e-8)` which artificially inflates distances. This needs refinement to properly handle same-position scenarios.

**Comprehensive Haplotype Configuration Scoring**:
The current transition scoring should be enhanced to evaluate both possible haplotype configurations for multiallelic sites:
- Configuration A: `0|1 1|0` (haplotype 0 gets first variant, haplotype 1 gets second)
- Configuration B: `1|0 0|1` (haplotype 1 gets first variant, haplotype 0 gets second)

For each configuration, calculate the full transition chain: `upstream_anchor → multiallelic_site → downstream_anchor` for both haplotypes, then select the configuration with the higher overall likelihood. This provides more principled haplotype assignment than current single-variant flipping approaches.

**Li-Stephens Emission Probability Integration** (MICRO and MICRO_DONOR modes):
The current multiallelic enforcement uses simplified emission scoring (heterozygote bonus) in MICRO and MICRO_DONOR modes even though it runs during `phaseWindow()` when full Li-Stephens HMM probabilities are available. The architecture is well-positioned for enhancement:

- **Current timing**: `enforce_sample()` is called during `phaseWindow()` when `haplotype_segment_*` objects exist with computed forward/backward probabilities
- **Missing access**: The HMM segment object is not passed to `enforce_sample()`, preventing access to `prob[]` (forward), `Alpha*` (backward), and emission probability calculations
- **Enhancement opportunity**: Pass the `haplotype_segment_*` object to `enforce_sample()` and use the same Li-Stephens emission model as the `IMPUTE()` function

This would replace the simple `het_bonus = 0.1` scoring in MICRO mode with true Li-Stephens emission probabilities: `P(assignment | donors) = Σ(donor_weight_k × forward_k × backward_k)`, providing more accurate multiallelic resolution using the same probability model as the core phasing algorithm. TRANSITION mode would remain unchanged as it uses only transition probabilities.

**Implementation**: Modify `enforce_sample()` signature to include HMM segment, extract probability arrays during enforcement window, apply Li-Stephens emission scoring in `compute_emission_score()` for MICRO and MICRO_DONOR modes.

**Memory Efficiency**: The HMM segment object can be passed by const reference (`const haplotype_segment_single& segment`) to avoid copying large probability arrays in memory. The arrays (`prob[]`, `Alpha[]`, `probSumK[]`) already exist in the segment object and can be accessed directly at their memory locations. This follows existing SHAPEIT5 patterns (e.g., `const variant_map& V`, `genotype& g`) and avoids memory overhead while providing read-only access to forward/backward probabilities for Li-Stephens emission calculations.

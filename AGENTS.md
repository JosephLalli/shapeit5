# AGENTS.md

Guidance for working on this repository with the Codex CLI agent. Focus areas: determinism, the multiallelic problem and codebase structure, and how to validate changes.

## Determinism & Reproducibility

- Do not pass any thread flag when you need deterministic output.
  - Using `--thread 1` (phase_rare) or `--threads 1` (some other tools/pipelines) still breaks determinism due to internal ordering and orientation effects.
  - Omit the flag entirely to run single‑threaded deterministically.
- Always set a fixed `--seed`.
  - Recommended values used in tests: `--seed 42` or `--seed 15052011` (see test scripts).

## Multiallelic Support via HMM Supersites

- Problem nature
  - At some genomic coordinates, multiple distinct ALT alleles exist across samples. Splitting these into multiple biallelic rows can yield biologically invalid haplotypes if a single haplotype ends up carrying more than one ALT at the same position.
  - Goal: ensure each haplotype carries at most one ALT per position, enforced inside the HMM so illegal configurations are unreachable, not fixed post‑hoc.

- Codebase structure (pre‑addition baseline)
  - PBWT, conditioning sets, and genotype storage are binary (0/1) and highly optimized.
  - Emission kernels (AVX2) select between two per‑hap emission values based on a boolean predicate per conditioning haplotype.
  - Transition model and MCMC schedule are separate and should remain unchanged.

- Supersite approach (summary)
  - Preprocess input: collapse split biallelic rows sharing (chr, pos) into one logical “supersite”. Create a compact, bit‑packed per‑hap code buffer: 0 = REF, 1..M = ALT codes.
  - Keep core math binary and AVX2‑optimized: biallelic sites run the current fast path unchanged.
  - For supersites, process one effective site per coordinate by selecting anchors (see below) and applying a per‑hap supersite mask at the anchor row; siblings are skipped for emission, PBWT, and transition accounting.
  - PBWT: remain binary and ignore supersites entirely.
  - Transitions: unchanged.
  - Output: re‑split supersites back into biallelic lines without assigning >1 ALT per haplotype.

### Data preparation and hygiene
- Donor panel (reference/scaffold): pre‑scrub so that each donor haplotype carries at most one ALT at any supersite. If a donor violates this, drop the offending variant(s) or the site from the panel.
- Target (unphased inputs): split multiallelic sites to biallelic rows, but do not encode more than one ALT per haplotype at a coordinate in observed data; missing is allowed. If any sample has more than two different ALT alleles at the same coordinate in the input, drop that site from the unphased dataset.

## Implementation Plan

1) Supersite detection and construction
- Group variant records by (chr, bp). If ≥2 ALTs, build one supersite.
- Assign stable ALT codes 1..M with a deterministic canonicalization policy (documented tie‑breakers; no data‑dependent randomness).
- For each conditioning haplotype, assign exactly one code (0 or k). Store codes in a tightly bit‑packed buffer with per‑site `bitwidth` and `panel_offset` metadata.
- Maintain mappings from original split records to their supersite for output reconstruction.

2) Per‑site haplotype code accessor
- Provide an accessor that, for a site and K conditioning hap indices, returns contiguous arrays of small codes for vector loads.
- Biallelic sites: keep existing fast bit extraction (0/1). Supersites: unpack just those K codes using `panel_offset` and `bitwidth`.

3) Emission kernel update (AVX2)
- Per‑hap anchors and supersite masks (v1):
  - At each supersite, select two haplotype‑specific anchors (one per chromosome copy) conceptually, but run a single physical anchor row; construct a per‑hap supersite mask that allows REF or hap‑specific anchor ALT and suppresses other ALTs (multiply by epsilon) at this site. Then run the existing AVX2 emission math unchanged. Sibling rows are skipped entirely.
  - For observed 0/k or k/k: anchor ALT = k. For observed 0/0: select anchor via PBWT consistency (ties → lowest ALT). For missing: same PBWT rule (future: weight with PL/GL).
- Optional (future): support an “hmm” anchor scorer that prefers the ALT with highest HMM‑consistent probability instead of PBWT counts.

4) PBWT update guard
- Do not modify PBWT internals. Skip PBWT selection/updates at all supersite rows (anchors and siblings); use only biallelic sites for PBWT.

5) Output handling
- Generate phased GT fields from final hap codes at supersites.
- Re‑split to biallelic records if required by downstream compatibility, enforcing “≤1 ALT per haplotype at a position”.

6) Instrumentation and rollout
- Add counters for number of supersites, PBWT‑skipped positions, and basic timing to quantify impact.
- Gate the feature behind a build/runtime flag for incremental rollout; default to off until validated.

## Testing & Validation

- Deterministic runs
  - Omit any thread flag; pass a fixed `--seed`.

- WGS PAR2 fixtures (CHM13 T2T, 1KGP)
  - Region: `chrX:153929053-154248138`
  - Map: `test/info/par2.gmap.gz`
  - Inputs: `test/wgs/target.*.1kgp_t2t.par2.bcf`
  - Scripts: `test/scripts/phase.wgs.{family,unrelated,haploid}.sh`

- Validation tools
  - `tools/check_oneallele <out.bcf>` to count per‑position violations.
  - Integration tests compare MD5 of headerless VCFs in `test/scripts/expected/`.
  
- Supersite‑specific tests
  - No‑multiallelic regression: on biallelic‑only inputs, outputs must be bit‑for‑bit identical to current builds, and performance within noise.
  - Supersite construction unit tests: given split records at the same bp, verify deterministic code assignment and that no haplotype receives >1 ALT code.
  - Emission kernel checks: scalar vs SIMD agreement on toy cases for 0/0, 0/k, k/k.
  - Integration correctness: trio/simulated truth at known multiallelic loci; verify zero one‑allele violations and assess switch error around supersites.

## Notes for Contributors

- Scope containment matters: restrict changes to preprocessing, the new accessor, the emission kernel’s supersite branch, the PBWT guard, and output writing.
- Keep biallelic fast paths untouched to preserve regression identity and performance.
- Determinism is mandatory for tests: omit thread flags; fix `--seed`.

### Anchor selection details (v1)
- Scorer: PBWT‑consistency using current K‑state donors (count donors carrying each ALT code per hap); ties broken deterministically (lowest ALT code).
- Stability: cache anchor choice per site per iteration so forward and backward passes use the same anchors.
- Ambiguity handling: for observed 0/k (unphased), masks allow REF and hap‑anchor ALT; for observed 0/0, anchor is chosen by the scorer but observed inputs remain 0/0.
- Masks: allowed (REF or hap‑anchor ALT) = 1.0; disallowed (other ALT) = epsilon (~1e‑6).
- Flags: gate feature behind `--enforce-oneallele-supersite`; anchor scorer selectable via `--supersite-anchor=pbwt|hmm` (default: pbwt).

### “Highest HMM probability” anchor (discussion)
- If/when using HMM‑probability instead of PBWT counts, we need to define:
  - the signal to score ALTs (e.g., previous‑locus Alpha/AlphaSum weights over K‑state donors),
  - when to evaluate (pre‑emission, cached per site per iteration), and
  - performance/determinism constraints. Until then, PBWT counts are the default.

## Multiallelic Imputation Plan

- v1 (anchor‑driven, binary per‑hap posteriors)
  - Goal: enforce “≤1 ALT per haplotype per position” without changing the AVX2 core or storage.
  - Anchor selection per hap (conceptually) per supersite, realized via a single anchor row and per‑hap masks. Siblings are skipped entirely.
  - At anchors, build a per‑hap supersite mask allowing REF or hap‑anchor ALT (weight=1.0) and suppress other ALTs (epsilon). Run existing AVX2 emission as‑is. IMPUTE remains binary (REF vs chosen ALT) and we accumulate ProbMissing as today.
  - Finalization behaves the same as biallelic: set each hap allele to ALT if aggregated posterior ≥ 0.5; otherwise REF. When writing, only the chosen ALT row is set; others remain 0.
  - Anchor scorer: PBWT‑consistency (ties → lowest ALT), cached per site per iteration for stability. Feature flagged via `--enforce-oneallele-supersite` (on/off) and `--supersite-anchor=pbwt|hmm` (default pbwt).

- v2 (true multiallelic per‑hap argmax over 0..M)
  - Goal: allow the HMM to impute which ALT (if any) each hap carries at a supersite by computing posteriors over codes {0..M}.
  - IMPUTE supersite branch at anchors:
    - Use AlphaMissing/AlphaSumMissing and per‑hap prob lanes to build K‑donor weights.
    - Bucket those weights into S_c for each code c ∈ {0..M} based on donor codes (accessor), per hap.
    - Normalize to P(code=c | data) per hap, store and accumulate in a per‑site structure (ProbMissingMulti) with offsets per supersite.
  - Finalization:
    - For each supersite and hap, choose argmax c ∈ {0..M} from the aggregated ProbMissingMulti and set that hap allele to c (0=REF, k>0=ALTk). Re‑split output: only the chosen ALT row gets a 1 for that hap.
  - Performance: restricted to supersite anchors; simple bucketed sums; biallelic baseline remains bit‑identical.

## TODO (Supersite Integration Follow-up)
- Build supersite-aware haplotype accessor returning packed codes for conditioning haplotypes.
- Extend AVX2 emission kernels (and scalar fallbacks) with a supersite branch driven by the accessor data.
- Skip PBWT updates at supersites when building conditioning sets.
- Reconstruct biallelic output records from supersite hap codes during VCF/BCF writing.
- Add targeted validation (unit + integration) to ensure supersite paths preserve determinism and correctness.

## Design Decision: Supersite Iteration Strategy

- Reason
  - We must integrate multiallelic supersites into the HMM while preserving the performance and stability of the AVX2-optimized forward/backward loops. The key choice is how the iteration over loci aligns to supersites without regressing hot-path behavior.

- Option A — Full supersite index refactor
  - Description: Replace per‑variant iteration with one entry per supersite across windows, PBWT selection, transitions, and storage.
  - Implementation sketch: Make `supersites[]` the canonical locus list; rebuild `window_set`, ambiguous/missing indices, transition offsets, and Alpha buffers on supersite indices; update PBWT grouping/selection and I/O mappings.
  - Pros/Cons: Clean semantics; high blast radius and regression risk in hotspots.

- Option B — Anchor‑variant within current traversal
  - Description: Designate one anchor row per supersite. In the HMM loop, only anchors execute supersite emissions; sibling rows are fast‑skipped. PBWT updates at siblings are also skipped.
  - Implementation sketch: Precompute `is_anchor` and per‑supersite sample genotype summaries (0/0, 0/k, k/k). In emission, guard on `is_anchor`; anchors load packed codes via the accessor and reuse existing AVX2 mask/blend. Add PBWT guard for supersites; re‑split in the writer from final supersite hap codes.
  - Pros/Cons: Minimal structural change; tiny predictable branch in hot path; satellites still visited but cheap.

- Option C — Logical supersite overlay
  - Description: Keep storage unchanged but iterate an anchor list instead of raw variant indices.
  - Implementation sketch: Build an anchor index array and rework windows, ambiguous/missing, and transitions to use it; maintain mappings back to physical rows for I/O.
  - Pros/Cons: Clean loops; still invasive for window/transition machinery.

- Decision for this branch
  - Implement Option B to minimize risk to AVX2 performance and limit code churn in core loops.
  - Working branch: `multiallelic-HMM-anchor-variant-implementation`.

### Practical Challenge and Solution (Binary Emissions Preserved)

- Challenge
  - Emission helpers (INIT/RUN/COLLAPSE for HOM/AMB/MIS) are inlined, AVX2‑optimized, and assume binary per‑hap inputs (`ah` is 0/1). Switching them to multi‑valued codes would require invasive SIMD rewrites across hot paths.

- Solution (keep emissions binary)
  - Per‑sample anchor selection: For each supersite at runtime, compute the sample’s supersite genotype summary (0/0, 0/k, k/k, or missing) and choose a single anchor variant row:
    - 0/k or k/k → anchor is the row for ALT k (matches “isAltK” using existing 0/1 bit).
    - 0/0 or missing → anchor defaults to the first ALT row.
  - Binary mask synthesis at anchor while keeping the same routines:
    - For 0/k and k/k: the anchor is the correct ALT row, so the existing boolean `ah` (from Hvar.get) already equals “code==k”; AVX2 emission code remains unchanged.
    - For 0/0: at the anchor, OR all sibling ALT bits across the supersite into the anchor row for K donors so any ALT donor is treated as ALT by the binary routines. This makes hom‑ref emissions penalize any donor carrying any ALT without changing SIMD math.
  - Skip siblings and PBWT: Non‑anchor rows are fast‑skipped (no DP update), and PBWT evaluation is disabled at all supersite rows to keep PBWT binary.

- Benefits
  - No changes to the core AVX2 emission math; only a small, predictable guard per site and a byte‑level OR preprocessing at anchors for hom‑ref.
  - Biallelic performance and bit‑identity preserved.

- Notes
  - The OR preprocessing only happens on supersite anchors and only for hom‑ref samples; supersites are expected to be rare, keeping overhead small.

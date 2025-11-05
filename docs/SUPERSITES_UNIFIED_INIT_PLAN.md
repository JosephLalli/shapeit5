# SHAPEIT5 phase_common — Supersites via Unified Emission API

This document describes:
- The problem we’re solving
- The idea to solve it
- Likely pitfalls
- A high‑level architecture sketch of the changes
- A detailed, step‑by‑step implementation plan

The goal is to make supersites drive the same HMM emission kernels (`INIT_HOM/INIT_AMB/INIT_MIS`) used by biallelic sites, while keeping multiclass logic local to supersites for missing data and projection back to split records.

---

## 1) Problem

- Multiallelic loci are split into multiple biallelic records in the input, which can yield biologically invalid phasing/imputation across splits (e.g., implying multiple alternate alleles on a single haplotype).
- Current supersite handling diverges from the biallelic HMM hot paths, leading to higher error rates at multiallelic loci and duplicated logic that is harder to reason about and optimize.
- The core HMM code assumes binary emissions (match vs mismatch) using packed bit matrices and AVX2‑friendly masks, while supersites need allele class identity (ALT1..ALT15) to maintain mutual exclusivity and to impute missing genotypes coherently.
- We want to:
  - Reuse the exact same HMM emission/transition kernels for all sites.
  - Treat supersites as “binary” at emission time by collapsing to a single best class per haplotype (winner‑take‑all, per supersite and hap), while preserving full multiclass distributions only where needed (missing supersite sampling and projection).

## 2) Idea (Solution Approach)

- Introduce a site‑neutral emission adapter that provides a per‑donor/per‑lane boolean “match” mask to the HMM.
  - Biallelic sites: mask comes from existing REF/ALT bits.
  - Supersite anchors: collapse each haplotype’s expected class to one class (winner); donors match if their class equals the lane’s expected class, otherwise mismatch (REF is just another non‑matching class).
- Refactor `INIT_HOM/INIT_AMB/INIT_MIS` to generic forms that consume a “match” mask and apply the same emission math as today (match → 1.0, mismatch → ed/ee), with transitions unchanged.
- Keep supersite multiclass logic local and minimal:
  - Only at missing supersite anchors, compute and store P(class | hap) into `SC` for sampling.
  - After sampling, project the chosen classes back to split records, enforcing exactly one ALT per haplotype at that locus (mutual exclusivity), and update packed bit matrices for PBWT and subsequent iterations.
- Gate DP to the supersite anchor (`global_site_id`); sibling split records skip DP but follow projection updates.

## 3) Likely Pitfalls

- Winner consistency:
  - Forward/backward must see the same collapsed expectations; compute lane expectations deterministically and avoid recomputing with divergent inputs.
  - Define a stable tie‑break rule (e.g., lowest class index) for equal scores in missing‑data sampling when argmax is used.
- Rare‑allele shortcuts:
  - Existing fast paths keyed on “any lane matched” or “rare mismatch pattern” must still work when driven by the match mask; ensure equivalent flags are computed from the mask.
- SIMD and alignment:
  - Emission adapters must not introduce unaligned loads or gathers; allocate per‑thread scratch with 32‑byte alignment and avoid heap churn in hot paths.
- Anchor gating:
  - Never run DP on sibling splits; ensure bookkeeping (counters/indices) remains consistent so sampling/projection touches all member splits correctly.
- PBWT / packed matrices consistency:
  - After sampling, update `H_opt_hap` bits to reflect the chosen class so PBWT state selection is consistent next iteration; then transpose.
- Mutual exclusivity on projection:
  - Ensure “exactly one ALT per hap” is enforced across all split records at a supersite; clear all others.
- Performance regressions:
  - Avoid per‑site allocations; cache donor classes and match masks; reuse buffers between forward/backward.
- Underflow/double precision fallback:
  - Keep existing underflow detection and fallback behavior identical.

## 4) High‑Level Architecture Sketch

- Site‑neutral emission adapter API:
  - Builds a `SiteView` describing the current locus (biallelic, supersite anchor, supersite sibling) and the lane expectations.
  - Produces a `MatchMask` (K donors × 8 lanes) with boolean equality (donor allele == lane’s expected allele/class).
- Generic INIT wrappers in HMM:
  - `INIT_HOM_generic`, `INIT_AMB_generic`, `INIT_MIS_generic` consume `MatchMask` (or bypass for MIS) and apply the same emission multiplication as current code.
  - Existing biallelic INITs become thin wrappers that build a trivial mask and delegate to the generic forms.
- Supersite specifics kept local:
  - Compute sample hap class codes (c0, c1) and lane expectations at anchor only.
  - Cache donor class codes in an aligned buffer.
  - For missing supersites, compute multiclass posteriors (`SC`) in backward and sample one class per hap; project to splits and update packed matrices.
- Windowing, transitions, Alpha/Beta storage, pruning, and PBWT flows remain unchanged.

## 5) Detailed Step‑by‑Step Implementation Plan

1. Add site‑neutral types (header shared by single/double HMM):
   - `enum SiteKind { Biallelic, SuperAnchor, SuperSibling };`
   - `struct SiteView { SiteKind kind; int locus; const SuperSite* ss; uint8_t lane_class[8]; uint8_t emit_kind; /* HOM/AMB/MIS */ };`
   - `struct MatchMask { aligned_vector32<uint8_t> by_donor_lane; /* size K*8 */ bool any_match_lane[8]; };`
   - Ensure 32‑byte alignment for `by_donor_lane`.

2. Introduce EmissionAdapter interface:
   - `class EmissionAdapter { public: void build_view(int abs_locus, SiteView& out) const; void build_match_mask(const SiteView&, MatchMask& out) const; };`
   - Implementation must not allocate in hot paths; store per‑thread scratch in `compute_job`.

3. Implement BiallelicAdapter:
   - `build_view`: classify locus (always `Biallelic`), fill `emit_kind` from sample genotype state, derive lane expectations from `amb_code` (ref/alt expectations per lane).
   - `build_match_mask`: for each donor k and lane h, set `by_donor_lane[k, h]` to 1 if donor’s REF/ALT bit equals lane’s expected allele, else 0. Populate `any_match_lane[h]`.

4. Implement SupersiteAdapter:
   - `build_view`:
     - Map `abs_locus` to supersite index; if none → `Biallelic`.
     - If sibling split (not `global_site_id`) → `SuperSibling` and return.
     - At anchor: compute sample hap class codes `c0`, `c1` (or MIS) using existing helpers and supersite metadata.
     - Classify `emit_kind`: MIS if both missing; HOM if `c0==c1`; else AMB.
     - Build `lane_class[8]`: lanes with amb bit=0 expect `c0`, amb bit=1 expect `c1` (for HOM, both equal; for MIS, value unused).
   - `build_match_mask` (only for anchor):
     - Load/calc `donor_class[k]` once from `packed_allele_codes` at the supersite’s `panel_offset` into aligned scratch.
     - Fill `by_donor_lane[k, h] = (donor_class[k] == lane_class[h])` (HOM/AMB) or leave unused (MIS).
     - Compute `any_match_lane[h]` for rare‑allele shortcuts.

5. Refactor HMM INITs to generic forms:
   - Add `INIT_HOM_generic(const MatchMask&, float mismatch_penalty)`.
   - Add `INIT_AMB_generic(const MatchMask&, float mismatch_penalty, const Ambiguous& amb)`.
   - Add `INIT_MIS_generic()` (unchanged behavior: set emission=1.0 for all donors/lanes).
   - Make existing `INIT_HOM/INIT_AMB/INIT_MIS` call the generic versions via a trivial `MatchMask` built from current bit logic (no behavior change for biallelic).

6. Integrate adapter in forward pass (`haplotype_segment_single::forward` and double):
   - For each locus in window:
     - Build `SiteView sv` via adapter.
     - If `sv.kind == SuperSibling`, skip DP for this split (maintain counters/bookkeeping if any) and continue.
     - If `emit_kind == MIS`, call `INIT_MIS_generic()`.
     - Else build `MatchMask mm` once and call `INIT_*_generic()` accordingly.
   - Keep transition math, Alpha storage, and SIMD batch layout unchanged.

7. Integrate adapter in backward pass:
   - Mirror forward: build `SiteView`; skip siblings; reuse/rebuild the same `MatchMask` for consistent emissions; call generic INITs in reverse order as today.
   - Underflow detection and precision fallback remain as is.

8. Missing supersite multivariant imputation (anchors only):
   - When `emit_kind == MIS` at a supersite anchor, compute SC per hap:
     - Accumulate by allele class: `sum[class_c][lane] += Alpha[k, lane] * Beta[k, lane] / AlphaSum[lane]` where `class_c` is donor’s class.
     - Normalize per lane so `sum_c SC[ lane, c ] = 1.0`.
   - Store into `compute_job::SC` at `ss.class_prob_offset`.

9. Sampling and projection:
   - Sampling stage reads `SC` at missing supersite anchors:
     - Sample one class per haplotype from the multiclass distribution (or argmax with stable tie‑break for deterministic mode).
     - Project to split records: for each ALT i in the supersite, set exactly one split bit per hap to 1 if chosen class equals `i+1`, else 0; if chosen class=REF (0), clear all ALT splits for that hap.
   - Update sample genotypes and `H_opt_hap` accordingly; then run `transposeHaplotypes_H2V()` for the next PBWT iteration.

10. Maintain rare‑allele optimizations:
    - Recreate existing fast‑path conditions (e.g., early returns) from `MatchMask.any_match_lane[]` so performance of biallelic code carries over to supersites.

11. Per‑thread scratch and alignment:
    - Add to `compute_job` when supersites are enabled:
      - `aligned_vector32<uint8_t> donor_class; // size K`
      - `aligned_vector32<uint8_t> match_mask; // size K*8`
      - `std::vector<float> SC;` and `std::vector<bool> anchor_has_missing` (existing supersite plumbing)
    - Allocate once; reuse for all loci in the window; no dynamic allocation in INIT/RUN/COLLAPSE paths.

12. Determinism and tie‑breaking:
    - Seed RNG for sampling from (sample_id, supersite_id, iteration) to make runs reproducible.
    - For argmax ties, pick the smallest allele class index.

13. Instrumentation and safeguards:
    - Keep `SHAPEIT5_DEBUG_UNDERFLOW` behavior intact.
    - Add optional `SHAPEIT5_TEST_TRACE=1` hooks at supersite anchors to dump lane expectations and a checksum of the match mask to aid debugging.

14. Validation checklist:
    - Unit/targeted tests:
      - Synthetic supersite with known donor classes: verify emission parity between biallelic and supersite paths when classes collapse to REF/ALT.
      - Missing supersite: check `sum_c SC[h, c] == 1.0` and projection sets exactly one ALT per hap across splits.
      - Mixed windows: ensure sibling splits skip DP and still update correctly after sampling.
    - End‑to‑end: compare main metrics and underflow logs before/after; monitor runtime to confirm no major regressions.

15. Feature flag and rollout:
    - Guard the adapter and generic INIT path behind `--enable-supersites` (or reuse the existing flag).
    - When disabled, code paths and performance remain identical to current biallelic behavior.

---

### Notes
- “Match vs mismatch” at emission time is intentionally binary, even at supersites. REF is treated as just another non‑matching class when the lane expects a specific ALT class.
- Multiclass identity is maintained only where it matters: missing supersite imputation and projection back to split records.
- All modifications preserve the Li‑Stephens transition/emission math and SIMD vectorization patterns; we are only unifying how emissions are fed in.


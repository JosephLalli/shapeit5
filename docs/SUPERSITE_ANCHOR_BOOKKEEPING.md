# Supersite Anchor Bookkeeping — Audit and Next Fixes

This note summarizes remaining issues after unifying emissions via the site‑neutral adapter, and the next steps to fix them.

## Current Status

- Emission adapters integrated in both single and double HMM paths.
- Sibling splits gated as MIS; anchors drive `SS_*` kernels.
- Supersite missing sampling projects one ALT class per haplotype back to split records.
- Build succeeds; tests run with user `LD_LIBRARY_PATH`.

## Failing Tests (high impact)

- Multivariant posterior normalization → posteriors sum to 0 at missing supersite anchors.
- Supersite parity/state fixtures → forward/backward probability and float/double state mismatches.
- No‑double‑counting segfault → likely index/skip math during supersite projection.

## Root Causes (hypotheses)

1) Missing‑index mismatch for SC normalization (most likely):
   - Backward `IMPUTE_SUPERSITE_MULTIVARIATE` uses `curr_rel_missing` to read `AlphaMissing/AlphaSumMissing`.
   - With sibling gating, the streaming missing index can diverge from the anchor index stored during forward.
   - Result: `AlphaMissing`/`AlphaSumMissing` mismatch → all‑zero accumulators and denominators.

2) Anchor parity vs split records:
   - Tests compare anchor‑only vs both‑splits emissions. With siblings gated as MIS, parity requires compensating the anchor emission or refreshing test expectations to atomic supersite semantics.

3) Float/double parity deltas:
   - Differences at collapse/transition scaling or half‑lane packing in `SS_*` vs `INIT_FROM_MASK`.

4) No‑double‑counting segfault:
   - Off‑by‑one/skip math when advancing `vabs/vrel` across supersite members during projection.

## Next Steps

- Implement explicit missing‑index mapping for supersite anchors:
  - During forward at an anchor MIS, record `missing_index_by_locus[abs_locus] = curr_rel_missing`.
  - During backward at the same anchor, pass the recorded index into `IMPUTE_SUPERSITE_MULTIVARIATE`.
  - Apply to both single and double paths.

- Re‑run tests; if SC normalization fixes, tackle parity and segfault next.

## Acceptance Criteria

- `test_forward_backward_comprehensive` Test 4 stops reporting zero posteriors; sums ≈ 1 per lane.
- No regressions in passing tests.


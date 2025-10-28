# Chat Summary: Multiallelic Supersite HMM Integration

## Context
To ensure biological validity (≤1 ALT allele per haplotype at a genomic position) without post‑processing, we integrate multiallelic support inside the HMM by collapsing split biallelic records at the same (chr, bp) into a single supersite. This preserves the repository’s binary PBWT/storage model and only adds minimal supersite‑specific data and logic.

## Design Highlights
- Supersites: Preprocess input to build one logical site per coordinate with an allele code per haplotype (0 = REF, 1..M = ALTk). Codes are stored in a compact, bit‑packed buffer separate from existing binary structures.
- Emissions (AVX2): At supersites, decode per‑hap codes to boolean predicates (e.g., code==0, code==k) and reuse the existing two‑lane mask/blend kernels to implement 0/0, 0/k, k/k behavior. Biallelic sites use the unchanged fast path. Optional per‑site GL lookup is supported via a tiny LUT on `(min(codeA,codeB), max(codeA,codeB))`.
- PBWT: Keep PBWT strictly binary; skip PBWT updates at supersites. Normal sites update PBWT as usual.
- Transitions: Unchanged. The HMM still advances one site at a time; that site may be a supersite.
- Output: Re‑split supersites to biallelic VCF lines while enforcing “≤1 ALT per haplotype per position”.

## Implementation Plan
1. Supersite detection and deterministic ALT code assignment (stable canonicalization policy).
2. Packed per‑hap code storage with per‑site `bitwidth` and `panel_offset`; `is_super_site[]` metadata and mapping from original records to supersites.
3. Per‑site accessor returning K codes for the conditioning haplotypes; biallelic path keeps current bit‑extract logic.
4. AVX2 emission branch for supersites: SIMD compares to derive predicates, reuse existing blending, apply mismatch factor for “other ALT”.
5. PBWT guard to skip updates at supersites only.
6. Output writer to reconstruct split lines from final hap codes.
7. Instrumentation and feature flag for staged rollout.

## Testing & Validation
- Determinism: omit any thread flag; pass a fixed `--seed`.
- No‑multiallelic regression: biallelic‑only input must be bit‑for‑bit identical to current builds.
- Unit tests: supersite construction (no hap gets >1 ALT), scalar vs SIMD emission equivalence on toy 0/0, 0/k, k/k.
- Integration: trio/simulated truth at known multiallelic loci; verify zero violations and assess switch error nearby.
- Performance: measure wallclock and RSS on realistic chunks; target single‑digit runtime overhead and ~≤10% memory increase.

## Current Status / Next Steps
- Align code to the plan: add preprocessing, accessor, AVX2 supersite branch, PBWT guard, and output rewrite.
- Add counters and logs for supersite frequency and PBWT‑skipped sites.
- Run validation suite; iterate on performance after correctness is confirmed.

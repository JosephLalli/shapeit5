# Supersite Anchor Bookkeeping Fixes (`fix-ss-segfault`)

## Overview

The `fix-ss-segfault` branch resolves memory corruption triggered by heterozygous supersites. The root cause was that sibling split records were still being counted as ambiguous, despite new control flow skipping them when building the `Ambiguous` buffer. The result was an under-sized buffer that later writes overran during pruning. The fixes re-centralise supersite state and ensure every subsystem treats multiallelic loci strictly via their anchor entry.

## Key Changes

- **Central supersite cache:** `genotype::build()` now scans each supersite once, caches flags (heterozygous, scaffold, fully missing), and exposes a reusable `SuperSiteContext` helper so callers can query “is this locus the anchor?” without inspecting raw split rows.
- **Anchor-only bookkeeping:** All `build()` passes (segment counting, ambiguous array population, diplotype masks) use the cached flags and skip sibling rows entirely. Ambiguous/missing counters now represent biological sites rather than split records.
- **HMM alignment:** Forward/backward models in both single- and double-precision modes convert anchor context into emissions and imputation behaviour. Supersite siblings are treated as homozygous and no longer inflate `curr_abs_ambiguous` or `curr_abs_missing`.
- **Window + prune updates:** Window segmentation, pruning heuristics, and merge remapping rely on the new context, ensuring ambiguous counts advance exactly once per supersite.

## Validation

- Rebuilt `phase_common` successfully via `make -C phase_common`.
- Pending: run the supersite-focused unit suite (`make -C tests`) and any pipeline regression tests once the current run completes.

## Next Steps

1. Monitor the ongoing test run for regressions.
2. Roll the branch into integration once supersite tests and the main phasing pipeline pass.

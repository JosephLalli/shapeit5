## Supersite WGS Transition Debug (2025-11-19)

### Changes applied
- Added transition trace hooks for backward HMM (`SHAPEIT5_TRANS_TRACE`, optional `SHAPEIT5_TRANS_TRACE_SAMPLE`) in both single and double precision paths to log per-boundary writes.
- Added strict bounds assertions on transition writes in `haplotype_segment_double::SET_FIRST_TRANS/SET_OTHER_TRANS`.
- Hardened backward sampling buffers: dynamic sizing of `currProbs` in `genotype::sampleBackward`, index bounds checks for transition reads and sampled indices, and abort-on-OOB diagnostics.

### Findings
- For sample `HG03805`, double-precision backward now writes transition slices for every segment boundary; early iteration crash is gone.
- The repeated `[TRANS_MISMATCH_ERROR]` lines are diagnostic noise: they compare against `G->n_transitions`, which is not updated after pruning/trim. The stored counts differ from the static expectation but do not indicate an actual write failure.
- Run now reaches finalization before aborting, indicating the earlier transition-buffer underruns were addressed. Final failure is a `std::vector<unsigned char>::operator[]` assert during finalization, suggesting a late-stage sampled index overrun (likely due to dipcount/offset mismatch after pruning).

### Next steps (recommended)
1. Stop or correct the noisy transition mismatch diagnostics: recompute `n_transitions` after pruning/trim, or gate the warning behind an opt-in flag.
2. Add a thin guard/log around DipSampled access in `genotype::sample*`/`make()` to capture segment, dipcounts, sampled index, and sample name when an index is out of range during finalization.
3. Consider recalculating transition buffer sizes and `n_transitions` after pruning/trim so sampling offsets stay consistent with pruned dipcounts.

### Status
- Small supersite expansion test: already passing (per prior context).
- Large WGS integration (`test/scripts/phase.chr22.wgs.sh`): still failing at finalization due to sampled index overrun; transition buffer underruns addressed. Debug instrumentation now in place to locate the remaining issue.

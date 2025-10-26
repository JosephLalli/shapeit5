# Chat Summary: Multiallelic Enforcement Implementation

## Context
We implemented multiallelic constraint enforcement for SHAPEIT5 to ensure biological validity (≤1 ALT allele per haplotype at each genomic position). Three algorithms were implemented with increasing complexity:
1. **Transition mode**: Fast, uses nearest phased heterozygous anchors
2. **Micro mode**: Enumerates valid assignments, scores with Li-Stephens model
3. **Micro-donor mode**: Uses PBWT donor haplotypes for enhanced scoring

## Initial Implementation Issue
The original design had per-sample enforcement (`enforce_sample()`) running **after** all HMM sampling in a separate loop. This caused a critical architectural problem:

**Problem**: Donor context (PBWT Kstates) was only available during the worker thread's HMM computation phase but was needed for micro-donor scoring that happened later.

## Solution: Unified Per-Sample Enforcement
We moved enforcement to run **immediately after sampling** within each worker thread's `phaseWindow()` function:

```cpp
// In phaseWindow() - after sampling step
if (oneallele_enforcer.enabled() && multiallelic_map.size() > 0) {
    oneallele_enforcer.enforce_sample(multiallelic_map, *G.vecG[id_job], V, 
                                     threadData[id_worker].Kstates,
                                     current_iteration_context, id_job);
    // Accumulate per-sample stats into epoch stats (thread-safe)
}
```

**Benefits**:
- Donor context (Kstates) available for micro-donor mode
- Violations corrected before any downstream operations
- Works for all three modes (transition/micro/micro-donor)
- Thread-safe statistics accumulation

## Statistics Architecture
- **OneAlleleStats**: Cumulative stats across all iterations
- **OneAlleleEpochStats**: Per-iteration stats (reset each iteration)
- **sample_epoch_stats_**: Per-sample stats (reset before each `enforce_sample()` call)
- Thread-safe accumulation using mutex-protected `accumulate_sample_stats()`


### Implementation Status
The code has:
- ✅ Batch `enforce_group_micro()` (donor-agnostic, used by old `enforce()`)
- ✅ Per-sample `enforce_group_micro()` (with donor context, used by `enforce_sample()`)
- Both implementations present in the .cpp file

## Current Status
- All three enforcement modes integrated - with no testing, we do not yet know if the integrated code is working as intended.
- Per-sample enforcement with donor context working - with no testing, we do not yet know if the integrated code is working as intended.
- Thread-safe statistics accumulation implemented - with no testing, we do not yet know if the integrated code is working as intended.
- Code ready for compilation and testing

## Next Steps
1. Run integration tests with all three modes
2. Validate statistics reporting
3. Performance profiling if needed
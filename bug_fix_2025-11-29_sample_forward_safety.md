# Bug Fix: `SAMPLE_DIPC_INVALID` Crash in `genotype::sampleForward`

## 1. Issue Description
The application was crashing with the error `[SAMPLE_DIPC_INVALID][forward] ...` during the phasing process. This error occurred because the `dipcode` selected during the forward sampling step was invalid (not present in the `Diplotypes` mask for the current segment).

### Root Cause
The `currProbs` vector in `genotype::sampleForward` was initialized with a fixed size of 64. The `rng.sample` function iterates through the entire vector provided to it. Due to floating-point precision issues, the accumulated probability sum (`sumProbs`) could be slightly less than the random target value. This caused `rng.sample` to iterate past the valid diplotypes (index >= `curr_dipcount`) into the zero-initialized tail of the vector, eventually returning an out-of-bounds index (e.g., 63). This index corresponded to an invalid diplotype, triggering the safety check and subsequent abort.

This bug affects both **biallelic** and **multiallelic** paths, though it was previously silent in the `main` branch (release 5.1) because the explicit `DIP_GET` safety check was absent. The silent failure likely led to incorrect phasing decisions (falling back to dipcode 0) rather than a crash.

## 2. The Fix
We modified `phase_common/src/objects/genotype/genotype_sweep.cpp` to explicitly resize the `currProbs` vector to `curr_dipcount` before passing it to `rng.sample`.

```cpp
// In genotype::sampleForward
curr_dipcount = countDiplotypes(Diplotypes[s]);
currProbs.resize(curr_dipcount); // <--- The Fix
```

This ensures that `rng.sample` only iterates over valid diplotypes. If the random target exceeds the sum of probabilities (due to precision), it will naturally clamp to the last *valid* element instead of running off into invalid memory.

## 3. Verification
*   **Reproduction:** A new unit test `tests/src/test_genotype_sampling_safety.cpp` was created to reproduce the vector sizing issue in isolation.
*   **Unit Test:** The new test passes with the fix applied.
*   **Integration Test:** The script `test/scripts/phase.chr22.wgs.sh` now completes the biallelic phasing run successfully.

## 4. Impact
This is a critical stability and correctness fix. It prevents crashes in debug builds and prevents silent logic errors in release builds where invalid diplotypes were previously being selected.

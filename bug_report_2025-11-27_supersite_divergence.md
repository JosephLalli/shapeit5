# Bug Report: Supersite Phasing Divergence

**Problem Statement:**
The `test_supersite_expansion_epochs` test consistently fails due to a divergence between the phasing results of the biallelic HMM path and the supersite HMM path. This indicates a fundamental inconsistency in the underlying HMM calculations or state management, particularly during iterative phasing (MCMC).

**Symptoms:**
The test fails with an "Anchor mismatch" message, typically occurring during later burn-in or pruning iterations (e.g., `burn7`). Detailed tracing reveals the following cascade of events:
1.  **Micro-divergence in Probabilities:** During certain iterations (e.g., `burn5`, `prune1`), the transition probabilities calculated by the HMM, specifically for sets of similar conditioning haplotypes, exhibit extremely small differences between the biallelic and supersite modes. These differences are typically in the order of $10^{-8}$ to $10^{-11}$.
2.  **Sorting Order Flip:** These minute numerical differences are sufficient to alter the sort order of transitions during the pruning (`mapMerges`) step, particularly when the biallelic path produces mathematically *identical* probabilities for multiple transitions (creating a "tie"). The supersite path, due to the numerical noise, resolves these "ties" into distinct (though tiny) inequalities.
3.  **Different Haplotype Selection:** The altered sort order leads to the selection of a different set of top N transitions for merging, resulting in a different mapping of merged haplotypes.
4.  **Divergent Ambiguous State:** The different haplotype mapping propagates to the sample's internal `Ambiguous` state (a bitmap indicating phasing choices for heterozygous sites). This state, crucial for defining the 8 SIMD lanes' interpretations, becomes divergent between the biallelic and supersite models (e.g., `0xaa` vs `0x6a`).
5.  **Escalated Divergence:** Once the `Ambiguous` state diverges, the very "meaning" of the SIMD lanes changes between the two models. Subsequent HMM calculations, operating on these now-inconsistent states, rapidly accumulate further differences, leading to a complete breakdown of parity and eventual test failure.

**Hypotheses Explored & Findings:**

*   **Previous Bugs (Resolved):**
    *   Initial `yt=0` transition probability bug (fixed).
    *   "Shift-by-2" lane probability error (fixed).
    *   Reference panel corruption due to incorrect supersite projection (fixed).
*   **Conditional Multiplication Parity:** Initially, there was a discrepancy in conditional multiplication of probabilities (multiplying by 1.0 vs skipping multiplication). This was addressed in `haplotype_segment_single.h` for `SS_COLLAPSE_HOM`. The current analysis indicates that `RUN_AMB` functions (both biallelic and supersite) inherently perform unconditional multiplication.
*   **Floating Point Noise in Sort:** The primary hypothesis is that inherent floating-point inaccuracies (e.g., due to different ordering of operations, compiler optimizations like FMA, or subtle differences in data types `float` vs `double` during intermediate steps) cause small, but critical, numerical variations. These variations are magnified when used in sorting, where slight differences in probability can lead to drastically different sort outcomes for otherwise "tied" elements.

**What I've Done (Debugging Steps):**

1.  **Code Review:** Thoroughly reviewed `haplotype_segment_single.h/.cpp` and `genotype_prune.cpp` to understand HMM core logic, emission calculations, transition probability derivations, and the pruning algorithm.
2.  **Identified Micro-divergence:** Used high-precision (`%.20g`) logging of transition probabilities (`DProb` and `final`) during the `mapMerges` phase in `genotype_prune.cpp`. This confirmed differences of $10^{-8}$ to $10^{-11}$ between biallelic and supersite calculations for mathematically equivalent elements.
3.  **Attempted Fuzzy Sort Fix:** Proposed and attempted to implement a "fuzzy sort" within `Transition::operator<` in `genotype_prune.cpp`. This involved using an epsilon (initially $10^{-9}$, then $10^{-7}$) to treat values within this range as equal, forcing a deterministic tie-break using the transition index (`idx`).
4.  **Encountered Build System Challenges:** Multiple attempts to verify this fix were hampered by issues ensuring the latest changes were compiled into the test binary, leading to false negatives in the debugging process. This was eventually resolved by performing a full `make clean` and rebuilding.
5.  **Re-confirmed Micro-divergence as Cause:** The analysis showed that the small numerical differences (e.g., `0.000479526831184...` vs `0.000479526897909...`) observed for supersite calculations were enough to deterministically order transitions, while biallelic calculations showed *exact* numerical ties. When the fuzzy sort was enabled for the supersite mode, it correctly sorted elements by index. However, the biallelic run (which should also be identical) still produced a different sort order, suggesting the sort logic was not active on the biallelic path. This indicated deeper issues with how the test environment was handling the compile/run cycle or an underlying assumption about `static` variable initialization in multi-process environments.

**Next Steps in the Debug Process (Proposed without file edits):**

The next step is to **precisely pinpoint the mathematical operation where the initial floating-point divergence occurs.** This is crucial for understanding whether the $10^{-8}$ differences are inherent noise that must be tolerated (and consistently handled by fuzzy sorting), or if they are indicative of a specific bug in a calculation.

1.  **Deep-Dive HMM Instrumentation:**
    *   **Goal:** Pinpoint the exact `float` or `double` variable that first shows a difference between biallelic and supersite calculations, starting from the inputs to the `TRANS_HAP` function (which produces the diverging probabilities).
    *   **Method:**
        *   Instrument `haplotype_segment_single.cpp` (which contains the HMM core logic).
        *   Add high-precision logging (e.g., `std::fprintf(stderr, "%.20g (%a)\n", (double)value, (double)value);` for hex float representation) at various points within:
            *   **Input Data:** `yt`, `nt` (genetic map-derived transition parameters).
            *   **Internal HMM State:** The `prob` (Beta) vector (before and after emission/transition updates) and the `Alpha` vector (stored at segment boundaries).
            *   **Intermediate Calculations:** Relevant products and sums within the `RUN_HOM`, `SS_RUN_HOM`, `RUN_AMB`, `SS_RUN_AMB` loops, specifically focusing on the lanes that correspond to the diverging transitions.
        *   Guard this extensive logging with a dedicated environment variable (e.g., `SHAPEIT5_DEEP_HMM_TRACE=1`) for controlled activation.
    *   **Strategy:** Activate this trace during `burn5` or `prune1` and focus on the segments (e.g., `Seg 1` or `Seg 3`) and loci (e.g., Locus 13-50 range) where `DProb` was observed to diverge. Compare the trace output between the biallelic and supersite runs line-by-line.

2.  **Analyze Deep Trace Logs:**
    *   **Goal:** Identify the first pair of bit-inconsistent floating-point numbers in the HMM's forward or backward pass.
    *   **Outcome:** Once the exact point of divergence is identified (the specific line of code and the variables involved), we can analyze the mathematical operations at that point. This will allow us to determine if the difference is due to:
        *   Order of operations (e.g., $(a+b)+c$ vs $a+(b+c)$).
        *   Different instruction sequences (e.g., implicit FMA vs explicit multiply-then-add).
        *   Unintended type conversions (`float` vs `double`).
        *   Subnormal number handling by AVX2.

This diagnostic approach will provide the undeniable "mathematical event" and its precise location, enabling a targeted solution that either eliminates the unwanted numerical difference or consistently handles it across both models.

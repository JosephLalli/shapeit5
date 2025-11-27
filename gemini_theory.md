Current Diagnosis: Probability Distribution Divergence

  The debugging session has revealed a critical divergence in the distribution of probabilities (prob array) between the
  biallelic and supersite HMM pathways, even when the overall sum of probabilities (probSumT) remains momentarily
  consistent. This divergence is evident from the first few iterations and ultimately leads to discrepancies in
  ambiguous masks and downstream haplotype sampling.

  Key Findings and Observations:

   1. Initial Probability Distribution Mismatch (Locus 16/32, Iteration 5):
       - At Locus 16 (Biallelic) and its corresponding Locus 32 (Supersite Anchor) in Iteration 5 (Burn5), a significant
         divergence in individual probability states (prob[0]) was observed:
           - Biallelic prob[0] (for conditioning haplotype 0, lane 0): 0.0039062500 (which is 1/256).
           - Supersite prob[0] (for conditioning haplotype 0, lane 0): 0.0000390625 (which is 1/25600, or 100 times
             smaller).
       - Crucially, at this same point, the probSumT (total probability mass across all states) for both paths was
         identical: 0.5050001144.
       - This indicates that while the total probability mass is conserved, its distribution across the 256 HMM states
         (32 conditioning haplotypes × 8 SIMD lanes) is drastically different. To compensate for prob[0] being 100x
         smaller in the supersite path, other prob[x] values must be significantly larger to maintain the same probSumT.

   2. Ambiguous Mask Divergence (Locus 11/22, Iteration 9):
       - Further downstream, at Locus 11 (Biallelic) and Locus 22 (Supersite Anchor) in Iteration 9 (Burn7), the
         ambiguous masks (amb_code/amb_mask) used for emission calculations were found to be different:
           - Biallelic amb_code: 0x65 (binary: 0110 0101)
           - Supersite amb_mask: 0x6a (binary: 0110 1010)
       - The difference is particularly striking in the lower 4 bits (lanes 0-3), where the bits are effectively
         inverted (0101 vs 1010). The upper bits (lanes 4-7) match.
       - This divergence in ambiguous masks directly leads to different emission probabilities for individual SIMD lanes
         in the RUN_AMB function, which in turn causes further divergence in the prob array distribution and eventually
         probSumT.

   3. Causality Chain Established:
       - The initial difference in prob distribution (observed as early as Locus 16/32 in Iteration 5) leads to subtle
         (or not-so-subtle) differences in the calculated posterior probabilities.
       - These differing probabilities influence the segment merging decisions made during pruning stages (e.g., in
         Prune1/Iteration 6).
       - Different merging decisions result in different combinations of parent haplotypes for merged segments, which
         then manifest as divergent Ambiguous masks stored in the genotype (G->Ambiguous).
       - Subsequently, during later iterations (e.g., Burn7/Iteration 9), these divergent Ambiguous masks are retrieved
         by the HMM, leading to significantly different emission probability calculations (0x65 vs 0x6a), which then
         cause the overall HMM likelihoods to diverge, ultimately resulting in the test failure.

  Root Cause Hypothesis:

  The primary root cause of the entire divergence chain is a fundamental difference in how the probability distribution
  across states is maintained or updated between the biallelic and supersite paths, even when the aggregate sum
  (probSumT) is momentarily consistent. This suggests one or both of the following:

   1. Panel Data Mismatch at Specific Loci: The underlying data (alleles) provided by the Hvar.get() function for
      biallelic paths and the ss_cond_codes (derived from panel_codes) for supersite paths may not be perfectly
      equivalent for the respective loci, leading to different match/mismatch patterns during emission calculation for
      individual states. The prob[0] difference of 100x while probSumT matches is a strong indicator that specific
      conditioning haplotypes are assigned drastically different probabilities.

   2. Subtle Differences in State Management/Indexing: Despite efforts to ensure parity, there might still be subtle
      differences in how states are indexed, aligned, or how probability mass is transferred or accumulated when
      transitioning between different variant types (biallelic vs. supersite anchor vs. supersite sibling), especially
      given the repeat_factor logic used in the test. The "Factor of 100" in prob[0] could indicate a scaling issue or
      an incorrect indexing that causes a specific haplotype's probability to be suppressed.

  Next Steps for Further Diagnosis (without altering code):

  To pinpoint the exact origin of the 100x prob[0] divergence at Locus 16/32, the next step would be to:

   1. Detailed Panel Code Comparison: Instrument the RUN_HOM/SS_RUN_HOM (or RUN_AMB/SS_RUN_AMB) functions at Locus 16
      (Biallelic) and Locus 32 (Supersite) in Iteration 5 to print:
       * The ag (sample allele) value.
       * The ah (conditioning haplotype allele) value for a few key k (conditioning haplotype index).
       * The ss_cond_codes[k] value for supersites.
       * The expected emission factor for each path.
      This would allow for a direct comparison of the input data used for emission calculation.

  The current diagnosis indicates that the problem is deeply rooted in the handling and representation of genotype data
  and conditioning haplotypes between the two HMM variants, which then cascades to affect ambiguous masks and ultimately

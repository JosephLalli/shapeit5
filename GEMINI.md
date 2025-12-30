# Gemini's Understanding of the SHAPEIT5 Repository

This document summarizes the key aspects of the SHAPEIT5 repository as understood by the Gemini agent.

## 1. Project Overview
**SHAPEIT5** is a suite of tools primarily focused on genetic phasing and imputation. Its core functionality revolves around implementing variations of the MCMC Li-Stephens model to determine the allele of origin for heterozygous variants and impute missing genotypes. It specifically supports biallelic sites, internally representing variants as 0 (reference) or 1 (not reference).

The repository is structured into several sub-projects, each addressing a specific aspect of the phasing workflow or related utilities:
*   **`phase_common`**: The core algorithm for phasing common variants using a full HMM. This was the primary focus of the investigation.
*   **`phase_rare`**: A mode for phasing rare variants onto existing common variant scaffolds using a simplified HMM.
*   **`ligate`**: A tool likely used for combining phased segments or regions.
*   **`switch`**: A tool probably used for calculating switch error rates, a metric for phasing accuracy.
*   **`simulate`**: A tool for generating simulated genetic data, including family structures (trio/duo phasing) and haploid samples, for testing and development.
*   **`xcftools`**: A collection of utilities for manipulating VCF/BCF files (e.g., concatenation, tag filling, viewing).

## 2. Core Algorithm: `phase_common`
The `phase_common` algorithm is an iterative MCMC approach with three stages: Burn-in, Pruning, and Main. Each iteration involves:
1.  **PBWT State Selection**: Efficiently selects `K` (approx. 200) "best" donor haplotypes per sample using the Positional Burrows-Wheeler Transform (PBWT) to reduce computational complexity.
2.  **Window Segmentation**: Divides the genome into overlapping windows for parallel processing.
3.  **Multi-threaded HMM Computation**: The core of the algorithm, performed independently for each sample in parallel.
    *   **Forward Pass**: Computes alpha values (probability of observations up to a variant given a state).
    *   **Backward Pass**: Computes beta values (probability of observations after a variant given a state), posterior probabilities, and imputes missing genotypes.
4.  **Sampling**: Selects the most likely phase configuration based on HMM posteriors.
5.  **Haplotype Update**: Updates the reference panel with newly sampled phase configurations.

### Key Optimizations & Features:
*   **AVX2 SIMD**: Heavy use of Intel AVX2 intrinsics (`_mm256_*`) to process 8 lanes (floats) or 4 lanes (doubles) simultaneously, significantly accelerating HMM computations.
*   **Precision Fallback**: Starts with single-precision (faster) and automatically retries with double-precision if numerical underflow occurs, ensuring stability.
*   **Supersites**: An extension to handle multiallelic variants atomically. It groups split multiallelic records at the same genomic position, ensuring biological validity during phasing and imputation. This involves specialized data structures, emission calculations, and post-HMM projection.
*   **Outer Product Seeding**: An enhancement to carry both row and column marginals at segment boundaries, improving accuracy, especially for supersites.

## 3. Data Structures & Representation
*   **Bit-packed Genotypes**: Genotypes are stored in a compact, bit-packed format (`genotype::Variants`) using macros like `VAR_GET_HET`, `VAR_SET_SCA`, etc., for efficient storage and manipulation.
*   **Graph-based Genotypes**: Each sample's genotype is represented as a directed acyclic graph (DAG) of segments, with diplotypes encoded as 64-bit masks.
*   **`Ambiguous` Array**: Stores bitmasks for heterozygous sites, defining the 8 different phase configuration hypotheses explored by the SIMD lanes.
*   **`conditioning_set`**: Manages the reference panel haplotypes (`H_opt_hap`, `H_opt_var`) and the selected PBWT neighbors (`indexes_pbwt_neighbour`).

## 4. Build System
*   **Makefiles**: The project uses a hierarchical Makefile system. The root `makefile` dispatches to sub-project Makefiles, which in turn include `common/makefile_common.mk`.
*   **Compiler**: `g++ -std=c++17`.
*   **Flags**: `-mavx2 -mfma` are standard, enabling AVX2 intrinsics. Optimization levels vary (`-O0` for default, `-O2` for static builds, `-g` for debug).
*   **Dependencies**: Relies on `Boost` libraries (iostreams, program_options, serialization) and `HTSlib` for VCF/BCF I/O.
*   **`aligned_vector32`**: Custom `std::vector` alias using `boost::alignment::aligned_allocator` to ensure 32-byte memory alignment for AVX2 operations.

## 5. Testing Strategy
The project employs a robust testing framework:
*   **Unit Tests**: Located in `tests/src/`, these are C++ executables designed to validate specific components or algorithms (e.g., `test_supersite_representation_parity`, `test_supersite_outer_product`).
*   **File-based Tests**: The `test/` directory contains numerous BCF files (`.bcf.csi`) and log files, serving as test cases and expected outputs for various scenarios.
*   **Tracing**: Debug logging mechanisms (e.g., `SHAPEIT5_DEBUG_UNDERFLOW`) are available to aid in debugging and verifying numerical stability and algorithmic parity.
*   **Parity Checks**: A significant focus is on ensuring that the supersite implementation produces results identical to or consistent with the biallelic implementation, especially at anchor loci.

## 6. Documentation
*   **`AGENTS.md`**: The primary source of truth for the project's architecture, algorithms, and recent developments, mirroring `.github/copilot-instructions.md`.
*   **Jekyll Site**: The `docs/` directory contains a Jekyll-based documentation website with detailed explanations, a glossary, and specific feature outlines (e.g., `.AGENT_markdowns/SUPERSITE_ANCHOR_BOOKKEEPING.md`).

## 7. Deployment & Usage
*   **Docker**: The `docker/` directory provides a `Dockerfile` and scripts (`build.sh`, `upload.sh`) for creating portable Docker images, facilitating deployment.
*   **Scripts & Tasks**: The `scripts/` and `tasks/` directories contain specialized scripts for running the phasing pipeline on large-scale datasets like the UK Biobank, indicating its use in significant academic research.
*   **External Data**: Relies on external genetic map files (e.g., `resources/maps/b38/*.gmap.gz`) for accurate recombination modeling.

## 8. Key Learnings & Insights
*   The project is a highly optimized, high-performance bioinformatics tool.
*   It leverages modern C++ features, SIMD instructions (AVX2), and careful memory management for efficiency.
*   The "supersite" feature is a complex but critical extension to handle multiallelic variants correctly, requiring significant architectural considerations and rigorous testing.
*   The development process appears to be mature, with a strong emphasis on testing, debugging, and detailed internal documentation.

This comprehensive understanding will enable effective assistance with future tasks related to the SHAPEIT5 project.

---

## 9. Detailed Breakdown of the `phase_common` Algorithm

This section provides a deeper dive into the structure and implementation of the `phase_common` algorithm.

### 9.1. High-Level Algorithmic Flow

The algorithm is driven by the `phaser::phase()` method in `phaser_algorithm.cpp`. It follows a well-defined MCMC scheme:

1.  **MCMC Iteration Loop**: The main loop iterates through a pre-defined number of "burn-in", "pruning", and "main" iterations.
    ```cpp
    // In phaser::phase()
    for (iteration_stage = 0 ; ... ) {
        for (int iter = 0 ; iter < iteration_counts[iteration_stage] ; iter ++) {
            // ...
        }
    }
    ```

2.  **PBWT State Selection**: At the beginning of each iteration, it calls `H.select()` to run the PBWT algorithm. This selects the `K` most relevant conditioning haplotypes (states) for each sample's haplotypes, which is a major performance optimization.

3.  **HMM Window Processing**: It then calls `phaseWindow()`, which orchestrates the HMM computation across all samples, typically in parallel.
    *   For each sample, the genome is broken into windows.
    *   For each window, the HMM forward-backward algorithm is run. The implementation uses a precision-fallback mechanism:
        ```cpp
        // In phaser::phaseWindow(int id_worker, int id_job)
        // Try single precision first
        haplotype_segment_single HS(...);
        HS.forward();
        outcome = HS.backward(...);

        // If underflow, retry with double precision
        if (outcome != 0) {
            haplotype_segment_double HS_double(...);
            // ...
            G.vecG[id_job]->double_precision = true;
        }
        ```
    *   After the HMM, depending on the MCMC stage, it either samples a new phase (`sample()`), merges high-confidence segments (`mapMerges`/`performMerges`), or stores the posterior probabilities (`store()`).

4.  **Haplotype Update**: After all samples have been processed for the iteration, the main `conditioning_set` (`H`) is updated with the newly sampled haplotypes from the `genotype_set` (`G`) via `H.updateHaplotypes(G)`. The `H` matrix is then transposed to prepare for the next PBWT selection.

### 9.2. File-by-File Implementation Role

| File(s)                                                     | Role in the Algorithm                                                                                                                                                           |
| ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `phaser/phaser_algorithm.cpp`                               | Orchestrates the main MCMC loop, dispatches HMM jobs to worker threads, and calls the major sub-system functions (PBWT selection, Haplotype update).                               |
| `containers/conditioning_set/conditioning_set_selection.cpp` | Implements the PBWT algorithm (`H.select()`) to efficiently select the `K` conditioning states for the HMM.                                                                       |
| `models/haplotype_segment_{single,double}.h`                | The computational core. Implements the forward-backward algorithm using AVX2 intrinsics. Contains the logic for `INIT`, `RUN`, and `COLLAPSE` operations for all variant types. |
| `objects/genotype/genotype_build.cpp`                       | Pre-processes the raw variant data for a sample into the graph representation (segments and diplotypes) that the HMM operates on.                                                 |
| `objects/genotype/genotype_managment.cpp`                   | Implements the `sample()`, `store()`, and `solve()` methods, which use the HMM posterior probabilities to determine the final phased haplotypes.                                    |
| `objects/super_site_builder.cpp`                            | Implements the logic to detect and group split multiallelic variants into `SuperSite` objects before the HMM iterations begin.                                                    |
| `phaser/phaser_finalise.cpp`                                | Contains the final steps of the pipeline, including calling `G.solve()` and `projectSupersites()`, and writing the final phased haplotypes to an output file.                     |
| `io/genotype_reader/genotype_reader_reading.cpp`            | Handles the initial reading of VCF/BCF files, parsing genotypes, and encoding them into the initial bit-packed format.                                                            |

### 9.3. Key Data Structures and Their Purpose

| Data Structure        | File(s)                               | Purpose                                                                                                                                                                                                                         |
| --------------------- | ------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `genotype`            | `objects/genotype/genotype_header.h`  | Represents a single individual's genotype data as a graph of segments. It holds the bit-packed `Variants` array, the `Ambiguous` masks, and the `Diplotypes` masks.                                                              |
| `conditioning_set`    | `containers/conditioning_set/*.h`     | Manages the full reference panel (`H_opt_hap`, `H_opt_var`) and the results of the PBWT selection (`indexes_pbwt_neighbour`), providing the `K` states for the HMM.                                                              |
| `variant_map`         | `containers/variant_map.h`            | Stores metadata for all variants, such as chromosome, position, and genetic map distance (`cm`), which is crucial for calculating transition probabilities.                                                                    |
| `hmm_parameters`      | `objects/hmm_parameters.h`            | A container for the core mathematical parameters of the HMM: transition probabilities (`t`, `nt`), emission probabilities (`ed`, `ee`), and effective population size (`Neff`).                                                   |
| `SuperSite`           | `objects/super_site_accessor.h`       | A struct that defines a multiallelic site, grouping multiple variant indices (`var_start`, `var_count`) that belong to the same biological site and defining the anchor variant (`global_site_id`).                               |
| `aligned_vector32<T>` | `utils/otools.h`                      | A `std::vector` using a Boost aligned allocator. This is used for all arrays involved in SIMD computations (e.g., `prob`, `probSumH`) to guarantee 32-byte memory alignment required by AVX2.                                     |

### 9.4. Critical Lines of Code and Concepts

*   **SIMD Lane Hypothesis Generation** (`objects/genotype/genotype_build.cpp`): This loop creates the `Ambiguous` bitmask that gives each of the 8 SIMD lanes a different phasing hypothesis to test. `n_unf` tracks the number of HETs already seen in the segment.
    ```cpp
    if (f_het) {
        for (unsigned int h = 0 ; h < HAP_NUMBER ; h ++) {
            bool allele = ((h>>n_unf)%2);
            if (allele) HAP_SET(Ambiguous[a1], h);
        }
        n_unf++;
    }
    ```

*   **Core HMM Transition Logic** (`models/haplotype_segment_single.h`): This shows the Li & Stephens model implemented with AVX2 intrinsics. The probability of the next state is a combination of staying with the same donor (`_prob * _nt`) or switching to any other donor (`_tFreq`).
    ```cpp
    // In RUN_HOM
    __m256 _prob = _mm256_load_ps(&prob[i]);
    _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq); // Fused Multiply-Add: prob = (prob * nt) + tFreq
    if (ag!=ah) _prob = _mm256_mul_ps(_prob, _mismatch);
    ```

*   **Supersite Anchor Gating** (`models/haplotype_segment_single.h`): This logic ensures that the main HMM dynamic programming is only run once for an entire multiallelic site—at its designated anchor variant. Sibling variants are skipped.
    ```cpp
    // In INIT_HOM, RUN_HOM, etc.
    int ss_idx = (*locus_to_super_idx)[curr_abs_locus];
    if (ss_idx >= 0) {
        const SuperSite& ss = (*super_sites)[ss_idx];
        if (curr_abs_locus != (int)ss.global_site_id) {
            // Sibling: true no-op or neutral initialization
            return;
        }
        // ... run full supersite HMM logic at anchor ...
    }
    ```

*   **Outer Product Seeding** (`models/haplotype_segment_double.h`): This function prepares the mixing probabilities for carrying over both row and column marginals from the previous segment, improving transition accuracy.
    ```cpp
    bool haplotype_segment_double::prepare_outer_product_mix(...) {
        // ...
        const __m256d stay = _mm256_set1_pd(nt / prev_total);
        const __m256d switch_vec = _mm256_set1_pd(yt / static_cast<double>(HAP_NUMBER));
        col_mix_lo = _mm256_fmadd_pd(prev_lo, stay, switch_vec);
        col_mix_hi = _mm256_fmadd_pd(prev_hi, stay, switch_vec);

        row_stay = nt / prev_total;
        row_switch = yt / static_cast<double>(n_cond_haps);
        return true;
    }
    ```
## 10. Supersite Phasing Accuracy Bug [RESOLVED]

This section documents a previously critical bug that affected the accuracy of the supersite phasing algorithm, which has since been resolved.

### 10.1. The Original Problem: Reference Panel Corruption

The core of the issue was a disconnect between the multiallelic nature of the supersite HMM and the strictly biallelic representation required by the PBWT algorithm and the main reference panel (`H_opt_hap`).

After the HMM made a phasing decision for a multiallelic supersite, the program failed to correctly "project" this decision back onto all the individual biallelic variants that constitute the supersite *during* the MCMC iterations. The `genotype::Variants` array was left in an inconsistent state, which then corrupted the main reference panel (`H_opt_hap`) when it was updated in the next iteration. This led to a progressive degradation of phasing accuracy.

### 10.2. The Resolution

A code review confirmed that this bug has been fixed. The codebase now contains the correct logic to prevent reference panel corruption. The fix consisted of two main parts:

1.  **Correct `projectSupersites` Implementation:** The function `genotype::projectSupersites` in `genotype_managment.cpp` was updated to correctly handle an arbitrary number of alternate alleles. It now properly translates the unified, multi-allelic phasing decision from a supersite's anchor variant back to the correct bit-packed representation for all underlying biallelic variants.

2.  **Correct Call Timing:** The `projectSupersites` function is now called at the end of both `genotype::sampleForward` and `genotype::sampleBackward` in `genotype_sweep.cpp`. This ensures that the projection occurs immediately after a new haplotype configuration is sampled and, crucially, *before* the main reference panel is updated.

These changes ensure that the reference panel remains consistent and biologically correct throughout the MCMC iterations, resolving the root cause of the accuracy degradation. A final call to `projectSupersites` in `phaser_finalise.cpp` also ensures the final output file is consistent.

## 11. Debugging the `test_supersite_expansion_epochs` Divergence

This section tracks the progress of debugging a persistent accuracy divergence between the biallelic and supersite HMM pathways, observed in the `test_supersite_expansion_epochs` test.

### 11.1. Initial Bug: `yt=0` Transition Probability [RESOLVED]

*   **Symptom:** The initial investigation revealed that for supersite anchors, the backward pass was using a transition probability `yt` of `0.0`. This effectively turned off the model's ability to account for recombination between haplotypes for supersites, causing a major divergence in probability calculations.
*   **Root Cause:** The `DEFERRED_INIT` block within the `backward()` method of `haplotype_segment_single.cpp` and `haplotype_segment_double.cpp` contained a logical flaw. When the backward pass started on a sibling variant, this block would incorrectly decrement `prev_abs_locus`. As a result, when the HMM reached the anchor, `prev_abs_locus` pointed to the adjacent sibling, which has a genetic distance of zero, leading to `yt=0`.
*   **Resolution:** The erroneous `prev_abs_locus--` line was removed from the `DEFERRED_INIT` block in both files. This corrected the sibling-skipping logic and resolved the `yt=0` issue.

### 11.2. Current Bug: "Shift-by-2" Lane Probability Error

*   **Symptom:** After fixing the `yt` bug, tests still fail. The `test_supersite_expansion_epochs` test now provides a clear diagnostic sign: all supersite sumH values are the same as biallelic values, but shifted by 2 positions. This indicates that while the total probability is conserved, it is being incorrectly distributed across the 8 SIMD lanes. `probSumH` (the array of per-lane probabilities) in the supersite path is a shifted version of the correct biallelic result.
*   **Current Hypothesis:** This is a logical indexing error, not a floating-point issue. The error is almost certainly located in the handling of ambiguous heterozygous sites (`AMB` path), specifically in how the 8-bit `amb_mask` is used to assign phasing hypotheses to the 8 SIMD lanes in `SS_RUN_AMB`. A simple off-by-two error or an incorrect bit-shift operation when processing this mask could explain the observed perfect 2-position shift.
*   **Next Steps:** The immediate goal is to pinpoint the source of this shift. This requires adding detailed, per-lane instrumentation to both `RUN_AMB` (biallelic) and `SS_RUN_AMB` (supersite) to trace how the `amb_mask` is interpreted and how the final emission vectors are constructed for each lane. Comparing the trace output from the two paths will reveal the exact line of code where the logic diverges.

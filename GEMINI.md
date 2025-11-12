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
*   **Tracing**: Extensive debug logging and tracing mechanisms (e.g., `SHAPEIT5_TEST_TRACE`, `SHAPEIT5_DEBUG_UNDERFLOW`) are built into the HMM to aid in debugging and verifying numerical stability and algorithmic parity.
*   **Parity Checks**: A significant focus is on ensuring that the supersite implementation produces results identical to or consistent with the biallelic implementation, especially at anchor loci.

## 6. Documentation
*   **`AGENTS.md`**: The primary source of truth for the project's architecture, algorithms, and recent developments, mirroring `.github/copilot-instructions.md`.
*   **Jekyll Site**: The `docs/` directory contains a Jekyll-based documentation website with detailed explanations, a glossary, and specific feature outlines (e.g., `SUPERSITE_ANCHOR_BOOKKEEPING.md`).

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
## 10. Critical Bug Report: Supersite Phasing Accuracy

This section details a critical bug affecting the accuracy of the supersite phasing algorithm, leading to a significant drop in performance compared to the standard biallelic algorithm.

### 10.1. The Problem: Reference Panel Corruption

The core of the issue lies in the disconnect between the multiallelic nature of the supersite HMM and the strictly biallelic nature of the PBWT algorithm and its main reference panel (`H_opt_hap`).

1.  **Biallelic PBWT:** The PBWT algorithm, used for selecting conditioning haplotypes (`K`), is designed to work on biallelic data (0 for REF, 1 for ALT). It is not aware of multiallelic supersites.
2.  **Multiallelic HMM:** The HMM correctly treats a supersite as a single entity, making a phasing decision for the entire site at its designated "anchor" variant.
3.  **Inconsistent State:** After the HMM sampling step (`genotype::sample`), the phasing decision is applied, but only the anchor variant's entry in the `genotype::Variants` array is guaranteed to be correct. The other constituent biallelic variants that form the supersite are not updated. This leaves the `Variants` array in an inconsistent state, where the biallelic representation does not match the true, phased multiallelic haplotype.
4.  **Reference Panel Corruption:** During each MCMC iteration, the `haplotype_set::updateHaplotypes` function is called. This function blindly copies the data from the (now inconsistent) `Variants` array into the main `H_opt_hap` reference panel.
5.  **K-inflation and Accuracy Loss:** This process pollutes the reference panel with nonsensical haplotypes. When the PBWT algorithm runs in the next iteration, it struggles to find good matches in the noisy panel, causing it to increase the number of selected states (`K`). This leads to a downward spiral of increasing `K` and decreasing phasing accuracy.

### 10.2. Design Intent of `projectSupersites`

The `genotype::projectSupersites` function was designed to be the crucial link between the multiallelic HMM and the biallelic PBWT.

*   **Purpose:** Its sole purpose is to "project" the unified phasing decision made at a supersite's anchor variant back onto all the individual biallelic variants that constitute the supersite.
*   **Mechanism:** For example, if the HMM determines a supersite's phasing is `allele 3 | allele 2`, `projectSupersites` is responsible for translating this into the correct bit-packed representation for the underlying biallelic variants (e.g., setting them to `0|0`, `0|1`, `1|0` in the `Variants` array).
*   **Goal:** By ensuring the `Variants` array is fully consistent *before* `updateHaplotypes` is called, it guarantees that the main reference panel remains clean and biologically correct, allowing the biallelic PBWT to function as intended.

### 10.3. Implementation Status and Flaws

A review of the codebase revealed two critical flaws that cause the bug:

1.  **Incorrect Call Timing:** The `projectSupersites` function is only called once, at the very end of the program in `phaser_finalise.cpp`, just before writing the output file. It is **never called during the MCMC iterations**. This is the primary reason for the reference panel corruption.
2.  **Buggy Implementation:** The existing implementation of `projectSupersites` in `genotype_managment.cpp` is itself flawed. It contains the comment `// For now, assume only ALT1 mapping` and is hardcoded to only handle the first alternate allele of a supersite. It cannot correctly project phasing for any site with more than two alternate alleles.

### 10.4. The Plan for Resolution

To fix this critical bug and restore supersite phasing accuracy, the following two steps must be taken:

1.  **Fix the `projectSupersites` function:** The logic in `genotype_managment.cpp` must be rewritten to correctly handle an arbitrary number of alternate alleles. The corrected version should use the `getSampleSuperSiteAlleleCode` helper function to determine the true, phased allele class for each haplotype and then iterate through the constituent variants to set their biallelic state correctly.
2.  **Call `projectSupersites` after sampling:** A call to `this->projectSupersites()` must be added at the end of both the `genotype::sampleForward` and `genotype_sampleBackward` functions in `genotype_sweep.cpp`. This will ensure that the projection occurs immediately after a new haplotype is sampled and before the `updateHaplotypes` function is called, thus preventing reference panel corruption.

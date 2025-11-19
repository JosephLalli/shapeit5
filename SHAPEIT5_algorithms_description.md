# A Deep Dive into the SHAPEIT5 Phasing Algorithms

This document provides a comprehensive technical description of the core phasing algorithms within the SHAPEIT5 software suite. It is built upon an analysis of the C++ codebase from the `where-are-the-errors` branch and transcripts of conversations with the development team. We will explore the data structures, mathematical models, and computational workflows for both the standard **biallelic** phasing pathway and the more complex **supersite** pathway for handling multi-allelic variants.

## Core Concepts: Phasing, Imputation, and the Li & Stephens HMM

At its heart, SHAPEIT5 is a tool for resolving haplotype phase and imputing missing genotypes. 

*   **Haplotype Phasing** is the process of determining which of the two parental chromosomes each allele at a heterozygous site belongs to. For example, if an individual is genotyped as `A/T` at a specific locus, phasing determines whether their two chromosomes (haplotypes) carry `A` and `T` respectively, or `T` and `A`. The output is two complete haplotype sequences for each individual.
*   **Imputation** is the process of inferring genotypes at sites that were not directly measured, using the statistical patterns of linkage disequilibrium (LD) found in a reference panel of known haplotypes.

SHAPEIT5 accomplishes this using a **Markov Chain Monte Carlo (MCMC)** framework. The core of this framework is an iterative algorithm that repeatedly samples haplotype configurations for each individual, conditioned on a reference panel of other haplotypes. Over many iterations, this process converges to a statistically likely set of phased haplotypes.

The engine driving the conditional sampling is a version of the **Li & Stephens (2003) Hidden Markov Model (HMM)**. In this model, a target individual's haplotype is conceptualized as an imperfect mosaic of haplotypes drawn from a reference panel. The HMM moves along the chromosome from one variant to the next. At each step, it can either:
1.  **Copy:** Continue copying the allele from the same reference haplotype it was previously copying.
2.  **Switch/Recombine:** Switch to copying from a different, randomly chosen haplotype in the reference panel.

The HMM also allows for **mutations**, where the allele on the target haplotype does not match the copied reference allele, which accounts for genotyping errors or rare variation. By calculating the probabilities of these events in a forward-backward manner, the algorithm can determine the posterior probability of any given phase configuration.

A key optimization in SHAPEIT5 is the use of the **Positional Burrows-Wheeler Transform (PBWT)**, as proposed by Richard Durbin (2014). Instead of running the HMM against the entire, massive reference panel (which could contain millions of haplotypes), the PBWT is used in each MCMC iteration to swiftly select a small subset of `K` (typically ~200) "conditioning haplotypes" from the panel that are most similar to the target individual's current haplotype estimates in a local window. This drastically reduces the state space of the HMM, making the problem computationally tractable without a significant loss in accuracy.

Furthermore, the HMM is accelerated through the use of **Single Instruction, Multiple Data (SIMD)** instructions, specifically Intel's AVX2 extensions. This allows the algorithm to process eight distinct phasing hypotheses (called "lanes") simultaneously, providing a substantial performance boost.

With this foundational understanding, we can now dissect the two primary algorithmic pathways in SHAPEIT5.

---

## 1. The Biallelic Algorithm

The biallelic algorithm is the fundamental phasing pathway in SHAPEIT5, designed to operate on genetic variants that have exactly two alleles: a reference allele (REF, coded as `0`) and a single alternate allele (ALT, coded as `1`).

### 1.1. Abstract Description

The abstract goal of the biallelic algorithm is to resolve the parental origin of alleles at all heterozygous sites for a given individual. It is an iterative refinement process that seeks the most probable phase configuration consistent with the observed genotypes and the patterns of variation in a reference panel.

Each MCMC iteration involves updating the phase for every individual in the study set. For a single individual, this update proceeds as follows:

1.  **State Selection:** The algorithm first identifies the `K` most relevant haplotypes from the global reference panel to act as "donors" or templates for the HMM. This is achieved by the PBWT, which finds haplotypes that share long stretches of identity with the individual's current haplotypes, effectively creating a bespoke, localized reference panel for the HMM calculations.

2.  **HMM Computation:** The chromosome is partitioned into windows, and within each window, a forward-backward HMM is executed.
    *   **Forward Pass:** The algorithm calculates the `Alpha` probabilities: the probability of observing the individual's genotype data up to a given variant, while being in each of the `K` possible HMM states (i.e., copying from each of the `K` donor haplotypes).
    *   **Backward Pass:** The algorithm calculates the `Beta` probabilities: the probability of observing the genotype data *after* a given variant, given that the HMM is in a particular state at that variant.

3.  **Posterior Decoding and Sampling:** By combining the `Alpha` and `Beta` probabilities at each variant, the algorithm can compute the posterior probability of any specific phase configuration. It then samples a new, complete phase for the individual's two haplotypes from this posterior distribution. This involves sampling a "path" through the graph of possible segment diplotypes.

4.  **Haplotype Update:** This is the final, crucial step of the MCMC iteration where the newly sampled phase is committed back to the main dataset. In the code, this is handled by the `conditioning_set` object (typically named `H`), which acts as the "global reference panel". After the `genotype::make()` method updates an individual's `Variants` array, the main `phaser::phase()` loop calls `H.updateHaplotypes(G)`, where `G` is the `genotype_set` containing all individuals. This method iterates through each individual and copies their two newly-phased haplotypes from their local `Variants` array into the main `H.H_opt_hap` matrix. This overwrite ensures that in the next MCMC iteration, other individuals will use this updated phasing information during the PBWT state selection, allowing the algorithm to progressively converge on a globally consistent solution.

This entire process is repeated for a specified number of "burn-in" iterations (to allow the MCMC to reach a stationary distribution), "pruning" iterations (where confident segments are merged to simplify the genotype graph), and "main" iterations (where final probabilities are collected).

### 1.2. Functional Deep Dive

The biallelic algorithm is realized through a sophisticated interplay of several key C++ objects and functions, primarily within the `phase_common` module.

#### Core Data Structures

*   **`genotype` object (`objects/genotype/genotype_header.h`):** This is one of the most important and highly optimized data structures, representing a single individual's genotype data. It is not a simple linear array but is structured as a graph of segments to optimize the HMM.
    *   **`Variants`:** A `std::vector<unsigned char>` that stores the core allele information in a bit-packed format. Its encoding is fundamental to the algorithm's memory efficiency and performance.
        *   **Bit-Packed Encoding:** Each byte in this vector stores the state for two adjacent variants. Each variant is allocated 4 bits (a nibble). The layout for a single variant's nibble is as follows (from most to least significant bit):

| Bit | Name | Purpose | Values |
| :-- | :--- | :--- | :--- |
| **3** | `HAP1` | Allele on Haplotype 1 | `0` = REF, `1` = ALT |
| **2** | `HAP0` | Allele on Haplotype 0 | `0` = REF, `1` = ALT |
| **1** | State Bit 1 | Genotype State | (See below) |
| **0** | State Bit 0 | Genotype State | (See below) |

        *   The two "State Bits" combine to define the genotype's nature:
            | State Bits (`B1B0`) | Code | Meaning |
            | :--- | :--- | :--- |
            | `00` | `HOM` | **Homozygous**. Both parental chromosomes carry the same allele. The `HAP0` and `HAP1` bits will be identical. |
            | `01` | `MIS` | **Missing**. The genotype was not measured. The `HAP0`/`HAP1` bits are used to store the imputed (guessed) alleles. |
            | `10` | `HET` | **Heterozygous**. The individual has one REF and one ALT allele. This is an ambiguous site that requires phasing. |
            | `11` | `SCA` | **Scaffold**. A heterozygous site that has been confidently pre-phased using pedigree information. |
    *   **Segments & `Lengths`:** The genotype is partitioned into **segments** at any site that is confidently homozygous (`HOM`). The HMM only needs to run on the ambiguous regions between these anchors. The `Lengths` member, a `std::vector<unsigned short>`, stores the number of variants in each of these segments.
    *   **`Diplotypes`:** A `std::vector<unsigned long>` where each element is a 64-bit mask corresponding to a segment. This mask represents the set of possible phasing solutions for that segment. Since the HMM uses 8 parallel lanes (hypotheses), there are 8 possible phase states for the first haplotype and 8 for the second, yielding 8x8 = 64 possible **diplotypes**. A `1` at a bit position `d` in the mask means that diplotype `d` is a valid hypothesis for that segment.
    *   **`Ambiguous`:** A `std::vector<unsigned char>` storing one 8-bit mask for each heterozygous site. Each bit in the mask corresponds to one of the 8 SIMD lanes. The bit's value (`0` or `1`) defines the phasing hypothesis for that lane, determining whether it tests the `REF|ALT` or `ALT|REF` orientation.

*   **`conditioning_set` object:** Located in `containers/conditioning_set/`, this object manages the reference panel.
    *   `H_opt_hap`: A `bitmatrix` storing the full reference panel, with haplotypes as rows and variants as columns.
        *   **The `bitmatrix` class (`containers/bitmatrix.h`):** This is a custom data structure designed for extreme memory efficiency when storing binary data. Instead of using a `bool` or `char` for each entry (which would consume at least one byte per bit), it packs 8 binary values into a single `unsigned char`. This reduces the memory footprint of the reference panel by a factor of 8, which is critical for handling large datasets with millions of haplotypes. Accessing an element at `(row, col)` involves bitwise arithmetic: the `get()` and `set()` methods calculate which byte in the underlying 1D array contains the target bit and then use bitwise shifts (`>>`) and masks (`&`, `|`) to read or modify that specific bit without disturbing the other 7 bits in the byte.
    *   `H_opt_var`: The transpose of `H_opt_hap`, used for efficient column-wise (per-variant) access.
    *   `select()`: This is the high-level method in `conditioning_set_selection.cpp` that orchestrates the selection of the `K` conditioning states (or "neighbors") for each sample haplotype, which is the fundamental optimization enabling the Li & Stephens HMM to run efficiently. It uses the Positional Burrows-Wheeler Transform (PBWT) to achieve this.
        *   **Conceptual Goal:** The HMM's complexity is determined by its number of hidden states. A naive implementation would use all `N` haplotypes in the reference panel as states, which is computationally infeasible for large panels. The goal of `select()` is to intelligently choose a small, fixed-size subset of `K` haplotypes (where `K` is much smaller than `N`) that are most likely to be good templates for the individual's own haplotypes in the current genomic window.
        *   **The PBWT Algorithm:** The PBWT is an algorithm that efficiently sorts a collection of binary sequences. As it moves from one position (variant) to the next, it maintains a sorted order of the haplotypes based on their sequence identity up to that point. A key output is the "divergence array," which records, for each haplotype in the sorted list, the position at which it last differed from the haplotype preceding it in the list. Haplotypes that are identical over long stretches will remain adjacent in the sorted list and have high divergence values.
        *   **Functional Implementation:**
            1.  The top-level `select()` method first stochastically chooses a set of loci (`sites_pbwt_selection`) at which the neighbor selection will actually occur. This is done by grouping candidate loci and randomly picking one from each group.
            2.  It then dispatches the core work to one or more threads, each executing the `select(int chunk)` method on a different portion of the chromosome.
            3.  The `select(int chunk)` method implements the PBWT. It iterates from locus to locus, maintaining the permutation array `A` (the sorted order of haplotype indices) and the divergence array `C`.
            4.  When the loop reaches a pre-selected site `l`, it calls `store(l, A, C)`.
            5.  The `store()` function performs the actual neighbor selection. For each of the individual's own haplotypes, it finds its current rank `h` in the sorted array `A`. It then "looks up" and "looks down" from `h` in the array, collecting the `depth` nearest neighbors in each direction. These neighbors are good candidates because their proximity in the sorted list implies they share a long prefix with the target haplotype.
            6.  The collected neighbors for both of an individual's haplotypes are pooled, and the final set of `K` states is stored in the `indexes_pbwt_neighbour` vector. This vector is initially laid out in memory for efficient writing during the PBWT pass.
            7.  Finally, `transposePBWTneighbours()` is called to reorder the `indexes_pbwt_neighbour` vector into a layout that is optimized for fast, contiguous memory access by the HMM during the subsequent `phaseWindow` computations.

#### Algorithmic Workflow and Key Functions

1.  **Orchestration (`phaser::phase` in `phaser/phaser_algorithm.cpp`):**
    The top-level `phase()` method orchestrates the entire MCMC process. It loops through the specified number of burn-in, pruning, and main iterations. In each iteration, it first calls `H.select()` to choose the conditioning states for all individuals. It then calls `phaseWindow()` to perform the HMM calculations and sampling. Finally, it calls `H.updateHaplotypes(G)` to refresh the reference panel with the newly sampled phases, followed by `H.transposeHaplotypes_H2V(false)` to prepare the data layout for the next PBWT selection.

2.  **HMM Computation (`haplotype_segment_single` in `models/haplotype_segment_single.h`):**
    This class is where the core forward-backward algorithm is implemented. An instance is created for each window of a single individual's chromosome.
    *   **`forward()`:** This method implements the forward pass. It iterates from the beginning to the end of a window. At each variant, it executes one of three main routines based on the genotype state (`HOM`, `HET`, `MIS`). For example, at a homozygous site, `RUN_HOM` is called. This routine updates the `prob` array, which holds the probability of each of the `K` states and 8 lanes. The update rule is a SIMD-optimized version of the Li & Stephens model:
        ```cpp
        // Simplified concept from RUN_HOM
        __m256 _prob = _mm256_load_ps(&prob[i]);
        // Update combines probability of staying with the same donor (multiplied by nt)
        // and switching to a new donor (adding the weighted average tFreq).
        _prob = _mm256_fmadd_ps(_prob, _nt, _tFreq); // prob = (prob * nt) + tFreq
        // Apply a penalty if the donor allele does not match the observed genotype.
        if (sample_allele != donor_allele) _prob = _mm256_mul_ps(_prob, _mismatch_penalty);
        _mm256_store_ps(&prob[i], _prob);
        ```
    *   At the boundary between two segments (i.e., after a run of heterozygous sites), the forward probabilities are stored in the `Alpha` vector. This is a critical step for connecting the HMM calculations across the chromosome's segments. The `Alpha` vector is a temporary data structure created within the `haplotype_segment` object for each HMM window. It acts as the memory of the forward pass. Structurally, it is a vector of vectors, effectively a 2D array `Alpha[segment_index][probability_index]`, where `segment_index` refers to the segments within the current window. The inner vector is a flat array of size `K * 8`, storing the raw forward probability for each of the `K` conditioning states and each of the 8 SIMD lanes. When the `forward()` pass finishes a segment, it copies its final probability array into the `Alpha` vector. The `backward()` pass then retrieves these stored values to correctly initialize its calculations at each segment boundary, allowing the algorithm to compute transition probabilities between segments and form a single, unbroken probabilistic chain across the entire window.
    *   **`backward()`:** This method works in reverse, from the end of the window to the beginning. It performs a similar set of calculations to compute `Beta` probabilities. As it proceeds, it combines the stored `Alpha` probabilities with the live `Beta` probabilities to compute the full posterior probabilities for transitions between the 64 possible diplotypes at each segment boundary. These transition probabilities are stored in the `CurrentTransProbabilities` vector passed in from the `compute_job`. For missing sites, it calls `IMPUTE()` to calculate the posterior probability of the alternate allele, which is stored in `CurrentMissingProbabilities`.

3.  **Sampling and State Update (`genotype` methods):**
    *   **`sampleForward()` / `sampleBackward()` (`objects/genotype/genotype_sweep.cpp`):** After the `backward()` pass has populated `CurrentTransProbabilities`, one of these methods is called. It walks through the segments of the genotype graph and uses the computed transition probabilities to perform a stochastic backtrack, sampling a single, most likely path of diplotypes through the segments. The result is stored in the `DipSampled` vector.
    *   **`make()` (`objects/genotype/genotype_managment.cpp`):** This is the final step where the sampled path is physically realized. For each segment, it takes the sampled diplotype from `DipSampled`, which specifies a pair of lanes (e.g., `hap0_lane`, `hap1_lane`). For each heterozygous site within that segment, it looks up the 8-bit `Ambiguous` mask. It uses the mask to determine the phase orientation implied by the chosen lanes (e.g., should `hap0` get REF and `hap1` get ALT, or vice versa?). It then uses the `VAR_SET_HAP0` and `VAR_SET_HAP1` macros to write these two bits directly into the correct 4-bit nibble within the `Variants` vector, overwriting the previous phase.

This detailed cycle of selection, HMM computation, sampling, and state-writing forms one complete MCMC update for the biallelic phasing algorithm.

---

## 2. The Supersite Algorithm

The supersite algorithm is a sophisticated extension designed to correctly handle multi-allelic variants, which are common in modern sequencing data but are represented as multiple biallelic records in VCF/BCF files. The primary challenge is ensuring that, for a single biological site, a phased haplotype contains at most one of the possible alternate alleles.

### 2.1. Abstract Description

The supersite algorithm mirrors the structure of the biallelic pathway but introduces a layer of abstraction to manage the complexity of multiple alleles.

1.  **Supersite Creation:** As a pre-processing step, the algorithm scans the variant list and groups all records that share the same chromosome and position into a `SuperSite` object. Within this group, one variant is designated the **anchor** (typically the one with the highest minor allele frequency), which will serve as the site's proxy during HMM calculations. All other variants at that position are its **siblings**.

2.  **State Selection (Anchor-Gated):** The PBWT state selection is modified to be "anchor-aware." When selecting the `K` conditioning states, it effectively ignores the sibling variants and only considers the genotypes at the anchor loci. This ensures the HMM's donor haplotypes are chosen based on the overall state of the multi-allelic site.

3.  **HMM Computation (Class-Based):** The HMM proceeds as normal, but its behavior is altered at supersites. The dynamic programming calculations are only performed at the anchor locus; sibling loci are treated as no-ops. The key innovation is the use of **4-bit class codes** instead of binary alleles. For a supersite with `N` alternate alleles, class `0` is REF, class `1` is the first ALT, class `2` is the second ALT, and so on.
    *   **Immutable State (`c0`/`c1`):** For a heterozygous supersite in a given sample (e.g., genotype `A/G` where `A` is `ALT1` and `G` is `ALT2`), the algorithm identifies the two canonical allele classes present. These are stored as `c0` and `c1` (e.g., `c0=1`, `c1=2`). This pair is considered the **immutable biological state** of the site for the sample and is used for all emission calculations within the HMM. This insight, derived from the developer conversations, was critical to fixing bugs related to reference panel corruption. The standard `Ambiguous` mask is reused to test the `c0|c1` vs `c1|c0` orientations across the 8 SIMD lanes.

4.  **Sampling (Ephemeral Choice `h0`/`h1`):** The sampling process remains structurally the same. A path is sampled through the segment graph. For a heterozygous supersite, this results in the selection of an orientation, which maps the chosen lanes to the immutable `c0`/`c1` classes. This choice for the *current iteration* is stored as an **ephemeral (`h0`, `h1`) class pair**.

5.  **Projection:** This is the most critical step unique to the supersite algorithm. After `make()` determines the ephemeral `h0` and `h1` classes, the `projectSupersites()` function is called. This function iterates over all sibling variants within the supersite. For the first haplotype, it sets the `HAP0` bit to `1` *only* for the specific variant (which could be the anchor or a sibling) corresponding to the `h0` class and clears the `HAP0` bit for all other variants within that supersite. It does the same for the second haplotype with `h1`. This guarantees that the final bit-packed representation in the `Variants` array is always biologically valid, with at most one alternate allele per haplotype.

6.  **Haplotype Update:** After the `projectSupersites()` function ensures the individual's `Variants` array has a biologically correct representation of the multi-allelic site, the update proceeds identically to the biallelic case. The main algorithm calls `H.updateHaplotypes(G)`, where the "global reference panel" (`H`, the `conditioning_set` object) is updated. This method copies the newly-phased and correctly projected haplotype data from the individual's `Variants` array into the central `H.H_opt_hap` matrix. This makes the improved, consistent haplotype available as a template for all other individuals in the next MCMC iteration.

### 2.2. Functional Deep Dive

The supersite pathway cleverly reuses the biallelic infrastructure by introducing a carefully designed abstraction layer.

#### Core Data Structures

*   **`SuperSite` object (`objects/super_site_accessor.h`):** A simple struct that defines a multi-allelic site. Key fields include:
    *   `global_site_id`: The VCF index of the anchor variant.
    *   `var_start`, `var_count`: The range of indices in the global `variant_map` that belong to this supersite.
    *   `panel_offset`: An offset used to locate this site's data within the `packed_allele_codes` vector.

*   **Global Supersite Maps:**
    *   `locus_to_super_idx`: A `std::vector<int>` that maps every variant locus in the dataset to its corresponding `SuperSite` index, or -1 if it's not part of a supersite. This is the primary mechanism for detecting supersites.
    *   `packed_allele_codes`: A `std::vector<uint8_t>` where the allele classes for all haplotypes at all supersites are stored in a compact format (2 classes per byte). `buildSuperSites` (in `objects/super_site_builder.cpp`) is responsible for generating this.

*   **`genotype` object extensions:**
    *   `supersite_class_pairs_base`: A `std::vector<uint8_t>` that stores the immutable `c0`/`c1` class pairs for each supersite in the sample. This is populated by `snapshotSupersiteBaseClasses()` before the HMM begins and is read by the HMM emissions logic.
    *   `supersite_class_pairs`: A `std::vector<uint8_t>` that stores the mutable, per-iteration sampled `h0`/`h1` classes. This is updated during sampling by `setSupersiteClassPair()`.
    *   `setSuperSiteContext()`: A method to attach pointers to the global supersite data structures (`super_sites`, `locus_to_super_idx`, etc.) to the `genotype` object, making them accessible during HMM calculations.

#### Supersite Emission Calculation: The "Match" Logic

A critical difference in the supersite pathway is how the HMM determines if a conditioning haplotype "matches" the sample's genotype. Unlike the simple allele comparison in the biallelic case, this is a sophisticated, real-time probabilistic process designed for maximum performance using SIMD instructions. The calculation happens inside the HMM forward pass, primarily within the `SS_RUN_AMB` and `SS_RUN_HOM` methods in `phase_common/src/models/haplotype_segment_single.h`.

The logic does **not** use a large, pre-computed table of emission probabilities as previously described. Instead, it calculates these probabilities on the fly.

**Part 1: Caching Conditioning Haplotype Codes**

As a minor optimization, the first time a supersite is encountered within an HMM window, a helper function `ss_load_cond_codes` is called. This function populates a temporary vector, `ss_cond_codes`, which serves as a local cache for the duration of the HMM window.

The process is as follows:
1.  The function iterates from `k = 0` to `K-1` for the `K` conditioning haplotypes selected for the current sample.
2.  For each `k`, it retrieves the global index of the conditioning haplotype.
3.  It then calls `unpackSuperSiteCode`, which uses the haplotype's global index and the supersite's specific offset (`ss.panel_offset`) to look up and unpack the correct 4-bit allele code from the global `panel_codes` array.
4.  The resulting 4-bit code is stored in `ss_cond_codes[k]`.

This creates a small, aligned cache containing only the allele codes for the `K` haplotypes relevant to the current HMM calculation, avoiding repeated lookups into the much larger global `panel_codes` array.

**Part 2: Real-time, SIMD-based Emission Calculation**

This is the core of the logic, executed at every supersite anchor for each of the `K` conditioning haplotypes. The most illustrative case is for a heterozygous site in `SS_RUN_AMB`.

1.  **Determine Per-Lane Expectations:** The algorithm first determines what allele each of the 8 SIMD lanes is testing. It retrieves the sample's two allele classes for the site (e.g., `c0` and `c1`) and uses the `G->Ambiguous` mask to build an `expected_class` array. For example, `expected_class[0]` might be `c0` while `expected_class[1]` is `c1`, and so on for all 8 lanes.

2.  **Perform SIMD Comparison:** For a given conditioning haplotype `k`, the code gets its 4-bit allele code from the cache: `ss_cond_codes[k]`. This value represents the allele class (0 for REF, 1 for ALT1, etc.) of that specific donor haplotype. It then uses AVX2 intrinsics to compare this single donor code against all 8 lane expectations simultaneously.
    *   The 8 expected-class codes are loaded into a 256-bit SIMD register (`exp_vec`).
    *   The single donor code, `ss_cond_codes[k]`, is broadcast into all 8 slots of another SIMD register (`donor_vec`).
    *   A single comparison instruction, `_mm256_cmpeq_epi32(donor_vec, exp_vec)`, is executed. This produces a result mask where the elements corresponding to lanes with a match are set to all `1`s, and mismatches are set to all `0`s.

3.  **Select Match/Mismatch Probability:** The final emission probabilities are generated with the `_mm256_blendv_ps` instruction. This powerful "blend" instruction takes three inputs: a vector of "mismatch" probabilities (all lanes set to `M.ed/M.ee`), a vector of "match" probabilities (all lanes set to `1.0f`), and the result mask from the comparison step. In a single operation, it builds the final `emit` register, selecting the match or mismatch value for each lane based on the mask.

This `emit` register, containing the 8 correct and lane-specific emission probabilities, is then multiplied by the forward-probability (`prob`) register to complete the HMM step. This real-time, SIMD-based approach provides extreme performance without the need for a large, pre-computed lookup table.

#### Algorithmic Workflow and Key Functions

1.  **Initialization (`phaser::initialise`):** Before the MCMC iterations start, the `buildSuperSites` function (in `objects/super_site_builder.cpp`) is called. This function is responsible for two critical setup tasks:
    *   **Supersite Discovery:** It analyzes the `variant_map`, groups variants by position, and creates the `SuperSite` objects and associated lookup maps (`locus_to_super_idx`).
    *   **Global Allele Code Generation:** For each supersite, it iterates through **every haplotype in the entire reference panel** (`H.n_hap`). For each haplotype, it determines which allele class it carries (0 for REF, 1 for the first ALT, 2 for the second, etc.) by checking the 1-bit values of the constituent biallelic variants. The resulting 4-bit allele code is then stored in the global `packed_allele_codes` vector. This vector, later accessed via the `panel_codes` pointer, acts as a master database of allele information for all supersites.
    *   Finally, `setSuperSiteContext()` is called for every `genotype`, and `snapshotSupersiteBaseClasses()` is called to capture the initial immutable `c0`/`c1` state for every heterozygous supersite in each sample.

2.  **The Abstraction Layer (`SiteEmissionAdapter` in `models/site_emission_adapter.h`):**
    This class is the linchpin that allows the `haplotype_segment_single` HMM core to remain agnostic about supersites. The main HMM loop doesn't contain `if (is_supersite)` logic. Instead, it calls `SupersiteEmissionAdapter::build_view()`.
    *   `build_view()` checks `locus_to_super_idx[abs_locus]`.
    *   If it's a normal site, it populates a `SiteView` object just like the biallelic case.
    *   If it's a supersite sibling, it sets `view.kind = SiteKind::SuperSibling`, telling the HMM to do nothing.
    *   If it's a supersite anchor, it sets `view.kind = SiteKind::SuperAnchor`. Crucially, it calls `G_->getSupersiteBaseClassPair()` to retrieve the **immutable `c0`/`c1`** and uses these, along with the standard `Ambiguous` mask, to determine the `view.lane_class` array of 4-bit codes for the 8 SIMD lanes.

3.  **HMM Computation (`haplotype_segment_single`):**
    The main `forward()` loop sees `view.kind == SiteKind::SuperAnchor` and calls `SS_RUN_AMB` or `SS_RUN_HOM`.
    *   `SS_RUN_AMB`: This function is the supersite equivalent of `RUN_AMB`. Instead of a simple `0/1` comparison, it compares the 4-bit donor allele code (from the `ss_cond_codes` cache) against the 4-bit expected class code in `view.lane_class[h]` for each lane `h`. The SIMD logic is more complex, often using `_mm256_blendv_ps` to select between match and mismatch penalties based on a comparison vector.

4.  **Sampling (`genotype::make`):**
    The logic for sampling a diplotype path remains the same. However, the `make()` function has a special branch for supersite anchors.
    *   When it encounters an ambiguous supersite anchor, it retrieves the immutable `c0`/`c1` from `supersite_class_pairs_base`.
    *   It then uses the sampled lanes (`hap0_lane`, `hap1_lane`) and the `Ambiguous` mask to decide the orientation, producing the ephemeral `h0` and `h1` classes for the current iteration.
    *   It stores this ephemeral choice by calling `setSupersiteClassPair(ss_idx, h0, h1)`.
    *   No bits are written to the `Variants` array at this stage. The function simply records the choice and moves on.

5.  **Projection (`genotype::projectSupersites()`):**
    This function is called from `sampleForward()` and `sampleBackward()` *immediately after* `make()`. This was a critical design decision to ensure consistency.
    *   It iterates through all supersites in the genotype.
    *   For each one, it retrieves the ephemeral `h0`/`h1` choice that was just stored by `make()` via `getSupersiteClassPair()`.
    *   It then loops from `0` to `ss.var_count - 1`. For each sibling variant in the supersite, it determines its own alternate allele class (e.g., sibling `i` corresponds to class `i+1`).
    *   It sets the `HAP0` bit for that sibling to `1` only if `h0` equals the sibling's class; otherwise, it clears it to `0`.
    *   It does the same for the `HAP1` bit with `h1`.
    This ensures that, for instance, if `h0` was class `2`, only the second sibling variant will have its `HAP0` bit set to `1`, and all other siblings will have their `HAP0` bit set to `0`, perfectly enforcing the biological constraint of a multi-allelic site.

Through this elegant separation of immutable biological state (`c0`/`c1`) from the ephemeral per-iteration choices (`h0`/`h1`) and a final, definitive projection step, the supersite algorithm successfully extends the powerful, highly optimized biallelic HMM framework to handle the full complexity of modern sequencing data.
<ctrl95>Here is a comprehensive description of the biallelic and supersite algorithms used in SHAPEIT5, based on the provided conversation logs and an analysis of the codebase.

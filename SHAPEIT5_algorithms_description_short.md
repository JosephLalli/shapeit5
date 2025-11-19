## The SHAPEIT5 Phasing Algorithm

SHAPEIT5 is a haplotype phasing and imputation tool that uses a Markov Chain Monte Carlo (MCMC) approach. The core of the algorithm is an iterative process that refines the phase of heterozygous genotypes by sampling from a conditional distribution of haplotypes. This process is performed for each individual in a dataset, using a reference panel of haplotypes to inform the phasing decisions.

The algorithm can be broken down into two main pathways: a **biallelic pathway** for sites with two alleles (a reference and a single alternate), and a **supersite pathway** for sites with more than two alleles (multiallelic sites). The supersite pathway is designed to be a direct, more complex, analogue of the biallelic pathway.

---

## 1. The Biallelic Algorithm

The biallelic algorithm is the fundamental phasing pathway in SHAPEIT5. It processes each variant with exactly two alleles independently.

### 1.1. Abstract Description

The goal of the biallelic algorithm is to determine, for each heterozygous site in a sample's genotype, which allele belongs to which of the two parental haplotypes. It does this iteratively. In each iteration, for each individual, it performs the following steps:

1.  **State Selection:** It efficiently selects a small subset of `K` "conditioning" haplotypes from a large reference panel. These `K` haplotypes are chosen because they are the most similar to the individual's current estimated haplotypes. This is done using the Positional Burrows-Wheeler Transform (PBWT), which dramatically reduces the computational complexity of the subsequent steps.

2.  **HMM Computation:** It uses a forward-backward Hidden Markov Model (HMM) to calculate the posterior probability of every possible phase configuration for the individual's genotype, given the `K` selected conditioning haplotypes. The HMM moves along the chromosome, and at each variant, it considers the probability of "copying" an allele from one of the `K` conditioning haplotypes versus "mutating" to a different allele.

3.  **Sampling:** Based on the posterior probabilities computed by the HMM, a new phase configuration for the individual's haplotypes is randomly sampled. This means that for each heterozygous site, a decision is made about which allele goes on which haplotype.

4.  **Haplotype Update:** The individual's haplotypes in the reference panel are updated with the newly sampled phase configuration.

This process (Steps 1-4) is repeated for a set number of MCMC iterations (burn-in, pruning, and main stages), allowing the phase estimates to converge towards the most likely configuration.

### 1.2. Functional Description

The biallelic algorithm is implemented through the interaction of several key C++ objects and functions:

*   **`phaser::phase()`:** This is the main orchestrating function in `phaser/phaser_algorithm.cpp`. It runs the main MCMC loop, which consists of burn-in, pruning, and main iterations. In each iteration, it calls `H.select()` and then `phaseWindow()`.

*   **`conditioning_set` object:** This object manages the full reference panel and the PBWT-selected conditioning states.
    *   `H_opt_hap`: A `bitmatrix` storing the main reference panel, with haplotypes as rows and variants as columns.
    *   `select()`: This method in `containers/conditioning_set/conditioning_set_selection.cpp` implements the PBWT algorithm to select the `K` conditioning haplotypes for each sample. The results (indices of the selected haplotypes) are stored in `indexes_pbwt_neighbour`.
    *   `updateHaplotypes()`: After sampling, this method updates the `H_opt_hap` matrix with the new haplotype data from the `genotype::Variants` array of each sample.

*   **`phaseWindow()`:** This function, within `phaser/phaser_algorithm.cpp`, dispatches the HMM computation work to multiple threads. Each thread processes a subset of the individuals, creating a `compute_job` for each.

*   **`genotype` object:** This object in `objects/genotype/genotype_header.h` is a complex and highly optimized representation of a single individual's genotype data. It is structured as a graph of segments.
    *   **Segments:** The genotype is not treated as a single monolithic sequence. Instead, it is broken down into **segments** at homozygous sites. The HMM runs independently between these homozygous sites. This is a key optimization.
    *   `Lengths`: A `std::vector<unsigned short>` that stores the number of variants in each segment.
    *   `Variants`: A bit-packed `std::vector<unsigned char>`. This is the core data store for the allele information.
        *   **Bit-Packed Encoding:** Each byte in this vector stores the state for two adjacent variants. Each variant gets 4 bits (a nibble) to store its information. The layout for a single variant's 4 bits is as follows (from most to least significant):
            | Bit | Name | Purpose | Values |
            | :-- | :--- | :--- | :--- |
            | **3** | `HAP1` | Allele on Haplotype 1 | `0` = REF, `1` = ALT |
            | **2** | `HAP0` | Allele on Haplotype 0 | `0` = REF, `1` = ALT |
            | **1** | State Bit 1 | Genotype State | (See below) |
            | **0** | State Bit 0 | Genotype State | (See below) |
        *   The two "State Bits" work together to define the nature of the genotype:
            | State Bits (`B1B0`) | Code | Meaning |
            | :--- | :--- | :--- |
            | `00` | `HOM` | **Homozygous**. The `HAP0` and `HAP1` bits will be identical. |
            | `01` | `MIS` | **Missing**. The genotype is unknown. |
            | `10` | `HET` | **Heterozygous**. The site is ambiguous and requires phasing. |
            | `11` | `SCA` | **Scaffold**. A heterozygous site phased with high confidence by pedigree data. |
    *   `Ambiguous`: A `std::vector<unsigned char>` storing an 8-bit mask for each heterozygous (`HET` or `SCA`) site. Each bit corresponds to one of the 8 SIMD lanes in the HMM. The bit's value (`0` or `1`) defines the phasing hypothesis for that lane (e.g., `0` = REF/ALT, `1` = ALT/REF).
    *   `Diplotypes`: A `std::vector<unsigned long>` storing a 64-bit mask for each segment. Each of the 64 bits represents a possible **diplotype**, which is a specific pairing of the 8 haplotype hypotheses for haplotype 0 and the 8 for haplotype 1 (8x8=64). The HMM forward pass calculates probabilities for each of the 8 haplotype states, and the backward pass combines these to calculate the posterior probability of each of the 64 diplotypes.
    *   `sampleForward()` and `sampleBackward()`: These methods in `objects/genotype/genotype_sweep.cpp` initiate the sampling process. They traverse the segments and use the HMM's posterior probabilities (`CurrentTransProbabilities`) to sample a diplotype for each segment. They then call `genotype::make()`.
    *   `make()`: This function in `objects/genotype/genotype_managment.cpp` takes the sampled diplotype path. For each ambiguous site, it uses the `Ambiguous` orientation mask and the sampled lanes to write the new phase directly into the `Variants` byte array using the `VAR_SET_HAP0` and `VAR_SET_HAP1` macros.

*   **`haplotype_segment_single`:** This class in `models/haplotype_segment_single.h` is the computational core of the HMM. It is a transient object created for each window of computation.
    *   It performs the forward (`forward()`) and backward (`backward()`) passes.
    *   The core HMM logic, such as `RUN_HOM` and `RUN_AMB`, is heavily optimized with AVX2 SIMD intrinsics to process 8 haplotype phasing hypotheses (lanes) in parallel.
    *   It stores the DP matrix in `prob` (a vector of floats for `K` states x 8 lanes) and stores boundary information in the `Alpha` vectors to pass between segments.

---

## 2. The Supersite Algorithm

The supersite algorithm extends the biallelic model to handle multiallelic variants. The key challenge is to ensure that for a given sample, only one of the possible alternate alleles is chosen per haplotype, maintaining biological validity.

### 2.1. Abstract Description

The supersite algorithm builds upon the same MCMC framework as the biallelic one, but with crucial modifications.

1.  **Supersite Creation:** Before the MCMC iterations begin, the algorithm identifies all variant records at the same genomic position and groups them into a **`SuperSite`** object. One variant is designated the **anchor**, and the others are **siblings**.

2.  **State Selection (with Anchor Gating):** The PBWT state selection is "gated" to only consider the anchor variants of supersites, ensuring the conditioning set is built based on the overall state of the multiallelic site.

3.  **HMM Computation (with Class-based Emissions):** The HMM runs only at anchor loci. The core difference is the use of 4-bit **class codes** instead of binary alleles (e.g., 0=REF, 1=ALT1, 2=ALT2).
    *   For a heterozygous supersite, the two canonical allele classes for that sample are called **`c0`** and **`c1`**. This pair is treated as **immutable** during HMM emission calculations. The same `Ambiguous` mask from the biallelic path is used to orient `c0` and `c1` for the SIMD lanes.

4.  **Sampling (with Ephemeral Choices):** During sampling, the HMM's posteriors are used to choose an orientation. This results in an **ephemeral** choice for the current iteration, called **`h0`** and **`h1`**, which will be one of the two immutable classes (`c0` or `c1`).

5.  **Projection:** After sampling, a critical **projection** step occurs. The `projectSupersites()` function ensures that for each haplotype, only the single bit corresponding to the chosen class (`h0` or `h1`) is set across all sibling variants, while all others are cleared. This guarantees a biologically valid bit-packed representation.

6.  **Haplotype Update:** The main reference panel is updated with this projected, valid haplotype state.

The key design principle is the separation of **immutable state** (`c0`/`c1`) from **mutable, per-iteration choices** (`h0`/`h1`), making the supersite path a true analogue of the biallelic one.

### 2.2. Functional Description

The supersite algorithm modifies several key data structures and functions:

*   **`buildSuperSites()`:** (in `objects/super_site_builder.cpp`) Called during initialization, this function scans the `variant_map` and creates `SuperSite` objects, building the `locus_to_super_idx` map that links every variant locus to its supersite index (-1 if not part of one).

*   **`SuperSite` object:** A struct that defines a multiallelic site, containing its `global_site_id` (the anchor locus), `var_start` and `var_count` (to identify member variants), and `panel_offset` (for finding its data in a packed format).

*   **`genotype` object extensions:**
    *   `supersite_class_pairs_base`: A `std::vector<uint8_t>` storing the immutable `c0`/`c1` class pairs for each supersite. This is snapshotted before the HMM by `snapshotSupersiteBaseClasses()`.
    *   `supersite_class_pairs`: A `std::vector<uint8_t>` storing the mutable, per-iteration sampled classes `h0`/`h1`. This is updated by `setSupersiteClassPair()` during sampling.
    *   `setSuperSiteContext()`: Attaches pointers to the global supersite data structures to the `genotype` object.

*   **`SiteEmissionAdapter`:** (in `models/site_emission_adapter.h`) This class acts as an abstraction layer. The main HMM loop calls `build_view()` to get a `SiteView` object.
    *   `SupersiteEmissionAdapter::build_view()`: This is the crucial function.
        *   It checks `locus_to_super_idx` to identify the site type.
        *   For an anchor, it sets `view.kind = SiteKind::SuperAnchor` and retrieves the **immutable `c0`/`c1`** from `G_->getSupersiteBaseClassPair()`. The per-lane expected classes (`view.lane_class`) are then determined by orienting `c0` and `c1` with the `Ambiguous` mask.
        *   For a sibling, it sets `view.kind = SiteKind::SuperSibling`, causing the HMM to skip DP calculations.

*   **`haplotype_segment_single` modifications:**
    *   The main HMM loop dispatches based on the `SiteView.kind`. At anchor sites, it calls supersite-specific routines (`SS_RUN_AMB`).
    *   `SS_RUN_AMB`: Instead of a simple `0/1` comparison, it uses a match mask based on the 4-bit `lane_class` for each SIMD lane against the 4-bit donor codes.

*   **`genotype::make()` modifications:** (in `objects/genotype/genotype_managment.cpp`)
    *   For an ambiguous supersite anchor, it resolves the ephemeral `h0` and `h1` classes by applying the `Ambiguous` mask to the immutable `c0`/`c1`.
    *   It calls `setSupersiteClassPair()` to store the sampled `h0`/`h1` for the current epoch.
    *   Crucially, **it does not directly write to `Variants` here.** The writing is deferred to `projectSupersites()`.

*   **`genotype::projectSupersites()`:** (in `objects/genotype/genotype_managment.cpp`)
    *   This function is called immediately after `make()` in `sampleForward()` and `sampleBackward()`.
    *   It iterates through all supersites, retrieves the just-sampled `h0`/`h1` pair, and then iterates through all member splits. For each split, it sets the haplotype bit if its class (`h0` or `h1`) matches the split's alternate allele class, and clears it otherwise. This enforces the one-ALT-per-haplotype rule and is the supersite equivalent of `VAR_SET_HAP0`/`VAR_SET_HAP1`.
<ctrl95>Here is a comprehensive description of the biallelic and supersite algorithms used in SHAPEIT5, based on the provided conversation logs and an analysis of the codebase.

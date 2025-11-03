## Project
SHAPEIT5 — `phase_common` super-site integration

## Source of Truth
- This AGENTS.md mirrors `.github/copilot-instructions.md`, which is the authoritative, detailed guide distilled from `SUPERSITE_CONVERSATION_SUMMARY.md`.
- If anything diverges, prefer `.github/copilot-instructions.md` and update this file to match.

## Context
- SHAPEIT5's algorithm is a program that employs two variations (phase_common and phase_rare) of the MCMC Li-Stephens model to take genetic variant calls and determine the allele of origin (right or left) for heterozygous variants
- When variant callers were unsure of the genotype and emitted a missing call(./.), SHAPEIT5 also imputes the most likely genotype.
- It only supports biallelic sites, as it internally represents variants as 0 (reference) or 1 (not reference)
- MCMC uses forwards/backwards algorithm to calculate "emission" probability (probability of 1; default 0.001 likelihood of genotype error) and "transmission probability" (likelihood of variant arising from different source haplotype)

## SHAPEIT5 algorithm details

### What is Phasing?
Phasing is the process of determining which alleles came from which parent. When sequencing DNA, we get genotype calls like "0/1" at heterozygous sites, but we don't know if the "1" allele is on the chromosome inherited from mom or dad. Phasing figures this out, creating two separate haplotypes (one maternal, one paternal).

For example, if you have three nearby variants:
- Variant 1: A/G (genotype 0/1)
- Variant 2: C/T (genotype 0/1)  
- Variant 3: G/A (genotype 1/0)

Phasing determines whether the haplotypes are:
- Haplotype 1 (maternal): A-C-A
- Haplotype 2 (paternal): G-T-G

Or some other combination.

### The Li & Stephens Hidden Markov Model (HMM)

SHAPEIT5 uses a probabilistic model based on the Li & Stephens (2003) algorithm. The key insight is that **your DNA sequence is a mosaic of pieces copied from other individuals in the population**, with occasional copying errors and recombination events that switch which individual you're copying from.

Think of it like this: imagine you're writing a term paper by copying paragraphs from 100 different sources. Most of the time, you copy exactly, but sometimes you make a typo (sequencing error). And every so often, you switch to copying from a different source (recombination event). If someone gave you the final paper and the 100 sources, you could probably figure out which paragraphs came from which source - that's essentially what the HMM does.

**Key HMM Components:**

1. **Hidden States (K states)**: Each state represents "currently copying from donor haplotype K". SHAPEIT5 typically uses K≈200 donor haplotypes selected via PBWT (explained below).

2. **Emission Probabilities**: The probability of observing your allele given you're copying from donor K:
   - If your allele matches the donor: probability ≈ 1.0 (correct copy)
   - If your allele doesn't match: probability ≈ 0.001 (sequencing error rate, default `ed/ee`)
   - If your genotype is missing (./.): probability ≈ 1.0 for all donors (uninformative)

3. **Transition Probabilities**: The probability of switching donors between adjacent variants:
   - `nt` (no transition): probability of staying with the same donor = 1 - recombination_rate
   - `yt` (yes transition): probability of switching to a different donor = recombination_rate / (K-1)
   - Recombination rate is based on genetic distance (centiMorgans) from a genetic map

4. **Forward-Backward Algorithm**: 
   - **Forward pass**: Computes α(k) = probability of observing all variants up to position i and being in state k
   - **Backward pass**: Computes β(k) = probability of observing all variants after position i given state k
   - Together: posterior probability = α(k) × β(k), normalized

### SHAPEIT5 Modes

SHAPEIT5 has two phasing/imputation modes:
- **phase_common**: Used on common variants (MAF > 0.001 by default). Uses the full HMM on all samples.
- **phase_rare**: Used on rare variants. Phases rare variants onto existing common variant scaffolds using a simplified HMM.

This documentation focuses on **phase_common**, which is the core algorithm.

---

## phase_common Implementation Outline

### High-Level Architecture

The `phase_common` algorithm consists of several major subsystems that work together:

```
Input VCF/BCF
     ↓
[1] Data Loading & Initialization
     ↓
[2] MCMC Iteration Loop (burn-in → pruning → main)
     ├─→ [3] PBWT State Selection (select K donor haplotypes per sample)
     ├─→ [4] Window Segmentation (divide genome into manageable chunks)
     ├─→ [5] Multi-threaded HMM Computation
     │        ├─→ Forward Pass (compute α values)
     │        └─→ Backward Pass (compute β values + imputation)
     ├─→ [6] Sampling (pick most likely phase configuration)
     └─→ [7] Haplotype Update (update reference panel with new phases)
     ↓
Output Phased VCF/BCF
```

### [1] Data Loading & Initialization (`phaser_initialise.cpp`)

**Key Data Structures:**
- `genotype_set G`: All samples being phased
  - Each `genotype` object stores one sample's data as a graph of segments
  - Variants encoded as bit-packed values: HOM/HET/MIS/SCAFFOLD (2 bits per variant)
  - Ambiguous (heterozygous) sites stored separately with orientation masks
  
- `conditioning_set H`: Reference panel haplotypes used as "donor" states
  - Stored as bit matrices: `H_opt_var[variant][haplotype]` and `H_opt_hap[haplotype][variant]`
  - Transposed views for cache-efficient access patterns
  
- `variant_map V`: Variant metadata (position, genetic distance, allele frequency)
  
- `hmm_parameters M`: HMM model parameters
  - Transition probabilities `nt[i]`, `yt[i]` pre-computed from genetic map
  - Emission error rates `ed`, `ee` (typically ed/ee ≈ 0.001)

**Supersite Data Structures** (when `--enable-supersites` is used):
- `std::vector<SuperSite> super_sites`: Multi-allelic site metadata
  - `global_site_id`: Anchor variant index (where HMM DP runs)
  - `panel_offset`: Byte offset into packed 4-bit allele codes
  - `var_start`, `var_count`: Member variant indices (CSR-style)
  - `n_alts`: Number of alternate alleles (1-15)
  - `class_prob_offset`: Offset into SC buffer for multivariant posteriors (Phase 3)
  - `n_classes`: Number of allele classes (C = 1 + n_alts)
- `std::vector<int> locus_to_super_idx`: Maps variant index → supersite index (or -1)
- `std::vector<int> super_site_var_index`: Flat array of member variant indices
- `std::vector<uint8_t> packed_allele_codes`: Compact 4-bit codes (2 per byte) for panel haplotypes

**Graph Representation:**
Each sample's genotype is stored as a directed acyclic graph (DAG) of segments:
- Each **segment** represents a region where the phasing is locally determined
- Each segment has a **diplotype** = set of compatible (haplotype0, haplotype1) pairs
- Diplotypes encoded as 64-bit masks where bit d=1 means diplotype d is possible
- Initial state: most segments have all 64 diplotypes possible (fully ambiguous)

### [2] MCMC Iteration Scheme (`phaser_algorithm.cpp::phase()`)

SHAPEIT5 uses an iterative MCMC approach with three stages:

1. **Burn-in** (default: 7 iterations): Sample new phase configurations to explore the state space
2. **Pruning** (default: 8 iterations): Sample + merge segments with high confidence (>95% probability) to reduce complexity
3. **Main** (default: 10 iterations): Sample + store posterior probabilities for final output

Each iteration follows this pattern:
```cpp
for (each iteration) {
    H.select();              // [3] PBWT: Select K conditioning states
    phaseWindow();          // [4-6] Window segmentation + HMM + sampling
    H.updateHaplotypes(G);  // [7] Update reference panel
    H.transposeHaplotypes_H2V();  // Prepare for next PBWT
}
```

### [3] PBWT State Selection (`conditioning_set_selection.cpp`)

**Purpose**: For computational efficiency, the HMM doesn't use all N haplotypes in the panel as states - instead, it selects K≈200 "best" donor haplotypes per sample using the Positional Burrows-Wheeler Transform (PBWT).

**How PBWT Works:**
The PBWT is a data structure that efficiently finds haplotypes that are **identical-by-descent (IBD)** in a local region. It works by:

1. **Prefix Array**: At each variant, maintain a sorted array A where haplotypes are ordered by their suffix (all variants from position i onward)
2. **Divergence Array**: Track where each haplotype last differed from its neighbors
3. **Selection**: For each sample's haplotype, find the K nearest neighbors in the PBWT array (= most similar locally)

**Key Insight**: Haplotypes that share long IBD segments are more likely to be good "copy sources" in the HMM, because they represent recent common ancestry.

**SHAPEIT5 Implementation Details:**
- Multi-threaded PBWT sweep (`conditioning_set::select()`)
- Divides genome into chunks for parallel processing
- At selected "evaluation sites", stores K nearest neighbors per haplotype
- IBD2 tracking: Detects identical-by-descent pairs (parent-child, twins) and bans them from being selected together

**Data Flow:**
```
H_opt_var (variant-major)
    ↓
PBWT forward sweep (partition by allele, track divergence)
    ↓
Select K neighbors at evaluation sites → indexes_pbwt_neighbour[]
    ↓
Transpose to haplotype-major layout for HMM access
```

### [4] Window Segmentation (`window_set.cpp`)

**Purpose**: The genome is too large to phase in one HMM pass, so it's divided into overlapping windows.

**Window Requirements** (for successful split):
- Span ≥4 segments in the genotype graph
- Contain ≥100 variants
- Span ≥genetic distance threshold (default: 0.0001 cM, effectively disabled)

**Recursive Split Algorithm**:
```cpp
bool split(min_length, left_index, right_index) {
    if (too_small) return false;
    
    split_point = random(middle_half_of_range);  // Randomize for MCMC diversity
    
    if (split(left_half) && split(right_half))
        return true;  // Successful recursive split
    else
        return [left_index, right_index];  // Keep as single window
}
```

**Window Metadata** (per window):
- `start_locus`, `stop_locus`: Variant range
- `start_segment`, `stop_segment`: Segment range in genotype graph
- `start_ambiguous`, `stop_ambiguous`: Heterozygous variant range
- `start_missing`, `stop_missing`: Missing variant range
- `start_transition`, `stop_transition`: Transition range in graph

### [5] Multi-threaded HMM Computation (`haplotype_segment_single.cpp`)

This is the computational core of SHAPEIT5. Each sample is processed independently in parallel.

**Per-Thread Data** (`compute_job` class):
- `std::vector<double> T`: Transition probabilities per window
- `std::vector<float> M`: Missing probabilities per window
- **Supersite additions** (Phase 3):
  - `std::vector<float> SC`: CurrentSuperClassPosteriors buffer storing P(class_c | hap) per missing supersite
  - `std::vector<bool> anchor_has_missing`: Per-supersite flag indicating all members missing for this sample
  - Read-only references to global supersite metadata (`super_sites`, `locus_to_super_idx`, `super_site_var_index`)

**Thread Management**:
```cpp
void phaseWindow() {
    for (each sample in parallel) {
        for (each window) {
            // Try single-precision first (faster)
            haplotype_segment_single HS(...);
            HS.forward();
            outcome = HS.backward(...);
            
            // If underflow, retry with double-precision
            if (outcome != 0) {
                haplotype_segment_double HS(...);
                HS.forward();
                outcome = HS.backward(...);
                mark_sample_as_needing_double_precision();
            }
        }
    }
}
```

#### 5a. Forward Pass (`haplotype_segment_single::forward()`)

**Purpose**: Compute α[k] = P(observations up to variant i, state=k) for all variants and all K states.

**Key Data Structures**:
- `prob[k*8 + h]`: Current forward probability for donor k, lane h (8 lanes for AVX2 SIMD)
- `Alpha[s]`: Stored forward probabilities at segment boundaries
- `AlphaMissing[m]`: Stored forward probabilities at missing sites (for imputation)
- `probSumH[h]`: Sum of prob across all donors for lane h
- `probSumT`: Sum across all lanes and donors

**Three HMM Operations** (inlined for each variant type):

1. **INIT**: Initialize first variant in window or segment
   ```cpp
   void INIT_HOM() {  // Homozygous reference or alternate
       for (each donor k) {
           emission = (donor_matches_sample) ? 1.0 : (ed/ee);
           prob[k] = emission;
       }
       probSumT = sum(prob);
   }
   
   void INIT_AMB() {  // Heterozygous (ambiguous orientation)
       for (each donor k) {
           // Use amb_code mask to handle both 0|1 and 1|0 orientations
           // across 8 SIMD lanes simultaneously
           emission_g0 = (donor==0) ? 1.0 : (ed/ee);
           emission_g1 = (donor==1) ? 1.0 : (ed/ee);
           prob[k] = blend(emission_g0, emission_g1, based_on_amb_code);
       }
   }
   
   void INIT_MIS() {  // Missing genotype
       for (each donor k) {
           prob[k] = 1.0;  // All donors equally likely
       }
   }
   ```

2. **RUN**: Propagate within a segment
   ```cpp
   bool RUN_HOM(rare_allele) {
       for (each donor k) {
           emission = (donor_matches_sample) ? 1.0 : (ed/ee);
           
           // Li & Stephens transition:
           // Stay with same donor (nt) or switch to any donor (yt)
           transition = nt * prob[k] + yt * probSumH[lane];
           
           prob[k] = emission * transition;
       }
       probSumT = sum(prob);
       
       // Optimization: skip transition if rare allele doesn't match
       return (sample_has_rare_allele);
   }
   ```

3. **COLLAPSE**: Handle segment boundaries (diplotype changes)
   ```cpp
   void COLLAPSE_HOM() {
       // At segment boundary, previous diplotype constrains which
       // donor haplotypes are compatible
       for (each donor k) {
           emission = (donor_matches_sample) ? 1.0 : (ed/ee);
           
           // probSumK[k] = sum of prob[k] across all diplotypes
           // that include this donor
           transition = nt * probSumK[k] + yt * (1.0 / K);
           
           prob[k] = emission * transition;
       }
   }
   ```

**AVX2 SIMD Optimization**:
All operations process 8 lanes simultaneously using Intel AVX2 intrinsics:
- `__m256` = 8 floats (32 bytes, aligned)
- `_mm256_load_ps`, `_mm256_store_ps`: Load/store aligned float vectors
- `_mm256_mul_ps`, `_mm256_add_ps`: Vectorized arithmetic
- `_mm256_fmadd_ps`: Fused multiply-add (prob * nt + tFreq)

**Lane Semantics (HAP_NUMBER=8)**:
The 8 lanes explore different phase configurations for the diploid sample:
- **Diplotype encoding**: 64 possible diplotypes encoded as (hap0, hap1) pairs where hap0, hap1 ∈ {0..7}
  - `DIP_HAP0(d) = d >> 3` extracts hap0 index (0-7)
  - `DIP_HAP1(d) = d & 7` extracts hap1 index (0-7)
- **Lane assignment**: Each lane h ∈ {0..7} tracks probabilities for specific haplotype indices
- **Biallelic heterozygous sites**: The `amb_code` (8-bit mask from `Ambiguous` array) determines per-lane orientation
  - If `HAP_GET(amb_code, h) == 0`: lane h "wants" haplotype 0's allele (REF or ALT)
  - If `HAP_GET(amb_code, h) == 1`: lane h "wants" haplotype 1's allele (REF or ALT)
  - `genotype::build()` computes `amb_code` based on heterozygous site patterns: `amb_code[h] = (h >> n_unf) % 2` for each HET site
  - With multiple HET sites, different lanes explore different phasing combinations
- **Multiallelic heterozygous sites (supersites)**: Each lane expects a specific **allele class** (not just REF/ALT)
  - `expected_class[h]` = which allele code lane h expects (0=REF, 1=ALT1, 2=ALT2, ...)
  - Computed from `amb_code`: `expected_class[h] = (amb_mask & (1<<h)) ? c1 : c0`
  - Where `c0`, `c1` are the two allele codes carried by the sample's haplotypes
  - Example: genotype `ALT2|ALT3` → c0=2, c1=3 → some lanes expect class 2, others expect class 3
- **Homozygous sites**: All lanes use identical emissions (no orientation ambiguity)
- **Missing sites**: Imputation probabilities computed per-lane, allowing different predictions for different phase configurations

**Example (2 biallelic heterozygous sites)**:
```
Target: [0|1, ./., 0|1] against conditioning haps [0,0,0] and [1,1,1]
Lane semantics after build():
  Lane 0: Explores hap0→cond_hap0, hap1→cond_hap0 (both from REF haplotype)
  Lane 3: Explores hap0→cond_hap1, hap1→cond_hap1 (both from ALT haplotype)
  Lane 1,2,5,6: Explore recombination states (mix of REF/ALT)
Missing site imputation:
  Lane 0: P(ALT) ≈ 0.001 (strongly prefers REF, consistent with cond_hap0)
  Lane 3: P(ALT) ≈ 0.999 (strongly prefers ALT, consistent with cond_hap1)
  Others: P(ALT) ≈ 0.5 (uncertain due to recombination)
```

**Example (multiallelic heterozygous supersite)**:
```
Target: ALT2|ALT3 at multiallelic site (codes: c0=2, c1=3)
Lane semantics:
  Lanes with amb_code bit=0: expect class 2 (ALT2)
  Lanes with amb_code bit=1: expect class 3 (ALT3)
Donor emissions:
  If donor carries class 2: high emission for lanes expecting class 2, low for lanes expecting class 3
  If donor carries class 3: high emission for lanes expecting class 3, low for lanes expecting class 2
  If donor carries REF or other ALT: low emission for all lanes
Result: Same diplotype sampling mechanism, but with multi-allelic class expectations per lane
```

**Storage**:
- At each **segment boundary**: Store `Alpha[s] = prob`, `AlphaSum[s] = probSumH`
- At each **missing variant**: Store `AlphaMissing[m] = prob`, `AlphaSumMissing[m] = probSumH`

#### 5b. Backward Pass (`haplotype_segment_single::backward()`)

**Purpose**: 
1. Compute β[k] = P(observations after variant i | state=k)
2. Compute posterior probabilities: α[k] × β[k]
3. Impute missing genotypes
4. Compute transition probabilities between segments

**Same Operations** (INIT/RUN/COLLAPSE) but in reverse order:
- Start at last variant, move toward first
- Use same emission probabilities
- Use backward transition probabilities (same as forward, but indexed differently)

**Missing Genotype Imputation**:
```cpp
void IMPUTE(missing_probabilities) {
    // Posterior = α[k] × β[k] / sum(α[k] × β[k])
    for (each donor k) {
        for (each lane h) {
            posterior = Alpha[k,h] * prob[k,h] / AlphaSum[h];
            missing_probabilities[h] += posterior;  // Accumulate across donors
        }
    }
    // Result: missing_probabilities[h] = P(allele=1 | data, lane h)
}
```

**Supersite Multivariant Imputation** (Phase 3, when `--enable-supersites` used):
For missing supersites, compute multivariant distribution over all allele classes:
```cpp
void IMPUTE_SUPERSITE_MULTIVARIANT(SC, ss, ss_idx) {
    // For each donor k, accumulate by their allele class
    for (k in donors) {
        class_c = ss_cond_codes[k];  // Which allele this donor carries
        sum[class_c] += Alpha[k] * Beta[k] / AlphaSum;
    }
    
    // Normalize: SC[offset + hap*C + c] = P(class_c | hap)
    for (c in 0..n_classes) {
        for (h in 0..7) {
            SC[offset + h*C + c] = sum[c][h] / total;
        }
    }
}
```
- **Mathematical Guarantee**: Sum over all classes equals 1.0 for each haplotype, ensuring mutual exclusivity
- **Runtime lookups**: 
  - Panel donor code: `unpackSuperSiteCode(panel_codes, ss.panel_offset, hap_idx)`
  - Sample code: `getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap)`
  - Multivariant posterior: `SC[class_prob_offset + hap*C + c]` = P(class_c | hap)
- **Anchor gating**: Run DP only at `SuperSite.global_site_id`; sibling split records skip emission/transition

**Transition Probability Computation**:
At segment boundaries, compute probability of each diplotype transition:
```cpp
P(diplotype_t | data) = 
    (AlphaSum[hap0_prev] * probSumH[hap0_curr] / AlphaSumSum) *
    (AlphaSum[hap1_prev] * probSumH[hap1_curr] / AlphaSumSum)
```

This gives a probability distribution over all valid diplotype transitions.

### [6] Sampling (`genotype_managment.cpp`)

**Purpose**: Use HMM posteriors to select the most likely phase configuration.

**Forward Sampling**:
```cpp
void sampleForward(transition_probs, missing_probs) {
    current_diplotype = sample_from_distribution(transition_probs[0]);
    
    for (each segment) {
        // Sample next diplotype based on transition probabilities
        next_diplotype = sample_from_distribution(
            transition_probs[segment]
        );
        
        current_diplotype = next_diplotype;
    }
    
    // Now make() to apply phasing to variants
    make(sampled_diplotypes, missing_probs);
}
```

**make() Function**:
Applies the sampled diplotype sequence to the variant array:
```cpp
void make(DipSampled, MissingProbs) {
    for (each segment s) {
        hap0 = DIP_HAP0(DipSampled[s]);  // Extract haplotype indices
        hap1 = DIP_HAP1(DipSampled[s]);
        
        for (each variant in segment) {
            if (heterozygous) {
                // Apply orientation from amb_code
                set_phase_from_Ambiguous[amb_idx];
            }
            if (missing) {
                // Sample from posterior probability
                allele0 = (random() < MissingProbs[hap0]) ? 1 : 0;
                allele1 = (random() < MissingProbs[hap1]) ? 1 : 0;
            }
        }
    }
}
```

**Pruning** (optional, in PRUN stage):
After sampling, identify segment transitions with very high confidence (>95%) and merge them, reducing the graph complexity for future iterations.

### [7] Haplotype Update (`conditioning_set_selection.cpp`)

After sampling new phase configurations, update the reference panel:

```cpp
void updateHaplotypes(genotype_set G) {
    for (each sample) {
        for (each variant) {
            H_opt_hap[sample_hap0][variant] = sample.hap0_allele;
            H_opt_hap[sample_hap1][variant] = sample.hap1_allele;
        }
    }
}
```

Then transpose from haplotype-major to variant-major layout for next PBWT iteration:
```cpp
void transposeHaplotypes_H2V() {
    // Block-wise transpose for cache efficiency
    H_opt_var = transpose(H_opt_hap);
}
```

### Performance Optimizations

1. **Precision Fallback**: Start with single-precision (faster), retry with double-precision on underflow
2. **AVX2 SIMD**: Process 8 lanes simultaneously for all HMM operations
3. **Memory Alignment**: All arrays 32-byte aligned for SIMD (`aligned_vector32`)
4. **Cache Efficiency**: Bit-packed storage, transposed matrix views, block-wise operations
5. **Multi-threading**: Independent samples processed in parallel
6. **PBWT**: Reduce states from N haplotypes to K≈200 most relevant
7. **Segment Pruning**: Merge high-confidence segments to reduce graph complexity

### Error Handling

- **Underflow Detection**: If `probSumT` becomes too small, HMM returns error code
- **Precision Escalation**: Automatically retry with double-precision
- **IBD2 Constraints**: Detect and prevent identical samples from being used as mutual donors
- **Validation**: Check for biologically impossible configurations (e.g., trio Mendelian conflicts)

### Debugging & Logging

- **Debug builds**: `make -C phase_common debug` adds `-g` and reduces optimization for easier stepping; run `gdb` or `valgrind` on binaries/tests
- **Logging**: Use `vrb.bullet()` and `vrb.title()` for structured console output. Inspect `Alpha`, `AlphaSum` (32‑byte aligned) in debugger for HMM state
- **Test tracing**: Set `SHAPEIT5_TEST_TRACE=1` to write verbose per‑locus forward/backward TSV traces to `tests/out/`
- **Underflow diagnostics**: Set `SHAPEIT5_DEBUG_UNDERFLOW=1` to append forward/transition underflow events to `logs/underflow.tsv` with sample, locus, cm, yt/nt, probSumT, probSumH[], and prior segment Alpha summaries

### Common Implementation Pitfalls

- **Don't** allocate inside `RUN_AMB` (or any hot path); hoist buffers to class members to avoid allocator churn and `posix_memalign` crashes
- **Don't** perform unaligned AVX2 loads; ensure 32‑byte alignment for `_mm256_load_ps` sources
- **Do** use `aligned_vector32<T>` for all HMM arrays consumed by AVX2 intrinsics

---

## phase_rare Implementation Outline

(This section is intentionally left brief - phase_rare is a separate algorithm not currently the focus of supersite development)

- Uses phased common variants as a scaffold
- Applies simplified HMM only to rare variants
- Lower computational cost than phase_common

---

## Supersite Extension Implementation Outline

### Motivation: Invalid Phasing and Imputation

**Problem 1: Invalid Phasing**
- Haplotypes at multiallelic sites can only have one ALT per haplotype
- Multiallelic variants are common, especially short tandem repeats
- Current workflow splits multiallelic into biallelic (e.g., ALT 2/3 → split1: 0/0, split2: 0/1, split3: 0/1)
- This allows biologically invalid phasing: (0|0, 1|0, 1|0) implies haplotype carries both ALT2 and ALT3

**Problem 2: Invalid Imputation**
- Missing multiallelic calls imputed independently per split
- Allows biologically invalid calls: (./. → 0/1, 1/1, 0/1) implies three alleles at diploid site

**Solution: Atomic Multiallelic Loci**

Supersites treat all split records at the same genomic position (`chr:bp`) as a **single HMM locus**, ensuring the entire multiallelic site is phased and imputed as one atomic unit with exactly one ALT per haplotype.

### High-Level Architecture

When `--enable-supersites` is enabled, the algorithm extends the standard phase_common pipeline:

```
Input VCF/BCF (split multiallelic)
     ↓
[1] Data Loading & Initialization
     ├─→ Standard structures (G, H, V, M)
     └─→ Supersite structures (super_sites, locus_to_super_idx, packed_allele_codes)
     ↓
[2] MCMC Iteration Loop (burn-in → pruning → main)
     ├─→ [3] PBWT State Selection
     ├─→ [3.5] buildSuperSites() ← NEW: Group split records, encode panel alleles
     ├─→ [4] Window Segmentation
     ├─→ [5] Multi-threaded HMM Computation
     │        ├─→ Forward Pass (anchor-gated supersite emissions)
     │        └─→ Backward Pass (multivariant imputation for missing supersites)
     ├─→ [6] Sampling (multivariant sampling + projection to splits)
     └─→ [7] Haplotype Update
     ↓
Output Phased VCF/BCF (mutual exclusivity guaranteed)
```

### [1] Supersite Data Structures

**Extensions to phase_common structures** (when `--enable-supersites` enabled):

- `std::vector<SuperSite> super_sites`: Multi-allelic site metadata
  ```cpp
  struct SuperSite {
      uint32_t global_site_id;      // Anchor variant (where HMM DP runs)
      uint32_t panel_offset;         // Byte offset into packed_allele_codes
      uint32_t var_start, var_count; // Member variant indices (CSR-style)
      uint8_t n_alts;                // Number of alternate alleles (1-15)
      uint32_t class_prob_offset;    // Offset into SC buffer (Phase 3)
      uint8_t n_classes;             // C = 1 + n_alts (cached)
  };
  ```

- `std::vector<int> locus_to_super_idx`: Maps variant index → supersite index (or -1)
- `std::vector<int> super_site_var_index`: Flat array of member variant indices
- `std::vector<uint8_t> packed_allele_codes`: 4-bit allele codes (0=REF, 1-15=ALT1-15), 2 codes per byte

**Per-Thread Extensions** (`compute_job`):
- `std::vector<float> SC`: Multivariant posteriors P(class_c | haplotype) for missing supersites
- `std::vector<bool> anchor_has_missing`: Flags supersites with all members missing

### [3.5] Supersite Preprocessing (`super_site_builder.cpp::buildSuperSites()`)

Called once per MCMC iteration after PBWT state selection.

**Key Operations**:
- **Grouping**: Identify split multiallelic records at same `(chr, bp)` position
- **4-bit Encoding**: Assign codes 0=REF, 1-15=ALT1-15 per conditioning haplotype (first matching ALT wins)
- **Packing**: Store 2 codes per byte at `panel_offset`
- **Chunking**: Sites with >15 alleles split into multiple SuperSite records
- **Lookup Tables**: Build `locus_to_super_idx` and flatten `super_site_var_index`

### [5a] Forward Pass: Anchor-Gated Supersite Emissions

**Supersite Detection & Anchor Gating**:
```cpp
int ss_idx = (locus_to_super_idx) ? (*locus_to_super_idx)[curr_abs_locus] : -1;
if (ss_idx >= 0) {
    const SuperSite& ss = (*super_sites)[ss_idx];
    if (curr_abs_locus != ss.global_site_id) {
        return;  // Sibling split: skip DP, maintain bookkeeping only
    }
    // Anchor: run DP once for entire supersite
}
```

**Sample Classification** (based on both haplotypes, not per-split):
```cpp
uint8_t c0 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap0);
uint8_t c1 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap1);

if (c0 == MISSING && c1 == MISSING) → MIS
else if (c0 == c1)                  → HOM (even if both are ALT5)
else                                → AMB
```

**Emission Computation** (uses 4-bit codes instead of binary alleles):
- **HOM**: `emission = (donor_code == c0) ? 1.0 : (ed/ee)`
- **AMB**: Per-lane expected allele from `amb_code`, compare to `donor_code`
- **MIS**: `emission = 1.0` (all donors equally likely)

**Implementation**: `SS_INIT/RUN/COLLAPSE_{HOM,AMB,MIS}()` functions cache donor codes and apply same Li & Stephens transition/emission math as biallelic.

### [5b] Backward Pass: Multivariant Imputation

**Standard Operations**: Same INIT/RUN/COLLAPSE as forward pass (see [5a]).

**Multivariant Imputation** (Phase 3 extension for missing supersites):

Computes P(class_c | Alpha, Beta) for all allele classes {REF, ALT1, ALT2, ...} in single distribution:

```cpp
void IMPUTE_SUPERSITE_MULTIVARIANT(SC, ss, ss_idx) {
    // Accumulate per class based on donor allele codes
    for (k in 0..n_cond_haps) {
        class_c = ss_cond_codes[k];  // Which allele this donor carries
        sum[class_c] += Alpha[k] * Beta[k] / AlphaSum;
    }
    
    // Normalize: SC stores P(class_c | haplotype)
    for (c in 0..n_classes) {
        for (h in 0..7) {
            SC[offset + h*C + c] = sum[c][h] / total[h];
        }
    }
}
```

**Guarantee**: ∑_c SC[h, c] = 1.0 for each haplotype → exactly one class sampled → mutual exclusivity.

### [6] Sampling: Multivariant Sampling & Projection

**Standard Sampling**: Non-missing supersites use standard biallelic sampling (see phase_common [6]).

**Multivariant Sampling** (for missing supersites):

Sample one allele class per haplotype from multivariant distribution, then project to split records:

```cpp
// Sample class_0 and class_1 from SC multivariant posteriors
float r0 = rng.getDouble();
for (c in 0..n_classes) {
    cumsum += SC[offset + hap0*C + c];
    if (r0 <= cumsum) { class0 = c; break; }
}

// Project: set exactly ONE split to ALT per haplotype
for (ai in 0..n_alts) {
    split_vabs = super_site_var_index[ss.var_start + ai];
    alt_class = ai + 1;  // ALT1=1, ALT2=2, etc.
    
    if (class0 == alt_class) VAR_SET_HAP0(split_vabs);
    else VAR_CLR_HAP0(split_vabs);
}
```

**Result**: Each haplotype has exactly 0 or 1 ALT across all member splits.

### Supersite-Specific Optimizations & Error Handling

**Performance**:
- 4-bit encoding: 50% memory reduction vs byte-per-code
- Anchor gating: DP runs once per biological site (not per split)
- Donor code caching: 10-20% speedup in SS-heavy windows

**Error Handling**: Same underflow detection and precision fallback as biallelic (see phase_common section).

**Design**: Minimal divergence from biallelic - reuse HMM arrays, transition math, SIMD patterns. Zero overhead when `--enable-supersites` is off.

---

## Supersite Coding Patterns (Mirror Biallelic)
- Function taxonomy: Extract `SS_INIT_*`, `SS_RUN_*`, `SS_COLLAPSE_*` with signatures/returns matching biallelic functions; dispatch from existing entry points after centralized classification.
- Centralized classification:
  - Compute `c0`/`c1` via `getSampleSuperSiteAlleleCode` for both haplotypes.
  - Route to MIS if both missing, HOM if `c0==c1`, else AMB.
- AMB lane semantics:
  - Build per-lane want_c0/want_c1 masks from `amb_code`; blend emissions per lane accordingly.
  - Do not split semantics by low/high 128-bit halves; half-splitting is allowed only for packing.
- Anchor gating: If `curr_abs_locus != ss.global_site_id`, treat as sibling and skip DP while maintaining bookkeeping (e.g., missing counters if applicable).
- Naming consistency: Match biallelic locals (`_prob`, `_emit`, `_nt`, `_tFreq`, `_mismatch`, `_factor`) and accumulators (`probSumH`, `probSumT`, `probSumK`).

---

## Debugging & Logging
- Debug builds: `make -C phase_common debug` adds `-g` and reduces optimization for easier stepping; run `gdb` or `valgrind` on binaries/tests.
- Logging: use `vrb.bullet()` and `vrb.title()` for structured console output. Inspect `Alpha`, `AlphaSum` (32‑byte aligned) in debugger for HMM state.
- Tests: set `SHAPEIT5_TEST_TRACE=1` to write verbose per‑locus forward/backward TSV traces to `tests/out/`.
- HMM underflow tracing: set `SHAPEIT5_DEBUG_UNDERFLOW=1` to append forward/transition underflow events to `logs/underflow.tsv` with sample, locus, cm, yt/nt, probSumT, probSumH[], and prior segment Alpha summaries.

## Common Pitfalls (Do/Don't)
- Don’t allocate inside `RUN_AMB` (or any hot path); hoist buffers to class members to avoid allocator churn and `posix_memalign` crashes.
- Don’t access supersite state without null/size checks (`super_sites`, `locus_to_super_idx`).
- Don’t perform unaligned AVX2 loads; ensure 32‑byte alignment for `_mm256_load_ps` sources.
- Do gate DP to anchors to prevent double counting at sibling splits.
- Do preserve AMB lane semantics from biallelic code; avoid half-lane logic divergence.

Alignment items for supersites:
- Don't overload `rare_allele` to signal supersite siblings (e.g., `rare_allele=2`). Use explicit supersite gating (`locus_to_super_idx`, `super_sites`) consistently in HMM paths.
- Don't start windows at supersite siblings when `--enable-supersites` is set; adjust window starts to the anchor or a non‑member locus.

---

## 8-Lane System: Correct Semantics (CRITICAL)

### What the 8 Lanes Actually Represent

**FUNDAMENTAL CONCEPT**: The 8 lanes are **NOT** 8 independent haplotypes or direct indices into diplotypes. Instead, they represent:

**8 different phase configuration hypotheses** for how heterozygous sites within a segment might phase together.

### Mathematical Foundation

#### Diplotype Encoding (64 possible states)
Each segment has a `Diplotypes[s]` bit mask encoding which of 64 possible diplotype configurations are valid:
- Each diplotype d ∈ {0..63} encodes a pair: `(hap0, hap1)` where hap0, hap1 ∈ {0..7}
- `DIP_HAP0(d) = d >> 3` extracts hap0 index (which of 8 lanes represents haplotype 0)
- `DIP_HAP1(d) = d & 7` extracts hap1 index (which of 8 lanes represents haplotype 1)
- **Key insight**: A diplotype maps **two of the 8 lanes** to the two haplotypes

#### Lane-to-Phase Pattern Mapping (Heterozygous Sites)
Within a segment containing n_het heterozygous sites (where n_het ≤ 3):
- Each lane h ∈ {0..7} represents a specific **binary pattern** of phase configurations
- For HET site i within the segment: lane h's phase is determined by bit (h >> i) & 1
- This creates 2^n_het unique phase configurations across the n_het sites

**Example with 2 heterozygous sites**:
```
Segment contains: [HET_site_1, HET_site_2]

Lane 0 (binary 000): HET1 phase bit=0, HET2 phase bit=0 → both sites prefer haplotype 0's allele
Lane 1 (binary 001): HET1 phase bit=1, HET2 phase bit=0 → HET1 prefers hap1, HET2 prefers hap0
Lane 2 (binary 010): HET1 phase bit=0, HET2 phase bit=1 → HET1 prefers hap0, HET2 prefers hap1
Lane 3 (binary 011): HET1 phase bit=1, HET2 phase bit=1 → both sites prefer haplotype 1's allele
Lanes 4-7: unused (only 2^2 = 4 configurations needed)
```

#### How Ambiguous Array Encodes This
In `genotype_build.cpp`, the `Ambiguous[amb_idx]` byte is constructed as:
```cpp
for (unsigned int h = 0; h < HAP_NUMBER; h++) {
    bool allele = ((h>>n_unf)%2);  // Extract bit (n_unf) from lane index h
    if (allele) HAP_SET(Ambiguous[amb_idx], h);
}
```
This means:
- For HET site at position n_unf=0: lanes {1,3,5,7} get bit set (pattern: 01010101)
- For HET site at position n_unf=1: lanes {2,3,6,7} get bit set (pattern: 00110011)
- For HET site at position n_unf=2: lanes {4,5,6,7} get bit set (pattern: 00001111)

**Semantic Interpretation**:
- `HAP_GET(amb_code, h) == 0` → lane h's hypothesis is that this HET site phases with haplotype 0
- `HAP_GET(amb_code, h) == 1` → lane h's hypothesis is that this HET site phases with haplotype 1

#### How Biallelic INIT_AMB/RUN_AMB Works Correctly
```cpp
unsigned char amb_code = G->Ambiguous[curr_abs_ambiguous];
for (int h = 0; h < HAP_NUMBER; h++) {
    g0[h] = HAP_GET(amb_code,h) ? M.ed/M.ee : 1.0f;  // If lane expects hap1's allele: mismatch for REF
    g1[h] = HAP_GET(amb_code,h) ? 1.0f : M.ed/M.ee;  // If lane expects hap1's allele: match for ALT
}
// Then for each donor k with allele ah:
__m256 _prob = _emit[ah];  // Use g0 if donor has REF (ah=0), g1 if donor has ALT (ah=1)
```

**Why this works for biallelic**:
- Biallelic sites have exactly 2 alleles: REF (0) and ALT (1)
- Donor allele ah ∈ {0, 1} directly indexes which emission vector to use
- Each lane h gets emission based on whether it expects hap0's allele or hap1's allele
- For a 0|1 genotype: some lanes explore 0|1 phase, others explore 1|0 phase

---

## CRITICAL BUGS IN SUPERSITE 8-LANE IMPLEMENTATION (Nov 2, 2025)

### **BUG #7: Fundamental Misinterpretation of Lane Semantics for Multiallelic Sites**

**Status**: ❌ **CRITICAL - BLOCKS ALL HETEROZYGOUS MULTIALLELIC PHASING**

**Location**: All `SS_*_AMB` functions in `haplotype_segment_single.h` and `haplotype_segment_double.h`

**The Error**:
The supersite code interprets `amb_code` as a direct lane-to-allele-class mapping:
```cpp
// INCORRECT INTERPRETATION:
uint8_t expected_class[HAP_NUMBER];
unsigned char amb_mask = G->Ambiguous[curr_abs_ambiguous];
for (int h = 0; h < HAP_NUMBER; ++h) {
    bool use_c1 = ((amb_mask >> h) & 1U);
    expected_class[h] = use_c1 ? c1 : c0;
}
```

**What this code assumes**:
- Each lane h directly represents a specific haplotype (0 or 1)
- `amb_mask bit h` tells us which haplotype lane h represents
- We can assign the allele class (c0 or c1) that lane h's haplotype carries

**The actual semantics**:
- Each lane h represents a **phase configuration hypothesis** across all HET sites in the segment
- `amb_mask bit h` tells us whether this HET site's phase configuration in lane h aligns with hap0 or hap1
- **We cannot assign a single "expected class" to a lane** because that lane represents a hypothesis about **multiple sites' phasing**

**Why Biallelic Works But Multiallelic Fails**:

For **biallelic** 0|1 genotype at a single HET site:
- Lanes with amb_code bit=0 explore hypothesis "this site is phased as hap0=REF, hap1=ALT"
- Lanes with amb_code bit=1 explore hypothesis "this site is phased as hap0=ALT, hap1=REF"
- Donor with allele=0 (REF) matches lanes expecting hap0's allele (bit=0), mismatches lanes expecting hap1's allele (bit=1)
- The `_emit[ah]` array correctly provides per-lane emissions because there are only 2 alleles

For **multiallelic** ALT2|ALT3 genotype (c0=2, c1=3):
- Lanes with amb_code bit=0 explore hypothesis "this site is phased as hap0=ALT2, hap1=ALT3"
- Lanes with amb_code bit=1 explore hypothesis "this site is phased as hap0=ALT3, hap1=ALT2"
- **But current code assigns**: lanes with bit=0 get expected_class=2, lanes with bit=1 get expected_class=3
- **This breaks when**:
  - Donor carries ALT1 (code=1): Should mismatch all lanes, but code might assign match if expected_class happens to be 1
  - Donor carries REF (code=0): Should mismatch all lanes (neither haplotype is REF), but code treats it as distinct from both expected classes
  - **The fundamental issue**: We're comparing 4-bit codes (0-15) when we should be considering **which allele the donor has relative to which allele the lane's hypothesis expects for that haplotype**

**Impact**:
- **All heterozygous multiallelic sites get incorrect emission probabilities**
- Lanes that should favor certain donor haplotypes instead favor wrong donors
- Phase configurations that match the true biological phase get lower probability than incorrect configurations
- **Result**: Extremely high error rates (likely >50%) for all HET multiallelic supersites

**Root Cause**:
The biallelic emission model cannot be directly extended to multiallelic sites because:
1. Biallelic: 2 alleles → can use binary indexing `_emit[ah]` where ah ∈ {0,1}
2. Multiallelic: C alleles → cannot use binary indexing, need per-lane class expectations
3. **The lane semantics are about phase patterns, not about individual haplotypes**

---

### **BUG #8: Multiallelic Sites Cannot Use Binary Phase Pattern Encoding**

**Status**: ❌ **ARCHITECTURAL - REQUIRES REDESIGN**

**Location**: `genotype_build.cpp`, Ambiguous array construction

**The Problem**:
The entire `Ambiguous` encoding system assumes **binary alleles**:
```cpp
for (unsigned int h = 0; h < HAP_NUMBER; h++) {
    bool allele = ((h>>n_unf)%2);  // BINARY: either 0 or 1
    if (allele) HAP_SET(Ambiguous[amb_idx], h);
}
```

**Why this is insufficient for multiallelic**:
- Encodes only 2 possible allele states per site (0 or 1)
- For ALT2|ALT3 genotype, need to encode which lanes expect class 2 vs class 3
- Current encoding loses information about **which specific ALT alleles** are involved
- Cannot represent phase patterns like "lane 0 expects ALT2, lane 1 expects ALT3, lane 2 expects ALT2, ..."

**Example Failure**:
```
Genotype at multiallelic site: ALT2|ALT3 (c0=2, c1=3)

Current Ambiguous encoding (binary):
  Lanes 0,2,4,6: amb_code bit=0 → assigned expected_class=2
  Lanes 1,3,5,7: amb_code bit=1 → assigned expected_class=3

What's missing:
  - No encoding of which diplotype configurations are valid
  - No way to represent that hap0 carries ALT2 and hap1 carries ALT3
  - Cannot distinguish between:
    * hap0=ALT2, hap1=ALT3 (actual genotype)
    * hap0=ALT3, hap1=ALT2 (swapped phase)
    * hap0=ALT2, hap1=ALT2 (homozygous, should be in HOM path)
```

**Architectural Limitation**:
The segment/diplotype/lane system was designed for **at most 3 binary (heterozygous) sites per segment** (2^3 = 8 lanes). For multiallelic sites:
- Need to track which of C allele classes each haplotype carries
- Cannot represent this in 8-bit amb_code
- Would need C^2 possible diploid configurations, not 2^n_het lane patterns

---

### **BUG #9: Diplotype Masks Don't Constrain Multiallelic Configurations**

**Status**: ❌ **ARCHITECTURAL - MISSING VALIDATION**

**Location**: `genotype_build.cpp`, Diplotypes array construction

**The Issue**:
The `Diplotypes[s]` bitmask is built based on **binary heterozygous site patterns**:
```cpp
for (unsigned int vrel = 0; vrel < Lengths[s]; vrel++) {
    bool f_het = VAR_GET_HET(MOD2(vabs+vrel),Variants[DIV2(vabs+vrel)]);
    if (f_het) {
        switch (n_unf) {
        case 0: Diplotypes[s] &= MASK_UNF0; break;  // Constrains to 32 of 64 diplotypes
        case 1: Diplotypes[s] &= MASK_UNF1; break;  // Constrains to 16 of 64 diplotypes
        case 2: Diplotypes[s] &= MASK_UNF2; break;  // Constrains to 8 of 64 diplotypes
        }
    }
    n_unf += f_het;
}
```

**What's happening**:
- MASK_UNF0 = 0x55AA55AA55AA55AAUL: constrains diplotypes based on 1st HET site
- MASK_UNF1 = 0x3333CCCC3333CCCCUL: further constrains based on 2nd HET site
- MASK_UNF2 = 0x0F0F0F0FF0F0F0F0UL: further constrains based on 3rd HET site
- These masks ensure valid lane-to-diplotype mappings for binary alleles

**What's broken for multiallelic**:
- Multiallelic HET sites are marked as VAR_GET_HET (indistinguishable from biallelic HET)
- The same binary masks get applied
- **But**: A diplotype (hap0=2, hap1=3) cannot be represented in this system!
  - `hap0=2` means "lane 2" but lanes represent phase patterns, not allele classes
  - The diplotype masks encode **which phase configurations are valid**, not **which allele classes are valid**
- **Result**: Invalid diplotype configurations get marked as valid, leading to incorrect sampling

**Example**:
```
Segment with 1 biallelic HET (0|1) and 1 multiallelic HET (ALT2|ALT3):

After applying masks:
  - MASK_UNF0 applied for first HET: 32 diplotypes remain
  - MASK_UNF1 applied for second HET: 16 diplotypes remain
  
But those 16 diplotypes encode phase patterns like:
  - Diplotype 0 (0,0): both sites phase with hap0
  - Diplotype 9 (1,1): both sites phase with hap1
  - etc.

None of these diplotypes encode "hap0 carries ALT2, hap1 carries ALT3"!
The diplotype system is orthogonal to the allele class system.
```

---

### **BUG #10: Segment Boundaries Don't Align With Multiallelic Constraints**

**Status**: ❌ **ARCHITECTURAL - PHASE CORRELATION LOSS**

**Location**: `genotype_build.cpp`, segment boundary logic

**The Problem**:
Segments are created based on:
1. Reaching 4 unfolded sites (HET or SCAFFOLD)
2. Reaching max variant count (65535)
3. Reaching max ambiguous count (255)

**For multiallelic sites**:
- Each multiallelic anchor counts as 1 HET site (if heterozygous)
- Segments can span multiple multiallelic sites
- **But**: The lane system cannot track phase correlations across different multiallelic sites

**Example Failure**:
```
Segment contains:
  - Variant 100: biallelic 0|1
  - Variant 101: multiallelic ALT2|ALT3 (supersite anchor)
  - Variant 102: multiallelic ALT1|ALT4 (different supersite anchor)

Lane 0's hypothesis:
  - For variant 100: hap0=0, hap1=1
  - For variant 101: ??? (cannot encode "hap0=ALT2, hap1=ALT3" in lane pattern)
  - For variant 102: ??? (cannot encode "hap0=ALT1, hap1=ALT4" in lane pattern)

The amb_code for variant 101 encodes a binary pattern (bit 0 = expect hap0's allele)
But "hap0's allele" is ALT2, which is not represented anywhere in the lane system!
```

**Architectural Mismatch**:
- Lanes represent **binary phase patterns** across HET sites
- Multiallelic sites require **C-way allele class assignments** per haplotype
- These two systems are fundamentally incompatible within the same segment

---

### **BUG #11: COLLAPSE Functions Don't Have Multiallelic Diplotype Context**

**Status**: ❌ **CRITICAL - INCORRECT SEGMENT TRANSITIONS**

**Location**: `SS_COLLAPSE_AMB` in `haplotype_segment_single.h`

**The Issue**:
At segment boundaries, COLLAPSE functions compute:
```cpp
__m256 _prob = _mm256_set1_ps(probSumK[k]);  // Sum across all diplotypes that include donor k
```

**What `probSumK[k]` represents** (from biallelic code):
```cpp
// In forward() when storing Alpha at segment boundary:
for (int k = 0; k < n_cond_haps; k++) {
    probSumK[k] = 0.0f;
    for (int h = 0; h < HAP_NUMBER; h++) {
        // Sum prob[k*8+h] for all lanes h that are valid in this diplotype configuration
        probSumK[k] += prob[k*8 + h];  // (Simplified - actual code checks diplotype masks)
    }
}
```

**Why this breaks for multiallelic**:
- `probSumK[k]` sums across lanes representing different **phase configuration hypotheses**
- For multiallelic, different lanes have different `expected_class[h]` values
- **Summing across lanes loses the allele class information**
- When transitioning to next segment, the collapsed probability has no memory of which allele classes were involved

**Example**:
```
Previous segment ends at multiallelic site with genotype ALT2|ALT3:
  - Lane 0-3 had expected_class = {2,3,2,3,...}
  - Lane 4-7 had expected_class = {3,2,3,2,...}
  - probSumK[k] = sum of all 8 lanes for donor k
  - This sum mixes lanes expecting different allele classes!

Next segment begins:
  - COLLAPSE function broadcasts probSumK[k] to all 8 lanes
  - All lanes get same initial probability regardless of allele class
  - **Loss of phase information from previous multiallelic site**
```

---

### **BUG #12: Missing Data Imputation Cannot Reconstruct Multiallelic Genotypes**

**Status**: ⚠️ **FUNCTIONAL - PHASE 3 ADDRESSED IMPUTATION BUT NOT PHASING**

**Location**: `genotype_managment.cpp`, `make()` function

**Partial Fix Status**:
Phase 3 (Oct 29, 2025) implemented multivariant imputation via `IMPUTE_SUPERSITE_MULTIVARIANT` which correctly:
- Computes P(class_c | haplotype) for all classes c
- Samples exactly one class per haplotype
- Projects to split records ensuring mutual exclusivity

**Remaining Issue**:
The **phasing** of non-missing heterozygous multiallelic sites still relies on the broken lane system:
```cpp
// In make() for non-missing heterozygous variants:
unsigned char hap0 = DIP_HAP0(DipSampled[s]);  // Which lane represents haplotype 0
unsigned char hap1 = DIP_HAP1(DipSampled[s]);  // Which lane represents haplotype 1

// For heterozygous sites, phase is determined by Ambiguous array
// But for multiallelic sites, this doesn't encode which allele each haplotype carries!
```

**The Gap**:
- Imputation (missing data) works correctly via multivariant posteriors
- Phasing (non-missing HET data) uses broken lane→allele class mapping
- **Result**: Non-missing HET multiallelic sites can still get incorrectly phased

---

## Summary of Architectural Issues

### Fundamental Incompatibility

The current architecture has an **irreconcilable conflict**:

1. **Segment/Lane System** (designed for binary alleles):
   - 8 lanes encode 2^3 = 8 possible binary phase patterns
   - Lanes represent hypotheses about how up to 3 HET sites phase together
   - Works perfectly for biallelic: allele ∈ {REF, ALT} maps to phase ∈ {hap0, hap1}

2. **Multiallelic System** (requires C-way allele encoding):
   - Need to track which of C allele classes each haplotype carries
   - C^2 possible diploid configurations (not 2^n_het phase patterns)
   - Cannot represent "hap0=ALT2, hap1=ALT3" in lane-based phase pattern system

### Why Supersite Phase 1-3 Couldn't Fix This

- **Phase 1-2**: Focused on atomic HMM treatment and data structures
  - Successfully treats multiallelic as single HMM locus ✓
  - Successfully prevents double-counting via anchor gating ✓
  - **But**: Inherited broken lane semantics from biallelic code ✗

- **Phase 3**: Implemented multivariant imputation
  - Successfully imputes missing multiallelic genotypes ✓
  - Guarantees mutual exclusivity via multivariant sampling ✓
  - **But**: Non-missing HET multiallelic still uses broken phasing ✗

### Impact on Error Rates

**Expected error rates by genotype class**:
- **HOM multiallelic** (e.g., ALT2|ALT2): Should work correctly ✓
  - Uses SS_*_HOM functions which compare donor code to single sample code
  - No lane ambiguity, straightforward emission computation
  
- **HET multiallelic, fully missing** (./. at all splits): Should work correctly ✓
  - Uses multivariant imputation (Phase 3 implementation)
  - Samples from correct posterior distribution
  
- **HET multiallelic, partially observed** (e.g., 0/1 at one split): **BROKEN** ✗
  - Classification sees different alleles in c0 vs c1
  - Routes to SS_*_AMB functions
  - Applies incorrect lane→expected_class mapping
  - **Extremely high error rates expected (>50%)**

- **HET multiallelic, fully observed** (e.g., 0/1 at split1, 0/0 at split2): **BROKEN** ✗
  - Same lane semantics issue as partially observed
  - **Extremely high error rates expected (>50%)**

---

## Supersite Implementation Status & Roadmap

**Done:**
- ✅ Data structures (SuperSite.class_prob_offset, compute_job.SC, anchor_has_missing)
- ✅ Backward pass multivariant computation (IMPUTE_SUPERSITE_MULTIVARIANT)
- ✅ Genotype projection (sample from multivariant, project to splits)
- ✅ Integration (phaser_algorithm passes SC, calls setSuperSiteContext)
- ✅ Clean compilation (both single and double precision)
- ✅ AMB lane semantics (per-lane expected class using G.Ambiguous mask)

**Newly Discovered Bugs (Unnecessary Divergence from Biallelic):**

**BUG #1: Duplicate SS_*_MIS functions**
- **Issue**: `SS_INIT_MIS()`, `SS_RUN_MIS()`, and `SS_COLLAPSE_MIS()` are identical to biallelic `INIT_MIS()`, `RUN_MIS()`, `COLLAPSE_MIS()`
- **Impact**: Code duplication, maintenance burden
- **Fix**: Have supersite MIS classification directly call biallelic MIS functions (missing data is representation-agnostic)
- **Status**: ✅ FIXED (Nov 2, 2025)
- **Implementation**: All 6 dispatcher switch statements in `haplotype_segment_{single,double}.h` (INIT/RUN/COLLAPSE for HOM/AMB) now call biallelic `INIT_MIS()`, `RUN_MIS()`, `COLLAPSE_MIS()` directly instead of `SS_*_MIS()` wrappers
- **Verification**: Clean compilation of phase_common and tests

**BUG #2: Inconsistent sibling handling between INIT and RUN/COLLAPSE**
- **Issue**: `INIT_HOM/AMB` at siblings call `INIT_MIS()`, but `RUN/COLLAPSE` just return early without state update
- **Impact**: Asymmetric behavior - INIT siblings get neutral probabilities, RUN/COLLAPSE siblings keep stale state
- **Fix**: Standardize - all functions now call appropriate MIS functions at siblings (uniform approach)
- **Status**: ✅ FIXED (Nov 2, 2025)
- **Implementation**: All RUN_HOM/AMB and COLLAPSE_HOM/AMB sibling gates now call `RUN_MIS()` and `COLLAPSE_MIS()` respectively instead of early return
- **Rationale**: Siblings at same chromosomal position have yt≈0, so MIS functions act as identity operations (copy forward anchor probabilities). This ensures correct bookkeeping and genetic distance tracking while having negligible computational cost.
- **Verification**: Clean compilation of phase_common and tests

**BUG #3: Missing bookkeeping updates at sibling loci**
- **Issue**: Sibling loci return early without updating `curr_abs_ambiguous`, `curr_abs_missing`, `AlphaMissing`, etc.
- **Impact**: Missing data imputation and segment indexing may be incorrect for supersites
- **Fix**: Add bookkeeping updates even when skipping DP (as documented in AGENTS.md patterns)
- **Status**: ⏳ TODO

**BUG #4: Unclear transition probability normalization in COLLAPSE**
- **Issue**: Biallelic `COLLAPSE_HOM` has comment "//Check divide by probSumT here!" suggesting uncertainty
- **Impact**: Unclear if formula is correct; hard to verify
- **Status**: ⚠️ **EXPERIMENTAL FLAG ADDED** - needs empirical validation
- **Investigation**: Mathematical analysis suggests current formula (no normalization) is likely correct because:
  - Switching represents uniform draw over K donors (independent of prev probability mass)
  - Total probability conserved: Σ_k[probSumK[k]*(nt/probSumT) + yt/K] = nt + yt = 1.0
  - However, asymmetry with RUN functions is visually suspicious
- **Experimental flag**: `SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION=1` enables alternative normalization
  - Default (env not set): `_tFreq = yt / n_cond_haps` (current behavior)
  - With flag: `_tFreq = (yt * probSumT) / n_cond_haps` (normalized like RUN)
- **TODO**: 
  - Run accuracy benchmarks with/without flag to determine optimal behavior
  - **Remove flag and hardcode optimal behavior before release**
  - Replace comment with clear mathematical explanation

**BUG #5: Duplicate emission logic in SS_COLLAPSE_AMB**
- **Issue**: Two separate 30+ line code paths for `ss_anchor_split_emissions` mode vs full-supersite mode
- **Impact**: Code duplication, maintenance burden
- **Fix**: Extract common transition code, parameterize only emission computation
- **Status**: ✅ FIXED (Nov 2, 2025) - Unified loop with precomputed emission vectors in both single and double precision
- **Implementation**: Refactored `SS_COLLAPSE_AMB` in both `haplotype_segment_single.h` and `haplotype_segment_double.h`
  - Single loop processes all donors with precomputed expected class array
  - Transition code (`probSumK[k] * nt + tFreq`) unified
  - Emission computation parameterized by mode (split vs strict 4-bit equality)
  - Reduced ~60 lines duplicate code to ~45 lines unified path per file
- **Verification**: Clean compilation with -O3 -mavx2 -mfma optimization flags

**BUG #6: Emission application divergence**
- **Issue**: Biallelic uses inline conditional multiplication, supersite precomputes emission vectors
- **Impact**: Different code paths for same operation; harder to verify correctness
- **Analysis**: The divergence exists because biallelic operates on binary (0/1) alleles while supersites operate on 4-bit codes (0-15) with per-lane semantics
- **Mathematical equivalence**:
  - Biallelic: `if (ag!=ah) _prob *= mismatch` ≡ `emit = (ag==ah) ? 1.0 : (ed/ee)`
  - Supersite: `emit = _mm256_blendv_ps(mis_f, match_f, match_mask)` ≡ `emit[h] = (donor_code==expected_class[h]) ? 1.0 : (ed/ee)`
- **Design decision**: **Keep both patterns** (Option C)
  - Biallelic: Optimized for simple binary comparisons, minimal branching overhead
  - Supersite: Required for per-lane class semantics, uses SIMD blend operations
  - Both patterns implement identical Li & Stephens emission probabilities
  - Rationale documented in code comments for reviewer clarity
- **Status**: ✅ DOCUMENTED (Nov 2, 2025) - No code change needed; patterns are mathematically equivalent and each optimized for their use case

**Known Implementation Debt:**
- Window starts may land on supersite siblings; ideally adjust to anchors or non-member loci when `--enable-supersites` is set
- Consider adding assertion to prevent window starts on siblings in debug builds

---

## TODO: Scaffolding Support for Supersites

### Overview
SHAPEIT5 supports trio/duo scaffolding (`--pedigree` option) where known parent-child relationships are used to phase heterozygous sites. This feature currently assumes biallelic variants and will require significant modification to support multiallelic supersites.

### Current Scaffolding Implementation (`genotype_mendel.cpp`)

**Trio Scaffolding Logic:**
```cpp
if (VAR_GET_HET(child)) {
    if (VAR_GET_HOM(father) && VAR_GET_HOM(mother)) {
        // Both parents homozygous → determine phase from parent alleles
        bool fath0 = VAR_GET_HAP0(father);
        bool moth0 = VAR_GET_HAP0(mother);
        if (fath0 != moth0) {
            VAR_SET_SCA(child);  // Set to scaffold (phased)
            // Copy parent alleles to child haplotypes
        }
    } else if (VAR_GET_HOM(father)) {
        // Only father homozygous → father's allele goes to one haplotype
    }
}
```

**Key Assumption:** Biallelic encoding where HET means 0/1 or 1/0, and HAP0/HAP1 bits directly encode the allele values.

### Pain Points for Multiallelic Support

#### 1. **Allele Code Mapping Across Family Members**

**Problem:** Parent and child may have different split representations of the same multiallelic site.

**Example Scenario:**
```
Position 1000: REF=A, ALT1=T, ALT2=G, ALT3=C

Father: A/T (codes: 0,1)
  - Split 1 (A/T): 0/1 → HET
  - Split 2 (A/G): 0/0 → HOM
  - Split 3 (A/C): 0/0 → HOM

Mother: G/G (codes: 2,2)
  - Split 1 (A/T): 0/0 → HOM
  - Split 2 (A/G): 1/1 → HOM
  - Split 3 (A/C): 0/0 → HOM

Child: T/G (codes: 1,2) - inherits T from father, G from mother
  - Split 1 (A/T): 0/1 → HET
  - Split 2 (A/G): 0/1 → HET
  - Split 3 (A/C): 0/0 → HOM
```

**Challenge:** Current code compares `VAR_GET_HOM(parent)` at each split independently, but we need to:
- Identify that this is the same multiallelic site across family members
- Translate allele codes between different individuals' split representations
- Determine inheritance pattern at the multiallelic level, not per-split

#### 2. **Mendelian Violation Detection**

**Problem:** Need to detect impossible inheritance patterns at multiallelic sites.

**Valid vs Invalid Inheritance:**

✅ **Valid (biallelic):**
- Father: 0/0, Mother: 0/0 → Child: 0/0 ✓
- Father: 0/1, Mother: 0/0 → Child: 0/0 or 0/1 ✓
- Father: 0/0, Mother: 0/0 → Child: 0/1 ✗ (violation)

✅ **Valid (multiallelic):**
- Father: ALT1/ALT1 (1,1), Mother: ALT2/ALT2 (2,2) → Child: ALT1/ALT2 (1,2) ✓
- Father: REF/ALT1 (0,1), Mother: ALT2/ALT3 (2,3) → Child: ALT1/ALT2 (1,2) ✓

✗ **Invalid (multiallelic):**
- Father: ALT1/ALT1 (1,1), Mother: ALT2/ALT2 (2,2) → Child: REF/ALT3 (0,3) ✗
- Father: REF/REF (0,0), Mother: REF/REF (0,0) → Child: ALT1/ALT2 (1,2) ✗

**Challenge:** Current biallelic logic uses simple HAP bit comparisons. Multiallelic requires:
- Checking if child's allele codes are subset of union(father_codes, mother_codes)
- Handling de novo mutations (child has allele not in either parent)
- Counting violations at the multiallelic level, not per-split

#### 3. **Scaffold Bit Setting Across Splits**

**Problem:** When scaffolding succeeds, need to mark ALL member splits as scaffolded, not just anchor.

**Current Behavior (biallelic):**
```cpp
VAR_SET_SCA(MOD2(v), Variants[DIV2(v)]);  // Set single variant
VAR_SET_HAP0(...);  // Set haplotype 0 allele
VAR_CLR_HAP1(...);  // Clear haplotype 1 allele
```

**Required Behavior (multiallelic):**
- Determine which allele code child inherits from each parent (e.g., child gets code 1 from father, code 2 from mother)
- Set SCA bit at **anchor** position
- Set HAP0/HAP1 bits at **appropriate member splits** based on allele codes
- Example: If child is 1|2, set HAP0 at split 1, set HAP1 at split 2, clear all others

**Challenge:** Need atomic scaffolding operation across all supersite members.

#### 4. **Partial Scaffolding Scenarios**

**Problem:** One parent may be heterozygous at different splits than the child.

**Example:**
```
Father: ALT1/ALT2 (1,2)
  - Split 1: 0/1 → HET
  - Split 2: 0/1 → HET

Child: ALT1/ALT3 (1,3)
  - Split 1: 0/1 → HET
  - Split 3: 0/1 → HET
```

**Current biallelic approach:** "If one parent is HOM, scaffold from that parent"

**Multiallelic challenge:**
- Father has ALT1 and ALT2; child has ALT1 and ALT3
- Child's ALT1 must come from father (successful constraint)
- Child's ALT3 cannot come from father → must come from mother
- Need partial scaffolding: constrain one haplotype, leave other ambiguous

#### 5. **genotype::build() Integration**

**Problem:** If we update anchor HET bits to reflect multiallelic heterozygosity, scaffolded sites need special handling.

**Proposed Anchor HET Update:**
```cpp
// After buildSuperSites(), update anchor HET bits based on multiallelic genotype
if (multiallelic_genotype_is_HET) {
    VAR_SET_HET(anchor);
}
```

**Scaffolding Interaction:**
- Scaffolding runs AFTER VCF reading but BEFORE genotype::build()
- If scaffold sets anchor to SCA, should we still set HET bit?
- Or should SCA take precedence (scaffolded sites are "known" phase)?

**Decision needed:** Should scaffold sites be excluded from HET bit update?

#### 6. **getSampleSuperSiteAlleleCode() Correctness After Scaffolding**

**Current Behavior:**
```cpp
uint8_t getSampleSuperSiteAlleleCode(const genotype* G, const SuperSite& ss, 
                                     const std::vector<int>& super_site_var_index, int hap) {
    for (uint32_t i = 0; i < ss.var_count; ++i) {
        int v_idx = super_site_var_index[ss.var_start + i];
        bool carries = (hap == 0) ? VAR_GET_HAP0(v_idx, v) : VAR_GET_HAP1(v_idx, v);
        if (carries) return static_cast<uint8_t>(i + 1);  // Return ALT index
    }
    return SUPERSITE_CODE_REF;
}
```

**After Scaffolding:**
- Scaffolding sets HAP0/HAP1 bits at specific splits
- `getSampleSuperSiteAlleleCode()` scans all splits to find which carries ALT
- **This should still work correctly** ✓

**But verify:** If child is scaffolded as 1|2:
- Split 1 has HAP0=1, HAP1=0 → hap0 scan returns code 1 ✓
- Split 2 has HAP0=0, HAP1=1 → hap1 scan returns code 2 ✓

#### 7. **Haploid Sample Handling**

**Current Code (`genotype_reader_reading.cpp`):**
```cpp
if (options.count("haploids")) {
    G.resetHaploidHeterozgotes(readerH.samples);
}
```

**Challenge:** Haploid samples can't be heterozygous at multiallelic sites either.
- Need to ensure haploid samples don't get marked as HET at supersite anchors
- Or, ensure scaffolding logic skips haploid samples at supersites

### Implementation Strategy

#### Phase 1: Detection and Validation
1. **Detect supersite trio inconsistencies**
   - Add supersite check in `scaffoldTrio()`: `if (locus_to_super_idx && (*locus_to_super_idx)[v] >= 0)`
   - For supersites, retrieve allele codes for all three family members
   - Validate Mendelian consistency at multiallelic level
   - Count violations separately (biallelic vs multiallelic)

2. **Short-term fix: Skip scaffolding at supersites**
   ```cpp
   for (int v = 0; v < n_variants; v++) {
       // Skip all supersite members
       if (locus_to_super_idx && (*locus_to_super_idx)[v] >= 0) continue;
       
       // Existing biallelic scaffolding logic
   }
   ```

#### Phase 2: Multiallelic-Aware Scaffolding
1. **Create `scaffoldSuperSiteTrio()` function**
   - Input: SuperSite, child genotype, father genotype, mother genotype
   - Retrieve allele codes for all three individuals
   - Compute valid inheritance combinations
   - Return: (success, child_hap0_code, child_hap1_code)

2. **Implement allele code translation**
   - Map parent allele codes to child's split representation
   - Handle cases where parent lacks certain ALT alleles

3. **Atomic scaffold operation**
   - Set SCA bit at anchor
   - Set HAP0/HAP1 bits at appropriate member splits
   - Clear HAP bits at non-inherited splits

#### Phase 3: Integration
1. **Update genotype::build() to handle scaffolded supersites**
   - Don't update HET bit for scaffolded sites (SCA takes precedence)?
   - Or, update HET bit but ensure build() logic handles SCA+HET correctly?

2. **Extend haplotype_set::updateHaplotypes()**
   - Ensure scaffolded supersite haplotypes are correctly written to panel

3. **Add comprehensive tests**
   - Trio with multiallelic sites (valid inheritance)
   - Mendelian violations at multiallelic sites
   - Mixed biallelic/multiallelic trios

### Open Questions
1. **How should scaffold file support multiallelic representations?**
   - Currently scaffold file provides biallelic phased genotypes
   - Would need to specify multiallelic allele codes in scaffold format

### Files Requiring Modification

1. **`genotype_mendel.cpp`**: Core scaffolding logic
   - `scaffoldTrio()`, `scaffoldDuoFather()`, `scaffoldDuoMother()`
   - Add supersite detection and handling

2. **`genotype_reader_reading.cpp`**: Scaffold file reading
   - May need to handle multiallelic scaffold representations

3. **`super_site_accessor.h`**: Add helper functions
   - `validateSuperSiteMendelianInheritance()`
   - `scaffoldSuperSite()`

4. **`phaser_initialise.cpp`**: Integration
   - Ensure scaffolding runs before `updateSuperSiteAnchorHetBits()` (if implemented)

5. **Tests**: Add trio scaffolding test cases
   - `tests/data/trio_multiallelic.vcf`
   - Validation scripts

### Priority

- **Short-term (P0)**: Add skip logic to prevent incorrect scaffolding at supersites
- **Medium-term (P1)**: Implement multiallelic-aware trio validation and Mendelian violation detection
- **Long-term (P2)**: Full multiallelic scaffolding support with allele code translation

**Status**: 📋 **TODO** - Not yet implemented; currently biallelic-only scaffolding would produce incorrect results at supersites

---

## Supersite Runtime Configuration

**CLI Flags:**
- `--enable-supersites`: Enable multiallelic phasing support (default: off)
- `--ss-anchor-split-emissions`: At supersite anchors, use biallelic split-site emission semantics (treat other-ALT donors as REF at the anchor split). Default: off (strict 4-bit class equality per lane)

**Environment Variables:**
- `SHAPEIT5_TEST_TRACE=1`: Enable verbose per-locus forward/backward TSV traces to `tests/out/`
- `SHAPEIT5_DEBUG_UNDERFLOW=1`: Enable HMM underflow diagnostics to `logs/underflow.tsv` (includes sample, locus, cm, yt/nt, probSumT, probSumH[], and prior segment Alpha summaries)
- **`SHAPEIT5_NORMALIZE_COLLAPSE_TRANSITION=1`**: **⚠️ EXPERIMENTAL - REMOVE BEFORE RELEASE** 
  - Enables alternative COLLAPSE transition normalization for Bug #4 testing
  - Default (not set): `_tFreq = yt / K` (current behavior, likely correct)
  - With flag: `_tFreq = (yt * probSumT) / K` (symmetric with RUN functions)
  - Purpose: Empirically determine which formula produces more accurate phasing
  - See Bug #4 in "Known Bugs" section for mathematical analysis

**Key Implementation Files:**
- `phase_common/src/models/haplotype_segment_single.{h,cpp}`: Float HMM core + supersite integration
- `phase_common/src/models/haplotype_segment_double.{h,cpp}`: Double HMM core + supersite integration
- `phase_common/src/models/super_site_macros.h`: Emission macro templates (INIT/RUN/COLLAPSE)
- `phase_common/src/objects/super_site_builder.{h,cpp}`: Multi-allelic grouping + panel encoding
- `phase_common/src/phaser/phaser_algorithm.cpp`: Window segmentation, threading, precision fallback
- `common/src/utils/otools.h`: `aligned_vector32`, toolbox externs
- `tests/makefile`: Unit test build rules, external object linking


## Side note re: compilation 
- This is a many cpu machine, and I recommend making with -j24 argument
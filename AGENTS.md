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

The 8 lanes represent potential haplotype orientations (for heterozygous sites) or are duplicated (for homozygous sites).

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
- Don’t overload `rare_allele` to signal supersite siblings (e.g., `rare_allele=2`). Use explicit supersite gating (`locus_to_super_idx`, `super_sites`) consistently in HMM paths.
- Don’t start windows at supersite siblings when `--enable-supersites` is set; adjust window starts to the anchor or a non‑member locus.

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
- **Fix**: Document whether COLLAPSE should normalize by `probSumT` or not (check Li & Stephens math)
- **Status**: ⏳ TODO - needs mathematical verification

**BUG #5: Duplicate emission logic in SS_COLLAPSE_AMB**
- **Issue**: Two separate 30+ line code paths for `ss_anchor_split_emissions` mode vs full-supersite mode
- **Impact**: Code duplication, maintenance burden
- **Fix**: Extract common transition code, parameterize only emission computation
- **Status**: ⏳ TODO

**BUG #6: Emission application divergence**
- **Issue**: Biallelic uses inline conditional multiplication, supersite precomputes emission vectors
- **Impact**: Different code paths for same operation; harder to verify correctness
- **Fix**: Unify approach - either precompute everywhere or use inline everywhere
- **Status**: ⏳ TODO - needs design decision

**Known Implementation Debt:**
- Window starts may land on supersite siblings; ideally adjust to anchors or non-member loci when `--enable-supersites` is set
- Consider adding assertion to prevent window starts on siblings in debug builds

---

## Supersite Runtime Configuration

**CLI Flags:**
- `--enable-supersites`: Enable multiallelic phasing support (default: off)
- `--ss-anchor-split-emissions`: At supersite anchors, use biallelic split-site emission semantics (treat other-ALT donors as REF at the anchor split). Default: off (strict 4-bit class equality per lane)

**Environment Variables:**
- `SHAPEIT5_TEST_TRACE=1`: Enable verbose per-locus forward/backward TSV traces to `tests/out/`
- `SHAPEIT5_DEBUG_UNDERFLOW=1`: Enable HMM underflow diagnostics to `logs/underflow.tsv` (includes sample, locus, cm, yt/nt, probSumT, probSumH[], and prior segment Alpha summaries)

**Key Implementation Files:**
- `phase_common/src/models/haplotype_segment_single.{h,cpp}`: Float HMM core + supersite integration
- `phase_common/src/models/haplotype_segment_double.{h,cpp}`: Double HMM core + supersite integration
- `phase_common/src/models/super_site_macros.h`: Emission macro templates (INIT/RUN/COLLAPSE)
- `phase_common/src/objects/super_site_builder.{h,cpp}`: Multi-allelic grouping + panel encoding
- `phase_common/src/phaser/phaser_algorithm.cpp`: Window segmentation, threading, precision fallback
- `common/src/utils/otools.h`: `aligned_vector32`, toolbox externs
- `tests/makefile`: Unit test build rules, external object linking

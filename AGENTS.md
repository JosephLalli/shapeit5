## Project
SHAPEIT5 — `phase_common` super-site integration

## Source of Truth
- This AGENTS.md mirrors `.github/copilot-instructions.md`, which is the authoritative, detailed guide distilled from `SUPERSITE_CONVERSATION_SUMMARY.md`.
- If anything diverges, prefer `.github/copilot-instructions.md` and update this file to match.
- Runtime binaries depend on locally installed Boost/HTSlib libraries; before running any tests or executable directly, export `LD_LIBRARY_PATH="$HOME/.linuxbrew/lib:/usr/local/lib:$LD_LIBRARY_PATH"` so the dynamic loader can resolve shared libraries (e.g., `libboost_iostreams.so.1.89.0`).

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

**Implementation**: Cache donor codes and apply same Li & Stephens transition/emission math as biallelic.

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

## Testing Conventions

```bash
make -C tests                      # Build all test binaries
make -C tests test-run             # Build + run all (sets LD_LIBRARY_PATH)
LD_LIBRARY_PATH=$HOME/.linuxbrew/lib tests/bin/test_supersite_hmm_states
LD_LIBRARY_PATH=$HOME/.linuxbrew/lib tests/bin/test_segment_boundary_multiallelic
```

- **Test harnesses** use `#define private public` to expose internals for state validation
- **`test_toolbox.cpp`** defines global objects (`rng`, `vrb`, etc.) via `_DECLARE_TOOLBOX_HERE` (instantiate exactly once)
- Synthetic data: no external fixtures yet (placeholders under `tests/data/`)
- `test_segment_boundary_multiallelic` now repeats the six-variant motif twice so the fixture always creates two segments; the transition-probability vector therefore contains one normalized block per segment (sum of entire vector ≈ `G.n_segments`, each block sums to 1.0). When debugging, print block boundaries rather than assuming the whole vector sums to 1.

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

## Supersite Implementation Status & Roadmap

**Done:**
- ✅ Data structures (SuperSite.class_prob_offset, compute_job.SC, anchor_has_missing)
- ✅ Backward pass multivariant computation (IMPUTE_SUPERSITE_MULTIVARIANT)
- ✅ Genotype projection (sample from multivariant, project to splits)
- ✅ Integration (phaser_algorithm passes SC, calls setSuperSiteContext)
- ✅ Clean compilation (both single and double precision)
- ✅ AMB lane semantics (per-lane expected class using G.Ambiguous mask)

**Previously addressed bugs/design antipatterns**

**Issue #1: Duplicate SS_*_MIS functions**
- **Issue**: `SS_INIT_MIS()`, `SS_RUN_MIS()`, and `SS_COLLAPSE_MIS()` are identical to biallelic `INIT_MIS()`, `RUN_MIS()`, `COLLAPSE_MIS()`
- **Impact**: Code duplication, maintenance burden
- **Fix**: Have supersite MIS classification directly call biallelic MIS functions (missing data is representation-agnostic)
- **Status**: ✅ FIXED (Nov 2, 2025)
- **Implementation**: All 6 dispatcher switch statements in `haplotype_segment_{single,double}.h` (INIT/RUN/COLLAPSE for HOM/AMB) now call biallelic `INIT_MIS()`, `RUN_MIS()`, `COLLAPSE_MIS()` directly instead of `SS_*_MIS()` wrappers
- **Verification**: Clean compilation of phase_common and tests

**Issue #2: Inconsistent sibling handling between INIT and RUN/COLLAPSE**
- **Issue**: `INIT_HOM/AMB` at siblings call `INIT_MIS()`, but `RUN/COLLAPSE` just return early without state update
- **Impact**: Asymmetric behavior - INIT siblings get neutral probabilities, RUN/COLLAPSE siblings keep stale state
- **Fix**: Standardize - all functions now call appropriate MIS functions at siblings (uniform approach)
- **Status**: ✅ FIXED (Nov 2, 2025)
- **Implementation**: All RUN_HOM/AMB and COLLAPSE_HOM/AMB sibling gates now call `RUN_MIS()` and `COLLAPSE_MIS()` respectively instead of early return
- **Rationale**: Siblings at same chromosomal position have yt≈0, so MIS functions act as identity operations (copy forward anchor probabilities). This ensures correct bookkeeping and genetic distance tracking while having negligible computational cost.
- **Verification**: Clean compilation of phase_common and tests

**Issue #3: Missing bookkeeping updates at sibling loci**
- **Issue**: Sibling loci return early without updating `curr_abs_ambiguous`, `curr_abs_missing`, `AlphaMissing`, etc.
- **Impact**: Missing data imputation and segment indexing may be incorrect for supersites
- **Fix**: Add bookkeeping updates even when skipping DP
- **Status**: Fixed (Nov 5, 2025)

**Issue #4: Transition probability normalization in COLLAPSE**
- **Issue**: Historical uncertainty on whether to scale `_tFreq` by `probSumT` in COLLAPSE
- **Resolution**: Fixed. COLLAPSE uses `_tFreq = yt / n_cond_haps` (uniform over K). No env flag.
- **Rationale**:
  - Switching is a uniform draw over donors; previous mass already captured in the stay term.
  - Conservation: Σ_k[probSumK[k]*(nt/probSumT) + yt/K] = nt + yt = 1.0
  - Matches biallelic and supersite paths; simplifies implementation.

**Issue #5: Duplicate emission logic in SS_COLLAPSE_AMB**
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

**Issue #6: Emission application divergence**
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
- **Verification**: Clean compilation with -O3 -mavx2 -mfma optimization flags

**Known Implementation Debt:**
- Window starts may land on supersite siblings; ideally adjust to anchors or non-member loci when `--enable-supersites` is set
- Consider adding assertion to prevent window starts on siblings in debug builds

---

## TODO (P3): Scaffolding Support for Supersites

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

**Haploid Challenge:** Haploid samples can't be heterozygous at multiallelic sites either.
- Need to ensure haploid samples don't get marked as HET at supersite anchors
- Or, ensure scaffolding logic skips haploid samples at supersites

### Trio/scaffolding Implementation Strategy

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


### Priority

- **Short-term (P0)**: Verify skip logic to prevent incorrect scaffolding at supersites
- **Medium-term (P1)**: Implement multiallelic-aware trio validation and Mendelian violation detection
- **Long-term (P2)**: Full multiallelic scaffolding support with allele code translation

**Status**: 📋 **TODO** - Not yet implemented; currently biallelic-only scaffolding would produce incorrect results at supersites

---

## Supersite Runtime Configuration

Note: some of these may be out of date.

**CLI Flags:**
- `--enable-supersites`: Enable multiallelic phasing support (default: off)
- `--ss-anchor-split-emissions`: At supersite anchors, use biallelic split-site emission semantics (treat other-ALT donors as REF at the anchor split). Default: on (parity mode; strict 4-bit class equality retained for storage/imputation)

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


## Side note re: compilation 
- This is a many cpu machine, and I recommend making with -j24 argument

---

## Supersite ↔ Biallelic Parity

### Problem Summary

- We want the supersite representation (split multiallelic records grouped by position) to be functionally equivalent to the biallelic representation at the level of the HMM’s anchor loci while retaining strict multiallelic information for imputation.
- Discrepancies were observed when comparing a 5‑variant biallelic setup to a 10‑variant supersite setup (5 anchors + 5 siblings):
  - Extra sibling steps previously introduced additional transitions and per‑step normalization. Fixes applied:
    - Siblings are now true no‑ops (no emission/transition math, do not advance `prev_abs_locus`).
    - Recombination at identical map positions is set to zero (`yt=0`) by returning `0` when `cm` ties.
  - Emission semantics mismatch at anchors:
    - Requirement: Keep strict 4‑bit class codes for storage/imputation, but present biallelic semantics to the HMM at anchors: “is the anchor ALT present?”
    - Solution in progress: At anchors with `ss_anchor_split_emissions` enabled, unify supersite and biallelic HMM by using the same mask‑based routines driven by the per‑site `Ambiguous` mask (per‑lane orientation) and the donor “anchor‑ALT” flags.
  - Despite these, normalized forward posteriors at anchors were still off in tests. Root cause under investigation: ensure the exact same mask and HMM code paths are used at anchors in both reps.

### Changes Implemented

- Transitions and siblings
  - `hmm_parameters`: return `yt=0` at identical genetic map positions; clamp only tiny positive distances.
  - HMM (double/single): siblings are no‑ops (no renormalization and no `prev_abs_locus` advancement).

- Anchor emissions (presentation) vs storage (strict class)
  - Tests enable `M.ss_anchor_split_emissions = true` to present biallelic semantics at anchors while strict class codes remain for Phase 3.
  - Double‑precision HMM now unifies anchor INIT/RUN/COLLAPSE under split semantics using the same biallelic mask‑based routines:
    - Build a per‑donor×lane `MatchMask` via `SupersiteEmissionAdapter::build_match_mask(..., /*use_anchor_split_semantics*/true, ...)`.
    - Call `INIT_FROM_MASK`, `RUN_FROM_MASK`, `COLLAPSE_FROM_MASK` at anchors.

- Instrumentation (tests)
  - Added anchor‑only, forward‑only windows and comparison of normalized per‑lane posteriors.
  - Dump lane expectations (from `Ambiguous`), donor flags (5‑var ALT vs SS anchor‑ALT), and mask matrices to diagnose mismatches.

### Rollout Plan

1) Source of truth for lane expectations at anchors
   - Add `amb_mask` to `SiteView` and set it in both biallelic and supersite `build_view()` at anchors.
   - Make `SupersiteEmissionAdapter::build_match_mask()` use `view.amb_mask` for expected per‑lane ALT when `use_anchor_split_semantics=true`.

2) Unify HMM code paths (single precision)
   - Mirror the double‑precision unification in `haplotype_segment_single.{h,cpp}` by adding `RUN_FROM_MASK` and `COLLAPSE_FROM_MASK`, and using them at anchors under split semantics.

3) Strong parity assertions (tests only)
   - At each anchor:
     - Assert mask matrix equality (per donor × lane) between 5‑var and 10‑var.
     - Assert normalized α equality per donor × lane (forward‑only at window end).
     - Run full backward and assert normalized α×β equality per donor × lane.
   - Gate strict assertions by an env flag for CI tuning.

4) Microcase and progressive expansion
   - Start with a tiny 3‑locus microcase (HOM–AMB–HOM) and 4 donors; prove parity, then progress to multi‑anchor scenarios.

5) Guardrails
   - Siblings remain true no‑ops; `yt=0` at identical positions.
   - Strict class codes preserved for Phase 3 imputation and sampling; only anchor emissions are reduced to biallelic semantics when configured.

### Acceptance Criteria

- At anchors, for windows ending at the anchor:
  - Per‑lane emission masks are identical between supersite and biallelic representations.
  - Normalized α and α×β posteriors match within tight tolerances.
  - Segment‑boundary transition distributions match to numerical precision.

- End results (phasing and multivariant imputation) agree across representations (tests already present: mutual exclusivity, SC normalization).

### Flags & Knobs

- `--ss-anchor-split-emissions` (hmm_parameters): enable biallelic presentation at supersite anchors (strict 4‑bit storage retained).
- `SHAPEIT5_TEST_TRACE=1`: emit trace lines (including `build_view` and anchor diagnostics) to help debug parity.

### Bug Retrospective & Fixes

#### Supersite SC Buffer Race Condition

**Issue**: Substantial accuracy reduction when supersites are enabled, manifested as K values (conditioning set size) increasing instead of decreasing during HMM iterations. Investigation revealed a race condition where multiple threads simultaneously modify `SuperSite.class_prob_offset` in `compute_job::make()`, causing memory corruption and heap scribbling that affects neighboring data structures like IBD2 tracks.

**Root Cause**: The `class_prob_offset` field in the `SuperSite` struct was being concurrently written by multiple worker threads without synchronization, leading to:
- Overlapping SC buffer writes between threads
- Memory corruption affecting adjacent data structures  
- Incorrect offset calculations for supersite class posterior storage
- K value divergence due to corrupted conditioning set management

**Fix**: Redesigned supersite SC buffer architecture to use thread-local storage:
- **Removed** `class_prob_offset` field from `SuperSite` struct to eliminate shared write target
- **Added** `supersite_sc_offset` vector to `compute_job` class for per-thread offset management
- **Updated** all `IMPUTE_SUPERSITE_MULTIVARIATE` calls to accept thread-local offset parameter
- **Modified** `setSuperSiteContext()` to use thread-local offset storage
- **Added** race condition verification assertions in debug builds

**Files Modified**:
- `phase_common/src/models/super_site_accessor.h` - Removed race condition target field
- `phase_common/src/objects/compute_job.{h,cpp}` - Added thread-local offset management
- `phase_common/src/models/haplotype_segment_{single,double}.{h,cpp}` - Updated function signatures
- `phase_common/src/objects/genotype/genotype_{header.h,managment.cpp}` - Updated context setting
- Multiple test files - Fixed compilation after struct field removal

**Result**: Eliminates memory corruption while preserving full supersite functionality and maintaining thread safety.

#### Supersite Expansion Parity

Summary of issues encountered while making the 10‑variant supersite setup (5 anchors + 5 siblings) parity‑equivalent to the 5‑biallelic setup at anchors, and the fixes applied:

- Sibling renormalization (clamping in test harness)
  - Issue: test code clamped identical cm to 1e‑7, causing `yt>0` between siblings and “effectively normalizing” forward state.
  - Fix: in tests, treat identical positions as `yt=0, nt=1`; clamp only tiny positive distances. Production already returns `0` for `dist<=0` and only clamps tiny positives.
  - Guardrail: siblings remain true no‑ops; do not advance `prev_abs_locus` (no unintended transitions/renormalization).

- Anchor emission semantics mismatch
  - Issue: supersite anchors used strict 4‑bit class equality vs. biallelic binary REF/ALT, producing differing α.
  - Fix: enable `ss_anchor_split_emissions` and unify to mask‑based INIT/RUN/COLLAPSE at anchors for both precisions. Expected per‑lane ALT comes from `Ambiguous` mask; donor flag is “is anchor ALT present?”. Strict class storage retained for Phase 3.

- Suspected lane inversion at locus 0 (INIT parity)
  - Observation: early on, first‑locus α appeared lane‑swapped.
  - Root cause/Resolution: after unifying mask semantics and checking parity of `Ambiguous` and donor anchor‑ALT flags, parity holds. Ambiguous mask is deterministic from `genotype::build()` via `n_unf`.
  - Instrumentation: test‑only `INIT_FROM_MASK` prints (mask, per‑lane α, probSumH) confirm correct mapping.

- Ambiguous index and mask parity at anchors
  - Issue: mismatched ambiguous indices or amb_mask could silently diverge emissions.
  - Fix: compute amb_idx consistently (skip supersite siblings) and assert amb_idx parity and full mask parity (verbose‑gated in tests).

- Dataset initialization pitfalls (toy data)
  - Issue: if siblings aren’t strictly 0|0 in sample/panel, or sibling cm aren’t identical, parity fails even with correct HMM code.
  - Fix: verbose‑gated assertions in the test enforce: identical cm for sibling pairs; sample siblings 0|0; panel siblings REF across donors.

- Forward path gating
  - Fixes: at window start, siblings INIT neutrally; within windows and at collapses, siblings set `update_prev_locus=false`. This preserves forward/backward probabilities across sibling → anchor with no change.


# Our current task:
# **"BUG #11": COLLAPSE Functions Don't Have Multiallelic Diplotype Context**

## Task: Preserve and use column-side marginals at segment boundaries in `phase_common`

### Goal
At genotype-graph segment boundaries in `phase_common`, stop broadcasting only the row marginal (`probSumK`) when seeding the next segment. Instead, also carry the column marginal (`probSumH`, 8 SIMD lanes) from the previous segment and seed via a row×column outer product with recombination mixing, then apply the (two-class) supersite emission mask `{a,b}` built from 4-bit allele IDs.

### Why
Broadcasting row marginals flattens the second-chain prior (column side) and erases per-lane class skew created by the previous segment—this is harmless-ish for biallelic ref|alt but can be materially wrong for supersite hets like alt1|alt2. Carrying `probSumH` fixes this with negligible overhead and unchanged vectorization.

---

## Codebase
- `main/phase_common/src/models/haplotype_segment_single.{h,cpp}`
- `main/phase_common/src/phaser/phaser_algorithm.cpp`
- (If needed) headers that define per-segment storage structs for forward states

Assume AVX2 `_mm256_*` intrinsics are used today; keep the same structure.

---

## High-level changes

1) Persist per-segment column marginals
   - Add `AlphaLaneSum`: per-segment 8-lane sums (column marginal) at segment end.
   - You already compute `probSumH` (lane sums) and `probSumT` (total). Persist them at the moment you finalize a segment in the forward pass.

2) Outer-product seeding in `COLLAPSE_*`
   - Replace broadcast seeding (`set1(probSumK[k])`) with:
     - `row_mix[k] = nt * (prevProbSumK[k] / prevProbSumT) + yt * (1.0f / K)`
     - `col_mix[8] = nt * (prevProbSumH[j] / prevProbSumT) + yt * (1.0f / LANES)`  (LANES=8)
     - `base(k, lanes) = row_mix[k] * col_mix(lanes)` (vector multiply)
   - Then apply the two-class emission mask for the supersite `{a,b}` using 4-bit donor allele IDs.

3) Backward consistency
   - In `SET_FIRST_TRANS` for the first segment in a window, compute initial diplotype weights from the product of marginals and apply the same `{a,b}` emission mask. (This mirrors forward seeding.)

4) Feature flag / fallback
   - Keep the current behavior (broadcast) as a fallback code path (e.g., if a compile-time flag is off, or if the previous segment didn’t produce valid sums).
   - Default: ON when supersites are active; otherwise keep current path for classic biallelic sites (identical result either way).

---

## Concrete edits

### 1) Header: store lane marginals per segment
File: `haplotype_segment_single.h`

Add members (aligned for AVX2):

// At top of class
static constexpr int LANES = 8;

struct alignas(32) Lane8 { float v[LANES]; };

std::vector<Lane8> AlphaLaneSum;   // per-seg column marginals (probSumH)
std::vector<float> AlphaTot;       // per-seg totals (probSumT) if not already present

Initialize/resize alongside `Alpha` / `AlphaSum` in the constructor or `init()`.

### 2) Forward: persist lane sums at segment end
File: `haplotype_segment_single.cpp`

In the forward pass, wherever you finalize a segment (right after computing `probSumH` and `probSumT` for that segment), persist:

// Assume seg_idx is the index of the segment just finalized
for (int j = 0; j < LANES; ++j) AlphaLaneSum[seg_idx].v[j] = probSumH[j];
AlphaTot[seg_idx] = probSumT;

This should live right where you already store `Alpha` / `AlphaSum` for the segment.

### 3) Build multiallelic (two-class) emission vectors
Still in `haplotype_segment_single.cpp`, add a helper that produces lane masks from donor 4-bit IDs for the next segment’s first locus (the supersite).

inline __m256 make_mask_eq_class(const uint8_t lane_cls[LANES], uint8_t a, float eps) {
    alignas(32) float m[LANES];
    for (int j=0;j<LANES;++j) m[j] = (lane_cls[j]==a) ? 1.0f : eps;
    return _mm256_load_ps(m);
}

You’ll need the donors’ allele class per lane at the supersite. If today you derive a 1-bit `amb_code`, generalize to a `uint8_t lane_cls[LANES]` fetched from the supersite metadata (each 4-bit ID ∈ {0..15}).

### 4) Replace broadcast in `COLLAPSE_HOM/AMB/MIS`
Find the three collapse routines in `haplotype_segment_single.cpp` (or the switch that calls them) where you currently do:

__m256 prob = _mm256_set1_ps(probSumK[k]);            // ← broadcast
// mixing:
prob = nt * prob / prevProbSumT + yt * (1.0f / K);
// emission vector multiply (2-class today)

Replace with outer-product seeding:

const float invK = 1.0f / K;
const float invL = 1.0f / LANES;
const float invT = 1.0f / prevProbSumT;

// Build column mix once per collapse
alignas(32) float col_mix_arr[LANES];
for (int j=0;j<LANES;++j)
    col_mix_arr[j] = nt * (prevProbSumH[j] * invT) + yt * invL;
const __m256 col_mix = _mm256_load_ps(col_mix_arr);

for (int k=0; k<K; ++k) {
    // Row mixing
    const float row_mix_k = nt * (prevProbSumK[k] * invT) + yt * invK;

    // Outer product prior (rank-1)
    __m256 base = _mm256_mul_ps(_mm256_set1_ps(row_mix_k), col_mix);

    // Emission: two-class {a,b} from 4-bit IDs
    // Inputs needed here:
    //   cls_k: 4-bit class ID for donor k at this supersite
    //   lane_cls[LANES]: 4-bit class IDs for each SIMD lane at this supersite
    __m256 emit;
    if (isHOM) {
        if (cls_k == a) emit = make_mask_eq_class(lane_cls, a, eps);
        else            emit = _mm256_set1_ps(eps);      // whole row penalized
    } else if (isHET) {
        if      (cls_k == a) emit = make_mask_eq_class(lane_cls, b, eps);
        else if (cls_k == b) emit = make_mask_eq_class(lane_cls, a, eps);
        else                 emit = _mm256_set1_ps(eps); // whole row penalized
    } else { // MIS
        emit = _mm256_set1_ps(1.0f);
    }

    __m256 out = _mm256_mul_ps(base, emit);
    _mm256_store_ps(&prob[k*LANES], out);
}

// Then recompute probSumH[LANES] and probSumT exactly as today (sum across rows and lanes).
// Keep existing underflow scaling / renormalization.

HOM/AMB/MIS specifics:
- HOM: a == b
- AMB: genuine het {a,b}, a != b
- MIS: emission is all ones; or keep current missing logic

### 5) Backward: `SET_FIRST_TRANS`
Where you set the initial diplotype distribution for the first segment, current code likely multiplies independent per-hap marginals already. Ensure you:
- Use both per-hap marginals (row and column) derived from the (backward) state at that boundary.
- Multiply by the same `{a,b}` emission mask to zero/penalize illegal pairs.

This keeps forward/backward consistent.

### 6) Thread the data
- Ensure `prevProbSumK` and `prevProbSumH` used in `COLLAPSE_*` refer to the previous segment’s stored values:
  - `prevProbSumK` is your existing per-row sums (e.g., `AlphaSum[seg-1]` or equivalent).
  - `prevProbSumH` is `AlphaLaneSum[seg-1].v`.
  - `prevProbSumT` is `AlphaTot[seg-1]`.
- Guard for `seg==0` (use neutral prior—your existing `INIT_*` path for the first segment).

---

## Tests

1) Three-donor toy (K=3, LANES=8 but populate only 3 lanes):
   - Prev segment supersite: target 0/1; donors at end: k1=0, k2=2, k3=1.
   - Next segment supersite: target 1/2; donors at start: k1=1, k2=0, k3=1.
   - Set `e=1e-4`, low recomb `y=1e-3`.
   - Expected: column marginal for class 2 (prev end) ≪ class 0/1; with outer-product, legal (1,2)/(2,1) cells start much smaller than broadcast; after emission+normalize, posterior differs from broadcast by >20% relative.

2) Biallelic sanity: ref|alt → ref|alt with balanced donors.
   - Outer-product and broadcast produce near-identical normalized distributions (relative diff < 1e-3).

3) Performance sanity: Same runtime (±1–2%) vs baseline on a real window (e.g., K≈64). No extra allocations in inner loop.

---

## Notes / constraints
- Keep AVX2 pattern: one `_mm256_set1_ps(row_mix_k)` per row, one vector multiply with `col_mix`, then emission vec multiply. No per-lane branching.
- `lane_cls[LANES]` must be available for the supersite (per donor lane). If not present, fall back to biallelic masks (current code path).
- Memory overhead: `AlphaLaneSum`: 32 bytes × (#segments) per individual. Negligible.
- Backward underflow: unchanged; still use your current scaling and double-precision fallback.

---

## Acceptance criteria
- Compiles and passes existing tests.
- New unit tests confirm: multiallelic supersite case diverges from broadcast in the expected direction; biallelic case remains effectively unchanged.
- No measurable slowdown and no increase in numerical underflow incidents.

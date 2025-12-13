## Project
SHAPEIT5 — `phase_common` supersite (multiallelic) integration

## Source of Truth
- This AGENTS.md mirrors `.github/copilot-instructions.md`, which is the authoritative, detailed guide distilled from `.AGENT_markdowns/SUPERSITE_CONVERSATION_SUMMARY.md`.
- If anything diverges, prefer `.github/copilot-instructions.md` and update this file to match.
- Runtime binaries depend on locally installed Boost/HTSlib libraries; before running any tests or executable directly, export `LD_LIBRARY_PATH="$HOME/.linuxbrew/lib:/usr/local/lib:$LD_LIBRARY_PATH"` so the dynamic loader can resolve shared libraries (e.g., `libboost_iostreams.so.1.89.0`).
- For a deep dive into missing-variant imputation (biallelic and supersite), see `.AGENT_markdowns/MISSING_LOGIC.687eb2d.md`, which documents the implementation as of commit `687eb2d` (later commits may have changed behavior slightly).

## Context
- SHAPEIT5 uses Li–Stephens (phase_common/phase_rare) to phase and impute genotypes.
- Phase_common is the full HMM; phase_rare is a scaffolded HMM for rare variants.
- Emissions use a default mismatch penalty `ed/ee ≈ 0.001`; transitions use `nt/yt` from the genetic map.

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

This documentation focuses on **phase_common** and its supersite support.

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

### [5] Multi-threaded HMM Computation (`haplotype_segment_single.*`)

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

### [7] Haplotype Update

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

### Debugging & Logging (General)

- **Debug builds**: `make -C phase_common debug` adds `-g` and reduces optimization for easier stepping; run `gdb` or `valgrind` on binaries/tests
- **Logging**: Use `vrb.bullet()` and `vrb.title()` for structured console output. Inspect `Alpha`, `AlphaSum` (32‑byte aligned) in debugger for HMM state
- **Test tracing**: `SHAPEIT5_TEST_TRACE=1` enables verbose, per‑locus console traces. Useful prints include:
  - Supersite adapter: `[SupersiteEmit] build_view(super)` with anchor id, c0/c1, and lane_class per lane
  - Emission masks: `build_match_mask locus=…` with donor raw codes and any_match_lane summary
  - Bial AMB path: `BIAL_RUN_AMB: locus=…` with g0/g1 per‑lane emission vectors and per‑donor emits (first few donors)
  - Supersite AMB path: `SS_RUN_AMB: locus=…` with expected classes per lane, donor codes, and strict class‑equality per‑donor emits
  - Cursor diagnostics: `[ss-amb-cursor]` range/advances; sibling bookkeeping via `[SupersiteSibling]`
- **Underflow diagnostics**: `SHAPEIT5_DEBUG_UNDERFLOW=1` appends forward/transition underflow events to `logs/underflow.tsv` with sample, locus, cm, yt/nt, probSumT, per‑lane probSumH, and prior‑segment Alpha summaries
- **Targeted supersite debug**: set both `SHAPEIT5_SUPERDEBUG_SAMPLENAME=<sample>` and `SHAPEIT5_SUPERDEBUG_BP=<anchor_locus>` to print compact, focused snapshots around a specific anchor (per‑split hap bits before/after projection, sampled h0/h1, etc.)
- **Packing trace**: `SHAPEIT5_TRACE_SUPERSITE_PACKING=1` prints panel packing access (unpack indices) for early calls to validate bounds and offsets
- **Supersite guard checks**: `SHAPEIT5_SUPERSITE_GUARDS=1` (default) enables bounds/consistency checks in supersite accessors; set `0` to disable

### Notes on legacy toggles

- `--ss-anchor-split-emissions` (parity mode) previously forced a temporary biallelic‑like emission at anchors. The final design always uses strict 4‑bit class‑equality at anchors and does not branch emissions. Keep parity mode off for production.

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

## Supersites: Intended Algorithm (Final Design)

Supersites collapse all split records at the same `(chr,bp)` position into a single HMM locus (the anchor). The HMM runs once per biological site and projects decisions back to member biallelic splits.

Key invariants and terminology:
- c0/c1 (immutable): The two supersite “allele classes” carried by the sample at that locus. They are 4‑bit codes: 0=REF, 1..N=ALT_i. These are snapshotted from the sample at supersite metadata build time and used only for emissions. They are never re‑derived from per‑split hap bits.
- h0/h1 (mutable, per‑epoch): The sampled classes chosen this epoch for haplotype 0 and 1, derived from the sampled lanes via the Ambiguous mask at the anchor (h0/h1 ∈ {c0,c1} for AMB; sampled from SC at MIS). After sampling, h0/h1 are projected to all member splits and do not affect emissions in subsequent epochs.
- Projection: Exactly one ALT per hap per biological site is set across member splits (or zero if class is REF), ensuring mutual exclusivity.

Supersite data:
- buildSuperSites() groups split records, chooses an anchor (global_site_id), and packs panel haplotype allele codes into a 4‑bit buffer (two codes per byte). Siblings are anchored to the same supersite id.
- Siblings never run their own DP; they are bookkeeping only (treated as MIS in HMM) because the anchor encodes the biological site.

Anchor emissions (strict class equality):
- For each anchor, build `lane_class[h]` by mapping c0/c1 through the Ambiguous mask at that locus (AMB → lanes alternate between c0 and c1; HOM → all lanes=c0).
- Donor match per lane is `donor_code == lane_class[h]` (4‑bit class equality). This single emission encodes both the multiallelic equality and the implied biallelic anchor bit (lane_class[h] == anchor_class).
- No “split/amb_mask” emission path is used at anchors.

Sampling and write‑back:
- The forward/backward algorithm builds the 8×8 per‑lane emission and transitions. sample() picks a diplotype (hap0_lane, hap1_lane).
- At an AMB anchor, derive h0/h1 from lanes via the Ambiguous mask into {c0,c1}. Equal lanes are allowed (temporary 0|0 or 1|1 at the anchor), mirroring bial behavior.
- Project h0/h1 to all member splits, setting exactly one ALT per hap where appropriate; clear others. The anchor’s HET/HOM flag remains unchanged by projection (as in the bial path).
- After `make(...)`, call `projectSupersites()`; then `H.updateHaplotypes(...); H.transposeHaplotypes_H2V(...)` for the next PBWT.

Missing anchors (multivariant imputation):
- Backward pass computes SC: `SC[offset + hap*C + c] = P(class_c | hap)` across C classes (REF+ALTs).
- sample() draws h0/h1 from SC independently (one class per hap); projection writes per‑split hap bits accordingly.

Snapshot timing:
- After building supersites (initialisation or metadata rebuild), set supersite context on each genotype and snapshot c0/c1. Use this immutable snapshot for all emissions until the next metadata rebuild.
- h0/h1 are stored per iteration and overwritten by the next sample().

PBWT / siblings:
- Siblings are masked from PBWT selection (only anchors are allowed in PBWT evaluation sites) to avoid double‑counting states.

Parity with biallelic:
- The bial path writes only hap bits at AMB and never toggles HET/HOM; supersites mirror this.
- Supersites allow equal lanes at AMB (temporary 0|0 or 1|1) to match bial’s transient behavior.
- Emission semantics at anchors are now unified (strict class equality) and do not depend on a parity toggle.

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

### Supersite Data Structures

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

### Supersite Preprocessing (`super_site_builder.*`)

Called once per MCMC iteration after PBWT state selection.

**Key Operations**:
- **Grouping**: Identify split multiallelic records at same `(chr, bp)` position
- **4-bit Encoding**: Assign codes 0=REF, 1-15=ALT1-15 per conditioning haplotype (first matching ALT wins)
- **Packing**: Store 2 codes per byte at `panel_offset`
- **Chunking**: Sites with >15 alleles split into multiple SuperSite records
- **Lookup Tables**: Build `locus_to_super_idx` and flatten `super_site_var_index`

### HMM Emissions at Supersite Anchors

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

**Sample Classification** (anchor view):
```cpp
uint8_t c0 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap0);
uint8_t c1 = getSampleSuperSiteAlleleCode(G, ss, super_site_var_index, hap1);

if (c0 == MISSING && c1 == MISSING) → MIS
else if (c0 == c1)                  → HOM (even if both are ALT5)
else                                → AMB
```

**Emission Computation** (strict class equality):
- Build lane_class[h] from immutable c0/c1 using the Ambiguous mask at the anchor.
- Emission per lane: 1.0 if `donor_code == lane_class[h]`, else `ed/ee`.
- MIS: use 1.0 (uninformative).

**Implementation**: Cache donor codes; compute class-equality mask per donor×lane; feed per‑lane emissions to the same Li–Stephens transition math as bial.

### Backward Pass: Multivariant Imputation

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

### Sampling: Lane Selection → h0/h1 → Projection

**Missing anchors**: use SC to sample classes, then project.

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

**Result**: Each haplotype has exactly 0 or 1 ALT across all member splits; anchor HET/HOM flag is not rewritten.

### Supersite-Specific Optimizations & Error Handling

**Performance**:
- 4-bit encoding: 50% memory reduction vs byte-per-code
- Anchor gating: DP runs once per biological site (not per split)
- Donor code caching: 10-20% speedup in SS-heavy windows

**Error Handling**: Same underflow detection and precision fallback as biallelic (see phase_common section).

**Design**: Minimal divergence from biallelic — reuse HMM arrays, transition math, SIMD patterns. Zero overhead when `--enable-supersites` is off.

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

### Test Output Format Standard

To ensure compatibility with `tests/run_tests.sh` automated test runner, all test binaries MUST follow this output format:

**Final success line (REQUIRED):**
```cpp
std::cout << "test_binary_name: PASS" << std::endl;
return 0;
```

**Example:**
```cpp
std::cout << "test_packing_format_diagnostic: PASS" << std::endl;
std::cout << "test_supersite_combine: PASS" << std::endl;
```

**For tests using TEST_RUN/TEST_EXIT macros:**
The macros handle this automatically - no manual PASS line needed.

**Words to AVOID in test output:**
- ~~"error"~~ in descriptive context (use "emission", "mismatch", "deviation" instead)
  - ✗ BAD: `Mismatch error rate: 0.01`
  - ✓ GOOD: `Mismatch emission rate: 0.01`
  - ✗ BAD: `However, if errors are observed in output`
  - ✓ GOOD: `However, if issues are observed in output`
- "failed" or "FAIL" outside of a proper test result line
- "✗" (cross mark) outside of a proper test failure line

**Allowed patterns that won't trigger false positives:**
- "error" is fine in proper context like variable names, technical terms
- The test runner specifically looks for failure indicators: `✗`, `failed`, `[FAIL]`
- It does NOT flag generic uses of "error" anymore (fixed Nov 2025)

**Rationale:** The test runner uses pattern matching to detect test results. Explicit failure markers (`✗`, `failed`, `[FAIL]`) indicate test failures, while the final `test_name: PASS` line indicates success.


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

This section tracks the development status of the supersite feature. While most individual components have been coded, critical bugs remain that prevent the feature from being production-ready.

**Feature Implementation Status:**
- ✅ **Implemented**: Data structures (compute_job.SC, anchor_has_missing, thread-local offset management)
- ✅ **Implemented**: Backward pass multivariant computation (IMPUTE_SUPERSITE_MULTIVARIATE)
- ✅ **Implemented**: Genotype projection (sampling from multivariant posteriors and projecting to splits)
- ✅ **Implemented**: Integration (phaser_algorithm passes context to genotype objects)
- ✅ **Implemented**: Dual precision support (compiles and runs in both single and double precision)
- ✅ **Implemented**: AMB lane semantics (correctly uses G.Ambiguous mask for per-lane phasing hypotheses)
- ✅ **Implemented**: Post-HMM projection approach (`projectSupersites` method)
- ✅ **Implemented**: Configurable emission semantics (`ss_anchor_split_emissions` parameter)
- 🚧 **Partially Implemented**: K inflation resolution. While the root cause within the HMM emission logic is fixed and passes unit tests, the issue persists in full integration runs.
- ✅ **Implemented**: Comprehensive test suite (unit tests for emission validation, parity, and K-inflation exist).
- ✅ **Implemented**: MAC-aware supersite builder. `buildSuperSites()` now accepts the PBWT MAC threshold, drops any sibling split whose MAC falls below that cutoff, and only emits a supersite when at least two ALTs survive. If the historical anchor is filtered out, the next surviving split becomes the anchor automatically; if ≤1 ALT survives, the locus stays biallelic/phase_rare just like the rest of the pipeline.

**Open Critical Issues:**
- ❗️ **Finalization Crash**: A reproducible out-of-bounds crash occurs at the end of a run when supersites are enabled, making the feature unusable.
- ❗️ **K-Inflation at Scale**: The conditioning set size (`K`) continues to grow in full integration tests, indicating a systemic issue that is not captured by unit tests.
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

**CLI Flags:**
- `--enable-supersites`: Enable multiallelic phasing support (default: off).

- `--ss-anchor-split-emissions`: A flag to control HMM emission semantics at supersite anchors. This flag allows switching between two different strategies for handling multiallelic variants.
  - **`false` (Default Mode)**: Represents the original and intended design for the supersite feature. It uses "strict" multiallelic semantics, where all equality checking for a supersite is handled at the anchor. During the HMM calculation, each alternate allele (e.g., ALT1, ALT2, ALT3) is treated as a distinct state using its 4-bit code. This provides the most accurate model of the underlying biology. After the HMM samples a phase, a separate function (`projectSupersites()`) is called to project the chosen multiallelic state back onto the correct biallelic split representation. This is the recommended mode for all production runs.
  - **`true` (Parity Mode)**: This mode was created to solve a specific development problem: how to validate the new, complex supersite logic against the existing, trusted biallelic algorithm. It forces the HMM to use simplified "biallelic" semantics *inside* the HMM calculation. For a given supersite anchor, the logic is reduced to a binary question: "Does the donor haplotype carry this anchor's specific ALT allele, or not?" All other alternate alleles at that site are treated as if they were the reference allele. This allowed for direct, apples-to-apples comparisons during testing, but is less accurate and is now retained only for specific debugging or validation scenarios.

**Environment Variables:**
- `SHAPEIT5_TEST_TRACE=1`: Enable verbose per-locus forward/backward TSV traces to `tests/out/`.
- `SHAPEIT5_DEBUG_UNDERFLOW=1`: Enable HMM underflow diagnostics to `logs/underflow.tsv`.

**Key Implementation Files:**
- `phase_common/src/models/haplotype_segment_single.{h,cpp}`: Float HMM core + supersite integration
- `phase_common/src/models/haplotype_segment_double.{h,cpp}`: Double HMM core + supersite integration
- `phase_common/src/objects/super_site_builder.{h,cpp}`: Multi-allelic grouping + panel encoding
- `phase_common/src/phaser/phaser_algorithm.cpp`: Window segmentation, threading, precision fallback
- `tests/makefile`: Unit test build rules, external object linking


## Side note re: compilation 
- This is a many cpu machine, and I recommend making with -j24 argument

## Current Status and Known Issues

The supersite extension is functionally complete in its core HMM logic and passes a comprehensive suite of unit tests. However, two critical issues have been identified during integration testing and full runs.

**1. Finalization Crash (vector out-of-bounds)**
- **Symptom**: A reproducible `std::vector` out-of-bounds assertion failure occurs during the final `G.solve()` call, specifically within `genotype::make(DipSampled)`.
- **Root Cause**: The issue is believed to stem from a mismatch between the genotype segment graph's variant counts (`Lengths[s]`) and the total number of variants (`n_variants`). Supersite sibling variants may be incorrectly counted during the initial `genotype::build` step, causing the total length of segments to exceed the bounds of the `Variants` array.
- **Status**: 🚧 **Not yet fixed.** This is a high-priority bug that prevents successful completion of runs with supersites enabled.

**2. K-Inflation in Integration Tests**
- **Observation**: While unit tests confirm that the HMM emission logic is correct and does not cause K-inflation, larger integration tests show that the conditioning set size (`K`) still steadily increases across MCMC iterations.
- **Root Cause**: The discrepancy suggests that while the core HMM algorithm is sound, there may be an issue in an upstream or downstream process that only manifests at scale. Suspects include the `PBWT` state selection or the `Haplotype Update` step, where subtle corruption could be introduced that is not captured by narrow unit tests.
- **Status**: 🚧 **Partially resolved.** The initial root cause in the HMM emission logic was fixed (see Implementation History), but the persistence of the issue in larger runs indicates a more complex, systemic problem.

---

## Implementation History and Bug Retrospective

The development of the supersite feature involved overcoming several significant technical challenges to ensure correctness and performance.

#### Supersite ↔ Biallelic Parity
A major focus was ensuring that the supersite HMM implementation was mathematically equivalent to the original biallelic implementation, especially at anchor loci. Early discrepancies were traced to several root causes:
- **Sibling Transitions**: Sibling variants (non-anchor members of a supersite) were incorrectly introducing extra recombination transitions. The fix was to make sibling loci true no-ops in the HMM, ensuring they do not advance the locus and have a recombination probability of zero.
- **Emission Semantics**: The original supersite emission logic used a direct 4-bit class comparison, which was not equivalent to the biallelic path's use of separate emission vectors for REF and ALT alleles. This was resolved by refactoring the supersite emission code to replicate the biallelic logic, ensuring identical emission patterns.
- **Race Conditions**: An early implementation had a race condition where multiple threads would write to a shared `class_prob_offset` field, causing memory corruption. This was fixed by moving to a thread-local storage model for this data.

These fixes were validated with a comprehensive test suite that now confirms numerical parity in unit tests.

---

## Next Phase: Advanced Optimization (Optional)

# **OPTIONAL ENHANCEMENT**: Column Marginal Preservation at Segment Boundaries

## Task: Preserve and use column-side marginals at segment boundaries in `phase_common`

### Goal
At genotype-graph segment boundaries in `phase_common`, stop broadcasting only the row marginal (`probSumK`) when seeding the next segment. Instead, also carry the column marginal (`probSumH`, 8 SIMD lanes) from the previous segment and seed via a row×column outer product with recombination mixing, then apply the (two-class) supersite emission mask `{a,b}` built from 4-bit allele IDs.

### Why
Broadcasting row marginals flattens the second-chain prior (column side) and erases per-lane class skew created by the previous segment—this is harmless-ish for biallelic ref|alt but can be materially wrong for supersite hets like alt1|alt2. Carrying `probSumH` fixes this with negligible overhead and unchanged vectorization.

---

## Supersite expansion parity status
- `tests/bin/test_supersite_expansion_epochs` (no missing): passing after segment cursor clamp + SC guard tightening.
- `tests/bin/test_supersite_expansion_epochs_w_missing`: still failing in repeat_factor=2. Symptom: burn1 anchor mismatch at bial locus 1 (bial 1|1 vs supersite 0|0) when two missing supersites are present. Guard assertions fixed and remain active; SC offsets are correct (ss_idx=1 → 1, ss_idx=6 → 25 for C=3) but SC for the first missing anchor is biased toward REF (e.g., hap0 0.30/0.70/0), while the second anchor shows ~0.82/0.18/0. Need to trace Alpha/Beta → SC mapping for multiple missing anchors in one window (check `missing_index_by_locus`, `AlphaMissing`, and donor class codes) to recover parity.

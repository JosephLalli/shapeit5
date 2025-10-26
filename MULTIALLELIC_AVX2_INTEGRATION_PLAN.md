# Multiallelic Super-Site Integration
## AGENTS.md

## Objective
Add native support for multiallelic loci to SHAPEIT5 (starting with `phase_common`) without post-processing. The solution must ensure that, at any given genomic position, each haplotype carries at most one alternate allele. This constraint must be enforced *inside* the HMM (forward/backward, emission likelihoods), not bolted on afterward.

We will do this by collapsing all biallelic-split variants that share the same (chr, bp) into a single "super-site" before HMM phasing. The HMM will then operate on one site per genomic coordinate instead of multiple split records.

Do **not** modify the PBWT internals to support multi-allele alphabets. PBWT remains binary; we will skip PBWT updates at super-sites.

---

## Key Concepts

### Super-site
A super-site is a single logical variant site representing all alternative alleles observed at a given (chr, bp).

- Allele alphabet for a super-site:
  - `0` = reference
  - `1..M` = the M distinct alternate alleles observed at that position
- Each haplotype is assigned exactly one code in `0..M`.
  - This enforces "≤1 ALT per haplotype at this position" by construction.
- A diploid sample genotype is expressed as a pair of codes for the two haplotypes:
  - `0/0`
  - `0/k`
  - `k/k`
  - Never `k/ℓ` with `k != ℓ` on the same haplotype.

We will phase super-sites in the HMM as single positions with multi-valued alleles, instead of multiple biallelic rows.

---

## High-Level Plan

1. **Preprocess input variants into super-sites**
   - Detect all sets of split biallelic records with the same `chr` and `bp`.
   - Collapse each set into one super-site.
   - Assign each unique ALT allele at that coordinate a compact integer code `1..M`.
   - For every reference haplotype in the conditioning panel, assign that haplotype exactly one code at this position (0 for reference or the appropriate ALT code).
   - Store these per-haplotype codes in a compact, bit-packed form.

2. **Expose a unified per-site accessor**
   - Create an accessor that, given a site ID and a set of conditioning haplotypes, returns an array of small allele codes (0..M) for those haplotypes.
   - For normal biallelic sites, this returns 0/1 (current behavior).
   - For super-sites, this returns one code per haplotype (0..M) by unpacking the bit-packed representation.

3. **Modify emission calculation in `phase_common`**
   - Update the AVX2 forward/backward emission kernel:
     - Detect if the current HMM site is a super-site.
     - For normal sites: keep existing 0/1 mask blend.
     - For super-sites: compute emission probability by comparing the two copied haplotype codes against the sample’s observed genotype at that site.
   - This guarantees we never emit biologically illegal combinations for that coordinate.

4. **Skip PBWT updates at super-sites**
   - PBWT in SHAPEIT5 is binary. Do not generalize PBWT to multi-allele.
   - For super-sites, we do not update PBWT / do not use the site to reshuffle conditioning haplotypes.
   - For normal sites, PBWT works unchanged.

5. **Do not modify transition logic**
   - The HMM transition model (recombination, copying model, MCMC schedule, etc.) remains unchanged.
   - We still advance one genomic “site” at a time. The only difference is that “site” may now be a super-site instead of a single biallelic record.
   - This avoids touching the hardest-to-validate math.

---

## Requirements / Invariants

### Biological correctness
- At any (chr, bp), output haplotypes must never show two distinct ALT alleles present on the same haplotype.
- A diploid genotype at a super-site must be consistent with:
  - hom-ref
  - het for one specific ALT allele
  - hom-alt for one specific ALT allele
- Genotypes like “ALT1 on hapA, ALT2 on hapA” must be impossible under the model, not just cleaned up afterward.

### Backward compatibility / safety
- If there are no multiallelic positions in the input, output must be bit-for-bit identical to current SHAPEIT5.
- Runtime slowdown in realistic data should be within single-digit %.
- Memory increase should be within ~10% for the relevant panel storage structures; less is better.
- PBWT logic in SHAPEIT5 must continue to work as-is for biallelic sites.

### Locality / isolation
- All code changes should be isolated to:
  1. variant preprocessing / haplotype panel preparation,
  2. new accessor for per-site haplotype allele codes,
  3. emission kernel changes, and
  4. a one-line guard around PBWT update.
- PBWT core implementation should not be edited.
- HMM transition code should not be edited.

---

## Implementation Steps (do these in order)

### 1. Super-site detection and construction
Add a preprocessing pass that runs before phasing for each chunk.

Input:
- Linearized list of per-site variant records after basic filtering (what `phase_common` currently consumes).
- For each record, we know: `chr`, `bp`, REF, ALT, per-haplotype allele presence (0/1 in the panel), per-sample genotype likelihoods or hardcalls.

Task:
- Group variant records by (`chr`, `bp`).
- If a group has only one ALT, keep it as a normal/biallelic site.
- If a group has ≥2 ALTs, create a super-site:
  - Assign ALT codes:
    - `code 1` → ALT1
    - `code 2` → ALT2
    - ...
  - For each haplotype in the reference/conditioning panel:
    - Determine which allele that haplotype actually carries at this position.
    - Assign exactly one code (0 for REF or k for ALTk).
    - If the haplotype appears to carry >1 ALT due to normalization artifacts, resolve this by:
      - picking a canonical ALT allele (longer allele / best sequence alignment),
      - or (fallback) treating ambiguous haplotypes as missing/REF for this site.
      - TODO: implement deterministic canonicalization policy and document it.

Output:
- An array `supersites[]`, each with:
  - `global_site_id`
  - `chr`, `bp`
  - `n_alt` (M)
  - `bitwidth = ceil(log2(M+1))` (typically 2–3 bits)
  - `panel_offset` (offset into packed haplotype code storage)
  - `is_super_site` flag

- A packed haplotype-code matrix:
  - For each super-site:
    - A tightly packed bitstream holding the allele code (0..M) for every haplotype in the panel.
    - Store in a single contiguous buffer, and remember `panel_offset` + `bitwidth` for decode.

Also:
- A mapping from original per-record indices → super-site index (for downstream lookups and for writing output VCF/BCF).

### 2. Per-site accessor for haplotype codes
Add a new accessor (to live alongside the existing haplotype panel access code, e.g. `Hvar` or equivalent):

```cpp
// Pseudocode interface
void getHaplotypeAlleleCodesForSite(
    uint32_t site_id,
    const uint32_t* cond_idx,     // indices of K conditioning haplotypes
    int K,
    uint8_t* out_codesA,          // allele codes for first haplotype copy
    uint8_t* out_codesB           // allele codes for second haplotype copy
);

Details:

For normal biallelic sites:
- Return codes 0 or 1 using the same fast bit-extract logic SHAPEIT5 already uses today.

For super-sites:
- Use `panel_offset` and `bitwidth` to unpack the compressed codes (0..M) for just those K haplotypes.
- Write them into `out_codesA` and `out_codesB` as aligned 8- or 16-lane blocks so the AVX load in the emission kernel is contiguous.
- This unpack can happen per HMM step; super-sites are expected to be a minority of sites.
- Keep this logic self-contained. The rest of SHAPEIT5 should not know bitwidth/packing details.

### 3. Emission kernel update in `phase_common`

Modify the AVX2 forward/backward code where emissions are currently computed using 0/1 logic.

Current behavior (simplified):
- For each site:
  - Get 0/1 “does haplotype carry ALT?” mask for K conditioning haplotypes.
  - Blend two emission values using that mask.

New behavior:
- Query `is_super_site[site_id]`.

If `false` (normal site):
- Run the existing code path unchanged.

If `true` (super-site):
- Call `getHaplotypeAlleleCodesForSite()` to get two arrays:
  - `A[i]` = allele code for haplotype A being copied at state i
  - `B[i]` = allele code for haplotype B being copied at state i
  - (0 = REF, k>0 = ALTk)
- Load those arrays into AVX2 registers as bytes/ints.
- Compute emission probability for the sample’s observed genotype at this site:

Rules:
- hom-ref (`0/0`): require `A==0 && B==0`
- het for ALT k (`0/k`): require `(A==0 && B==k) || (A==k && B==0)`
- hom-alt k (`k/k`): require `A==k && B==k`

Vectorize:
- Broadcast the sample’s `(g0,g1)` for this site.
- Compute equality masks (byte-wise compares) and combine with AND/OR to produce a boolean “match” mask per state.
- Convert that mask to an emission probability:
  - If using hard calls: match → high prob, mismatch → low prob (floor).
  - If using genotype likelihoods (GLs): use `(min(A,B), max(A,B))` as lookup key into a tiny per-site LUT of `P(data | genotype)`.

- Multiply these emission probabilities into the forward/backward DP the same way the current code multiplies `_emiss`.

This enforces that the HMM never prefers diploid configurations that violate the biology (e.g. it will never favor a state where haplotype A looks like “ALT2 and ALT5 simultaneously,” because A only ever has one code in the first place).

### 4. PBWT update guard

In the PBWT update logic (where SHAPEIT5 updates conditioning haplotype order / divergence using the current site’s alleles), add:

```cpp
if (!is_super_site[site_id]) {
    // existing PBWT update code
}

### Rationale

PBWT in SHAPEIT5 is binary (0/1). We are **not** teaching it multi-alphabet behavior here.

Skipping PBWT updates at super-sites means we are not using those sites to refine conditioning sets.  
This slightly weakens local conditioning in repetitive or STR-heavy regions, but keeps PBWT code untouched.

---

### 5. Output handling / writing phased VCF/BCF

When writing output:

- **For normal biallelic sites:** unchanged.

- **For super-sites:**
  - Generate phased GT fields that assign exactly one ALT per haplotype.
  - Reconstruct per-haplotype allele choice from the selected `(A,B)` codes in the final phased haplotypes.
  - Write back into split-biallelic representation if the final output is still expected to be biallelic-per-line.
  - **IMPORTANT:** The output re-splitting must never assign two different ALTs to the same phased haplotype at the same bp.

This step ensures the final VCF remains compatible with downstream tools that expect biallelic records, while preserving biological correctness.

---

### Validation / Testing Tasks

**No-multiallelic regression test**
- **Input:** dataset with only simple biallelic SNPs.  
- **Expected:** bit-for-bit identical output vs current SHAPEIT5. Runtime and RSS within noise.

**Unit test: super-site construction**
- Given several split records at same bp with different ALTs, confirm:
  - They collapse into one super-site.
  - Codes are assigned consistently to haplotypes.
  - No haplotype is assigned >1 ALT code.

**Unit test: emission kernel**
- For toy `K=8` conditioning haplotypes, hand-construct:
  - hom-ref, het ALT k, hom ALT k sample genotypes.
- Check that vectorized emission matches scalar reference logic.

**Integration test: phasing correctness**
- Trio or simulated truth set with a known multiallelic locus.
- Confirm we never phase ALT1 on paternal haplotype and ALT2 on paternal haplotype at the same coordinate.
- Compare switch error vs baseline SHAPEIT5 (should improve or remain unchanged at super-sites).

**Performance test**
- Phase a realistic chunk (e.g. 5–10 Mb) with a few % multiallelic sites.
- Measure:
  - wallclock runtime
  - peak RSS
- **Overhead targets:**
  - Runtime overhead: single-digit percent.
  - Memory overhead: <~10%.

---

### Non-Goals (Out of Scope)

- Do **not** implement multi-alphabet PBWT. PBWT must stay binary.
- Do **not** modify transition probability code, recombination maps, or MCMC schedule.
- Do **not** expand HMM state space to enumerate allele combinations.  
  State count must remain driven by the conditioning haplotypes.
- Do **not** attempt to model >2 different ALT alleles on the same haplotype at a site.  
  That is biologically forbidden and must remain unreachable.

---

### Deliverable Definition

After implementation:

- SHAPEIT5 (`phase_common`) can take input where multiple biallelic records share the same genomic position.  
- Internally it will:
  - Collapse them into one super-site.
  - Phase that super-site as a single HMM position using n-ary allele codes.
  - Guarantee each haplotype is assigned at most one ALT at that position.

- Output haplotypes will **never** contain two distinct ALT alleles on the same haplotype at the same coordinate.
- PBWT code remains unchanged except for a guard skipping updates at super-sites.
- Runtime and memory stay practical at UKB-scale.

This defines the **target behavior** for multiallelic support.  
Future work may explore multi-alphabet PBWT to allow multiallelic sites to influence conditioning haplotype selection,  
but that is explicitly **not required** for this milestone.

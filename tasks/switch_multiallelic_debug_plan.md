# Switch Multiallelic Debug Plan

## Goal
Identify why `switch` reports ~20% error rates on multiallelic files that should match the reference (aside from a few missing variants).

## Minimal Repro
- Run on a tiny region (100-200 sites).
- Keep the same truth/estimation files as the failing case.
- Ensure both files are indexed and queried with the same region string.

## Step 1: Per-Site Diagnostics
- For each synced record, log:
  - `pos`, `ID`, `n_allele` for truth/estimation/frequency
  - `REF/ALT` strings for each file
  - counts of missing calls in truth and estimation
  - counts of estimated calls after `MissingEst`/`Estimated` gating

## Step 2: Alignment Assertions
- Abort early if any of the following differ across the synced records:
  - `pos`
  - `ID` (if present)
  - full allele string (`REF` + comma-joined `ALT`)
- Log the first mismatch with record details.

## Step 3: GT Parsing Validation
- For a handful of samples per site, print:
  - raw GT allele indices from truth and estimation
  - decoded `REF/ALT` alleles for those indices
- Confirm allele indices are within bounds and that missing GTs are treated as missing.

## Step 4: Switch Logic Sanity
- Build 2-3 handcrafted multiallelic het sites with known phase.
- Verify:
  - orientation detection (match vs swapped)
  - `Checked`/`Errors` only set when orientation is known
  - `Missing`/`MissingEst`/`Estimated` gating behaves as expected

## Step 5: Iterate and Scale
- Apply fixes found in Steps 1-4.
- Re-run minimal region until error rates are sensible.
- Scale to the full region and confirm the expected low error rate.

## Notes
- Mendel logic currently skips multiallelic loci; this should not affect switch error computation but may affect Mendel-based outputs.
- MAC for multiallelic sites is computed as `min(sum(ALT_AF), 1 - sum(ALT_AF))` (or equivalent via AC/AN).

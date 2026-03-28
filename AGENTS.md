# SHAPEIT5 Coding Guide

## Scope and sync
- This file mirrors `.github/copilot-instructions.md`; update both together.
- Source of truth is the code. Treat existing markdown as context only.

## Project snapshot
- C++17/AVX2 haplotype phasing toolkit with modules: `phase_common`, `phase_rare`, `ligate`, `switch`, `simulate`, `xcftools`.
- `phase_common` is the main phasing engine; supersites are optional via `--enable-supersites`.

## Build and run
- Top-level build: `make` (runs each module makefile).
- Module builds: `make -C phase_common`, `make -C phase_rare`, etc.
- Tests: `make -C tests` to build, `make -C tests test-run` to execute selected tests.
- Runtime libs: tests and binaries expect `LD_LIBRARY_PATH=$HOME/.linuxbrew/lib:/usr/local/lib:$LD_LIBRARY_PATH`.
- Toolchain: g++ C++17 with AVX2 (`-mavx2 -mfma`), HTSlib, Boost (`program_options`, `iostreams`, `serialization`).

## Entry points and key code
- `phase_common/src/main.cpp`: CLI entry; defines `_DECLARE_TOOLBOX_HERE`.
- `phase_common/src/phaser/phaser_*`: option parsing, IO init, MCMC loop, supersite metadata rebuild.
- `phase_common/src/models/haplotype_segment_{single,double}.*`: HMM forward/backward (8-lane SIMD).
- `phase_common/src/models/site_emission_{types,adapter}.h`: unified biallelic vs supersite emission view.
- `phase_common/src/objects/genotype/*`: genotype graph, sampling, supersite projection.
- `phase_common/src/objects/super_site_builder.*` and `phase_common/src/models/super_site_accessor.h`: supersite metadata, class packing.
- `phase_common/src/containers/conditioning_set/*`: PBWT state selection (supersite-aware).

## HMM core (current behavior)
- Li-Stephens forward/backward HMM with `HAP_NUMBER=8` SIMD lanes.
- Emission mismatch multiplier is `M.ed/M.ee` (defaults `0.0001/0.9999`).
- Transition probabilities come from map distance via `getForwardTransProb` / `getBackwardTransProb`.
- Forward/backward loops consume `SiteView`/`EmitKind` from the emission adapters so biallelic and supersite anchors share the same DP path; supersite siblings are bookkeeping-only.
- Underflow handling: float path retries in double; samples track `double_precision`.

## Supersites (current implementation)
### Data model
- `SuperSite` holds anchor `global_site_id`, `n_alts` (1..255), `var_start/var_count`, `panel_offset/panel_span_bytes`, `n_classes=1+n_alts`, `rare_code_mask`.
- `locus_to_super_idx` maps variant -> supersite index; `super_site_var_index` is the flat member list.
- Panel codes are 1 byte each (0=REF, 1..n_alts=ALT1..ALTn), packed 1 per byte in `packed_allele_codes`.

### Build and refresh
- `buildSuperSites` creates one supersite per multiallelic record (n_alts > 1) and stores one-byte class codes per haplotype; anchor is the multiallelic locus.
- `rebuildSupersiteMetadata` runs at init and after each haplotype update; applies PBWT anchor mask/redirect and `M.markSuperSiteSiblings`.
- `updateSuperSiteAnchorEncoding` fixes anchor HET flags for multi-allelic hets.

### Classes and projection
- `resolveSupersiteClasses` reads observed `c0/c1` via `getSupersiteObservedGt` (canonicalized there) for downstream use.
- `snapshotSupersiteObservedGts` stores immutable `c0/c1` for emissions; `setSupersitePhasedGt` stores sampled `h0/h1`.
- `genotype::projectSupersites` re-applies sampled `h0/h1` after each `sample()`.

### Emissions and gating
- `SupersiteEmissionAdapter` builds a `SiteView` per locus; only anchors run DP, siblings are no-op bookkeeping.
- Anchor emissions use immutable `c0/c1` and the Ambiguous mask to set per-lane expected class; donor match is strict class equality.
- Supersite donors are cached in `ss_panel_matrix` (cond_idx snapshot) to keep forward/backward consistent.

### Missing and imputation
- `compute_job` marks `anchor_has_missing` by anchor MIS and allocates SC with per-supersite offsets.
- Backward pass calls `IMPUTE_SUPERSITE_MULTIVARIATE` at missing anchors; `missing_index_by_locus` maps forward missing slots.
- `genotype::make` samples one class per hap from SC (ALT-only CDF, REF implicit) and updates the anchor plus `ss_phased_gts`.

### PBWT
- `conditioning_set` masks non-anchor supersite members from evaluation and uses 8-bit class codes at anchors.
- `--no-supersite-pbwt` disables supersite-aware PBWT even when supersites are enabled.

## Indexing semantics (supersites)
- `Variants` and `Lengths` are per-variant (include siblings).
- `Ambiguous` and `Lengths_bio` are per biological site; siblings do not advance ambiguous or segment-locus cursors.
- Members are not guaranteed contiguous; always iterate via `super_site_var_index`.

## Development guardrails
- Hot HMM paths require 32-byte alignment (`aligned_vector32`) and no heap allocations.
- For supersite emissions, use `getSupersiteObservedGt` (immutable) rather than re-deriving from `Variants`.
- `_DECLARE_TOOLBOX_HERE` must appear in exactly one .cpp per binary (`phase_common/src/main.cpp`, `tests/src/test_toolbox.cpp`).

## CLI options worth remembering
- `--input`, `--region`, `--output`, `--thread`, `--map`, `--reference`, `--scaffold`, `--mcmc-iterations`.
- Supersites: `--enable-supersites`, `--no-supersite-pbwt`.
- Legacy toggles: `--revert-buffer-fix`, `--restore-legacy-min-transitions` (incompatible with `--enable-supersites`).

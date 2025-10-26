# AGENTS.md

Guidance for working on this repository with the Codex CLI agent. Focus areas: determinism, current feature state, and how to validate changes.

## Determinism & Reproducibility

- Do not pass any thread flag when you need deterministic output.
  - Using `--thread 1` (phase_rare) or `--threads 1` (some other tools/pipelines) still breaks determinism due to internal ordering and orientation effects.
  - Omit the flag entirely to run single‑threaded deterministically.
- Always set a fixed `--seed`.
  - Recommended values used in tests: `--seed 42` or `--seed 15052011` (see test scripts).

## Multiallelic One‑Allele Enforcement: Current State

- Common pipeline (`phase_common`)
  - Flags: `--enforce-oneallele`, `--oneallele-mode {transition|micro|micro-donor}`, `--oneallele-stats <file>`
  - Integrated after each sampling step and in a final pass; all three modes are production‑ready.

- Rare pipeline (`phase_rare`)
  - Flags: `--enforce-oneallele-rare`, `--oneallele-rare-mode {pp-basic|pp-enhanced|sparse-transition|sparse-micro}`, `--oneallele-rare-stats <file>`
  - Enforcement is wired into Step 3.5 after merging rare variants back onto the scaffold.
  - Stats type is unified between the phaser and the rare enforcer module.

## Testing & Validation

- Deterministic runs
  - Omit any thread flag; pass a fixed `--seed`.

- WGS PAR2 fixtures (CHM13 T2T, 1KGP)
  - Region: `chrX:153929053-154248138`
  - Map: `test/info/par2.gmap.gz`
  - Inputs: `test/wgs/target.*.1kgp_t2t.par2.bcf`
  - Scripts: `test/scripts/phase.wgs.{family,unrelated,haploid}.sh`

- Validation tools
  - `tools/check_oneallele <out.bcf>` to count per‑position violations.
  - Integration tests compare MD5 of headerless VCFs in `test/scripts/expected/`.

## Notes for Contributors

- Rare pipeline’s one‑allele enforcement is handled by `RareOneAlleleEnforcer` and uses the shared stats struct defined in the phaser class.
- If adding new enforcement metrics, extend the unified stats struct in `phase_rare/src/phaser/phaser_header.h` and use it across modules.

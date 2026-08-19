# Changelog

All notable changes to CLADE are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [0.1.0] — Unreleased

### Added
- Initial public packaging of the six-stage validation framework as `src/clade/`, with a `clade validate` CLI.
- Real, tested implementations of Stages 1–5 (structure-corrected association, clonal-recurrence, effect-direction, phylogenetic temporal ordering, matched-neighbor comparison), refactored from the original case-study analysis scripts into importable, unit-tested functions.
- `clade.provenance.phenotype` — encodes the fix for a real phenotype-file identity error caught during the case study (see `docs/reproducibility/phenotype_correction.md`).
- 26-test suite (unit + integration) using synthetic fixtures only — no dependency on the real 3,261-genome dataset.
- GitHub Actions CI (tests + lint).
- Full documentation set: framework specification, detailed analysis log, phylogenetic-alignment-sensitivity limitation writeup, open-items tracker, per-category results documentation.

### Renamed
- Public-facing project identity changed from **EpiCollapse** to **CLADE** (Clonal-Lineage-Aware Detection of Epistasis). Historical analysis logs and file names from the EpiCollapse era were preserved, not rewritten — see `docs/repository_audit.md`.

### Known limitations
- Stage 4 (temporal ordering) is not fully robust to core-genome alignment choice — see `docs/limitations/phylogenetic_alignment_sensitivity.md`. This is a real, open finding, not yet resolved.
- The `src/clade/` package has been verified end-to-end against synthetic fixtures, not yet re-run against the full real case-study dataset (the original case-study numbers, reported in both manuscripts, were computed by the pre-packaging analysis scripts, which remain the authoritative source for those specific numbers until CLADE's packaged version is cross-validated against them).

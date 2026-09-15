# Changelog

All notable changes to CLADE are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased] — implementation review

A technical/scientific review of the implementation against the manuscript
methodology. No manuscript result was recomputed and no reported number was
changed; where the code and the manuscript disagree, the discrepancy is flagged
rather than silently resolved. See `docs/review/implementation_review.md`.

### Fixed — scientific correctness

- **Missing evidence could promote a candidate to `convergent`.** The
  disposition rule returned CONVERGENT whenever no stage had *explicitly*
  failed, so a candidate with only a Stage 1 result — Stages 3, 4 and 5 never
  run — was classified as convergent across all six stages. Every stage field
  is now tri-state and convergence requires explicit positives at Stages 1, 3,
  4 and 5. A new `insufficient_evidence` disposition covers "nothing
  contradicts, but the evidence base is incomplete".
- **`unresolved` was unreachable from the CLI.** The verdict required a Stage 5
  concordance flag the CLI never computed (hard-coded to `None`). Stage 5 now
  reports effect direction (`Enriched_in_resistant`), and the DnaA-style
  disagreement pattern the manuscript reports as a headline finding is now
  reproducible end-to-end.
- **Screened-out candidates were reported as never tested.** A candidate failing
  the Fisher/FDR prevalence pre-filter received `not_individually_retested` —
  CLADE's phrase for a coverage limitation, not a result. Stage 1 bookkeeping is
  now explicit about which of "excluded", "no association", "no valid estimate"
  and "tested" applies.
- **The uncorrected screen could satisfy Stage 1.** The Fisher/FDR pre-filter is
  not structure-corrected and can now only ever produce a *negative* Stage 1
  result, never a positive one.
- **Skipped stages were recorded as negative.** `stage5_significant` was set to
  `False` when Stage 5 had not been run at all.
- **NaN genotypes were reconstructed as a third character state** in Fitch
  parsimony, and silently counted as absence in Stage 5. Non-binary states are
  now refused; genuinely uncalled data stays missing throughout.
- **A tree whose leaf labels did not match any sample ID produced a clean table
  of zero gains.** Stage 4 now requires a minimum overlap and reports it.
- **VCF no-calls became reference calls.** `build_matrix.py` coded `.` and
  `./.` as 0, asserting absence where the position simply could not be
  genotyped, and made a candidate position absent from the VCF look like
  universal absence. Both are now preserved as missing.
- Stage 3 treated a coefficient of exactly zero as "wrong direction", turning a
  null effect into explicit negative evidence.
- Stage 4's `Pct_post_resistance` is NaN when no origin was reconstructed; that
  is undefined, not a failure, and no longer contributes a negative verdict.
- Stage 5 distance ties were resolved by input column order, so the chosen
  neighbour depended on the caller's sample ordering. Ties now break on sorted
  sample ID and are counted.
- Stage 2's inner join silently dropped carriers with no lineage assignment from
  the carrier count.

### Added

- `clade.io.validation` — input validation for every entry point. Malformed
  input now fails with a named cause and exit status 2.
- `--stage1-results`: supply kinship-corrected results (e.g. pyseer output),
  the framework's stated preferred Stage 1 method, which previously had no
  route into the CLI.
- `--stage6`: feed recorded manual corroboration findings into the verdict.
- `--unique-neighbors`: one-to-one Stage 5 matching as a pseudoreplication
  sensitivity check; neighbour reuse is now reported by default.
- `--ambiguity-resolution`, and `temporal_ordering_sensitivity()`, which runs
  Stage 4 under both DELTRAN- and ACCTRAN-like tie-breaking and flags candidates
  whose verdict flips.
- Machine-readable output via `--output-tsv`; the evidence table now carries
  per-stage `pass`/`FAIL`/`not run` status, the missing stages, and the reason
  for every verdict.
- `tests/fixtures/synthetic_cohort/` — a 60-genome cohort with four planted
  candidates whose correct dispositions are known by construction, plus a
  deterministic generator. Software validation data only.
- Test suite grown from 31 to 108, including scientific controls and
  input-validation tests.

### Changed

- `classify_candidate` returns a `DispositionResult` (verdict + reasons +
  per-stage status) rather than a bare enum; `classify_candidate_enum` gives the
  old return type.
- `StageResults.stage5_concordant` is replaced by `stage5_enriched_in_resistant`,
  and `stage2_independent_origins` / `stage4_total_gains` were added.
- Missing `firthlogist` no longer aborts the run; Stage 1 is recorded as
  indeterminate with a clear explanation.
- `pyyaml` is declared as a dependency — it was imported by
  `clade.provenance.manifest` but never declared, so a fresh install could not
  import `clade.provenance` at all.

### Known limitations (unchanged, now documented and instrumented)

- Stage 5 pairs are not independent when one susceptible genome anchors many
  pairs; McNemar's p-value is anticonservative to that degree. The default
  behaviour is unchanged (it matches the published analysis) but the reuse is
  now measured and reported.
- Stage 4 remains dependent on tree topology and on ambiguity resolution.

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

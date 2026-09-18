# Changelog

All notable changes to CLADE are documented here. Format loosely follows [Keep a Changelog](https://keepachangelog.com/).

## [Unreleased] — statistical defects found by measurement, and a fully automatic pipeline

Everything below was found or verified by actually running the framework
against real data, not by inspection — each entry says what broke, how it was
caught, and how it was confirmed fixed.

### Fixed — scientific correctness

- **`firthlogist`'s `bse_` (standard error) is unreliable at this covariate
  dimensionality, not just its `pvals_`.** An earlier fix replaced the
  library's profile-likelihood p-value with a manually-computed Wald p-value
  from `coef_`/`bse_`. Measured directly against an independent, unpenalised
  likelihood-ratio fit on 37 real candidates: `bse_` varies only 1.3x across
  all of them (0.120–0.156) while the true standard error varies 163,000x
  (0.144–23,531), and the two are essentially uncorrelated (r=-0.18). A
  100-permutation null-calibration run using the Wald statistic found 16–20.5
  of 37 candidates significant under a true null, where Benjamini-Hochberg
  should give approximately zero. Stage 1 now computes its own
  likelihood-ratio test via `statsmodels.Logit`, using no `firthlogist` output
  for inference; `firthlogist`'s coefficient is used for direction only when
  the library is importable, and is no longer required — see the next entry.
- **Stage 1 could not run at all without `firthlogist` installed.** The
  coefficient still came from an unconditional `FirthLogisticRegression`
  import, so on any machine without it (Python >=3.11, where the package
  cannot build) every candidate failed with an opaque `ImportError`. Found by
  running the packaged CLI on this machine's own Python 3.12 environment.
  `_independent_fit()` is now the required computation (`statsmodels`, no
  `firthlogist` dependency); the library is used only as an optional
  refinement when present.
- **Stage 5's matched-pair test used one control genome for 1,491 of 2,492
  pairs on the real case-study cohort** — 83% of matches were exact distance
  ties, and the deterministic tie-break sent them all to the same handful of
  controls. McNemar's test assumes independent pairs; this reproduced the
  manuscript's published p=4.2e-11 exactly under the flawed method, confirming
  the mechanism. Fixed two ways: ties are now spread across every equally-near
  control rather than always the first (raised unique controls 61 -> 302 on
  the real cohort), and a new conditional-logistic-regression analysis
  (`stage5_conditional.py`, stratified by control genome via a Cox fit) is
  available as the defensible statistic where reuse remains high.
- **A single reconstructed Stage 4 origin was accepted as evidence of temporal
  ordering.** One event is trivially 0% or 100% "post-resistance" and cannot
  establish a pattern; several case-study candidates were reported at "100%
  post-resistance" on exactly one origin. `stage4_supports()` now requires
  `min_origins` (default 2) and returns not-evaluable below it.
- **`stage5_supports()` conflated "underpowered" with "contradicted".** A
  non-significant Stage 5 result — correctly signed but not significant — was
  mapped to the same explicit-contradiction state as a significantly
  wrong-signed one. Found while verifying the conditional-Stage-5 fix above by
  running the disposition code directly: a candidate with three of four
  required stages passing and a genuine non-significant (not contradictory)
  Stage 5 was coming back `REJECTED` instead of `INSUFFICIENT_EVIDENCE`. Now
  only a significant, wrong-signed result counts as a contradiction.
- **The CLI matched Stage 1's `Status` field against the literal string
  `"OK"`.** Adding engine provenance to that same string (to record whether a
  candidate's coefficient came from `firthlogist` or the independent fit)
  silently broke the match for every candidate. Moved provenance to its own
  `Engine` field; the CLI's check now uses the actual signal
  (`Firth_p is not None`), which also correctly stopped collapsing a
  quasi-separated candidate (a handful of cases in the minority cell) into
  "no estimate" the same as a completely separated one — a real, previously
  silent gap between two candidates the framework should have told apart.

### Added

- **A fully automatic pipeline mode.** `--candidates` and `--resistance_gene`
  are now optional. Omitted, `DETECT_RESISTANCE_GENE` picks the acquired-AMR
  gene closest to 50% cohort prevalence from the cohort's own AMRFinderPlus
  calls (most statistical power for every downstream stage; reported in full
  and overridable, not a validated phenotype assignment), and
  `DISCOVER_CANDIDATES` generates an unranked, genome-wide candidate list
  directly from `core.vcf` — every biallelic core-genome SNP becomes one
  candidate, deliberately making no ranking claim, unlike the project's
  original invalid MI/DCA discovery step. Verified against real data: the
  gene-detector picked a 51.3%-prevalence gene over a 76.7%-prevalence one
  from 300 real samples' real AMR calls (the point, not a bug); the
  candidate-discoverer emitted 5,783 candidates from the real 44MB
  `core.vcf`, including HTZ92_2925 M21L's own position, with no candidate
  list supplied.
- **`final_verdict.py` / `FINAL_VERDICT.md`** — the pipeline's last word,
  stated once: `CONVERGENT`, `NOMINATED FOR FOLLOW-UP (not confirmed)`, or
  `NO CANDIDATE`. The three-way rule mirrors
  `clade.classification.disposition` exactly. Verified against the real
  case-study evidence table (correctly nominates HTZ92_2925 M21L and no one
  else) and two synthetic edge cases (forced all-rejected, forced convergent).
- **`--from_matrix`** — an entry point for cohorts already past genome
  processing (genotype matrix, phenotype, lineages, tree, distance matrix
  already extracted), which is the state this project's real case-study
  analysis was actually run from. Runs the identical validation modules as
  the raw-reads path; only genome processing is skipped. Verified end to end
  against the real 3,254-genome cohort (2m20s, not stubbed): reproduces every
  published number to the decimal, including HTZ92_2925 M21L's conditional
  Stage 5 p=0.340 against a naive p=4.19e-11 that matches the manuscript
  exactly.
- **`stage5_conditional.py`, `fwer_permutation.py`, `benchmark_single_test.py`**
  as optional Nextflow stages (`--run_stage5_conditional`,
  `--run_fwer_permutation`, `--run_benchmark`) — the conditional matched-pair
  re-analysis above; a lineage-preserving permutation estimate of the whole
  six-stage screen's family-wise error rate (0.02–0.05 across 100
  permutations on the real cohort, at or below nominal, measured after the
  Stage 1 fixes above — the screen is calibrated even though its first-fitted
  Stage 1 was not); and a labelled reimplementation of a conventional
  single-test pipeline for comparison (16–17 apparent hits on the real
  cohort, several separated, zero CLADE-convergent).
- **`interactive_stage6.py`** — the one stage CLADE does not automate.
  Refuses to run without a real terminal rather than producing silent empty
  output; records y/n/unknown per candidate to `stage6.json` after every
  answer (safe to quit and `--resume`); explicit that its findings can adjust
  bookkeeping but cannot promote a `REJECTED` or `INSUFFICIENT_EVIDENCE`
  verdict to `CONVERGENT`.
- CI: a second `tests.yml` job installs with no `firthlogist` at all and
  confirms it is genuinely absent before running the suite — the previous
  single job always installed the `firth` extra and so could never have
  caught the no-firthlogist regression above. `lint.yml` now also checks
  `workflows/nextflow/bin/*.py`. `nextflow-stub-test.yml` gained two jobs:
  full auto-discovery with every optional stage enabled, and `--from_matrix`
  with every optional stage enabled — both replicated locally before being
  trusted, which caught a real "Process requirement exceeds available CPUs"
  failure the `--from_matrix` job would otherwise have hit on every run
  (that entry point bypasses `-profile test`'s resource ceiling and
  inherited the top-level 8-cpu default).

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

# CLADE

**Clonal-Lineage-Aware Detection of Epistasis** — a six-stage framework for validating candidate compensatory mutations in clonally reproducing bacterial populations.

This repository is the CLADE software package: the framework code, its CLI, and its test suite. It validates a candidate-locus list supplied to it (genotype matrix + phenotype + lineage assignments + optional phylogeny/distance matrix) — it does not claim to solve population-structure-safe candidate *discovery*.

## Overview

Antimicrobial resistance often carries a fitness cost, and resistant bacterial lineages that acquire secondary, compensatory mutations offsetting that cost are expected to be favored by selection. Screening genome collections for these mutations is confounded by clonal population structure: because bacteria reproduce clonally, a mutation that arose once in a successful resistant lineage will be statistically associated with resistance across every descendant of that lineage — whether or not it does anything compensatory at all.

## Why CLADE?

A single structure-corrected association test reduces this problem but does not solve it. In development, CLADE was validated against a real 3,261-genome *bla*OXA-23 / *Acinetobacter baumannii* cohort, and the result made the point directly: the dataset's single strongest candidate by kinship-corrected regression alone (a DnaA synonymous variant, recurring across 47 lineages, q=5.1×10⁻³) failed on two further checks — its effect ran the wrong direction, and it mostly predated resistance acquisition rather than following it. Under a single-test pipeline, this candidate would have been reported as the strongest finding. CLADE exists to catch exactly this.

## The six-stage framework

| Stage | Checks | Catches |
|---|---|---|
| 1. Structure-corrected association | Kinship-corrected regression (`pyseer`), Firth regression as fallback | Associations that vanish once population structure is accounted for |
| 2. Clonal-recurrence | Single lineage vs. many (MLST cross-tab) | "Real recurring pattern" vs. "happened once, inherited ever since" |
| 3. Effect direction | Is the association enriched (positive) in resistant genomes? | Frequently-skipped check that can disqualify a strong-looking Stage 1 result outright |
| 4. Temporal ordering | Fitch-parsimony ancestral reconstruction | Candidates that predate resistance (lineage markers, not compensation) |
| 5. Matched-neighbor comparison | McNemar's test, real patristic distance | Disagreement between population-wide and local evidence |
| 6. External corroboration | Literature / homology search (manual — not automated by this package, see `clade.validation.stage6_corroboration`) | Statistically convergent candidates lacking a known mechanism |

## Evidence convergence

A candidate is reported as **convergent** only if it passes every applicable stage. Every candidate that enters the framework is classified, not just the survivors:

| Disposition | Meaning |
|---|---|
| `convergent` | Passes every stage a valid test could be run for |
| `unresolved` | Stages genuinely disagree (e.g. population-wide evidence says no, local matched-pair evidence says yes) |
| `rejected` | Consistent negative evidence across ≥2 independent methods |
| `clonal_artifact` | Statistical signal fully explained by single-lineage concentration |
| `uninformative` | Near-fixed across the cohort regardless of any statistic |
| `not_individually_retested` | A coverage limitation, not a scientific verdict |

See `src/clade/classification/disposition.py` for the exact rule, and its docstring for the honest caveat on what it can and can't fully capture for borderline candidates.

## Installation

```
git clone <repo-url> && cd CLADE
python -m venv .venv && source .venv/bin/activate
pip install -e ".[firth]"
```

`firthlogist` (Stage 1's non-kinship fallback) requires Python <3.11 and `scikit-learn<1.6` — a real dependency conflict hit and resolved during development (`firthlogist` depends on a private scikit-learn API removed in 1.6+). Use `environment.yml` for a conda setup with this already pinned.

## Quick start

```
clade validate \
  --genotypes tests/fixtures/example_genotypes.csv \
  --phenotype tests/fixtures/example_phenotype.tsv \
  --lineages tests/fixtures/example_lineages.tsv \
  --candidates geneA geneB \
  --tree tests/fixtures/tiny_tree.nwk \
  --output evidence_table.md
```

## Input format

- **Genotypes**: CSV, `Sample_ID` + one 0/1 column per candidate locus.
- **Phenotype**: TSV, sample ID + a single 0/1 column — derive this from a genotype source with `clade.provenance.phenotype.derive_phenotype_from_amr_matrix` rather than maintaining it by hand (see below).
- **Lineages**: TSV with `sample` and `ST` columns (MLST output).
- **Tree**: Newick, optional (enables Stage 4).
- **Distances**: TSV pairwise patristic distance matrix, optional (enables Stage 5).

## Output format

A Markdown evidence table, one row per candidate, grouped by disposition — every candidate that entered CLADE appears, whatever the result.

## Reproducibility

- `pyproject.toml` / `environment.yml` / `requirements.txt` pin exact version constraints.
- 26 tests (unit + integration), synthetic fixtures only — CI never depends on a real dataset. Run locally with `pytest tests/`.
- `ruff check src/ tests/` for linting; both run in CI on every push (`.github/workflows/`).

## A note on phenotype provenance

`clade.provenance.phenotype` exists because of a real, caught error during development: a hand-maintained phenotype file was used through an earlier round of analysis before being discovered, by direct value comparison, to represent a *different* variable than its filename claimed. `derive_phenotype_from_amr_matrix` and `verify_phenotype_provenance` encode the fix: derive phenotypes programmatically from their generating source, and verify any externally-supplied phenotype file against that source before trusting it.

## Limitations

- Validated on one organism, one resistance determinant, one real dataset so far — generalization to other organisms or population structures is untested.
- No wet-lab validation has been performed on any candidate CLADE has reported as convergent; it establishes statistical/phylogenetic convergence, not causal proof.
- **Stage 4's dependency on tree topology is less bias-proof than might be assumed.** A direct comparison (SNP-only vs. full-alignment core-genome trees, built independently from the identical genome set) found the two disagree substantially in topology — not just branch length — even among high-confidence splits. This is exactly why CLADE requires convergence across multiple stages rather than trusting any one stage, including Stage 4, in isolation.
- Stage 6 (external corroboration) is a manual process by design — `clade.validation.stage6_corroboration.CorroborationRecord` gives it a structured place to be recorded, but does not automate literature or homology search.

## Citation

See [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).

---

## Repository map

| Directory | Contents |
|---|---|
| [`src/clade/`](src/clade/) | The framework: `io/`, `phylogeny/` (Fitch parsimony), `validation/` (Stages 1–6), `classification/` (disposition rule), `reporting/` (evidence-table generation), `provenance/` (phenotype derivation/verification), `cli.py` |
| [`tests/`](tests/) | `unit/` (per-module tests, hand-verified synthetic fixtures) and `integration/` (end-to-end CLI smoke test) |
| [`configs/`](configs/) | `example.yaml` (runnable, matches `tests/fixtures/`) and `case_study.yaml` (documents the real-dataset configuration; not runnable without external data) |
| [`.github/`](.github/) | CI workflows (tests, lint), issue templates, PR template |

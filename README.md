# 🧬 CLADE

**Clonal-Lineage-Aware Detection of Epistasis** — a six-stage framework for validating candidate compensatory mutations in clonally reproducing bacterial populations.

[![Tests](https://github.com/biocoderep/clade/actions/workflows/tests.yml/badge.svg)](https://github.com/biocoderep/clade/actions/workflows/tests.yml)
[![Lint](https://github.com/biocoderep/clade/actions/workflows/lint.yml/badge.svg)](https://github.com/biocoderep/clade/actions/workflows/lint.yml)
[![Nextflow pipeline](https://github.com/biocoderep/clade/actions/workflows/nextflow-stub-test.yml/badge.svg)](https://github.com/biocoderep/clade/actions/workflows/nextflow-stub-test.yml)
[![Container](https://github.com/biocoderep/clade/actions/workflows/container.yml/badge.svg)](https://github.com/biocoderep/clade/actions/workflows/container.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)
[![Nextflow DSL2](https://img.shields.io/badge/nextflow-DSL2-23aa62.svg)](workflows/nextflow/)

This repository is the CLADE software package: the framework code, its CLI, and its test suite. It validates a candidate-locus list supplied to it (genotype matrix + phenotype + lineage assignments + optional phylogeny/distance matrix) — it does not claim to solve population-structure-safe candidate *discovery*. The Nextflow pipeline can generate a candidate list automatically (every biallelic core-genome SNP, unranked) when none is supplied, but that is a documented brute-force fallback, not a statistically validated discovery method — see [`workflows/nextflow/README.md`](workflows/nextflow/README.md).

<details>
<summary>📑 <strong>Jump to a section</strong></summary>

- [Overview](#-overview)
- [Why CLADE?](#-why-clade)
- [The six-stage framework](#-the-six-stage-framework)
- [Evidence convergence](#-evidence-convergence)
- [Installation](#-installation)
- [Quick start](#-quick-start)
- [From raw reads: the Nextflow pipeline](#-from-raw-reads-the-nextflow-pipeline)
- [Input / output format](#-inputoutput-format)
- [Reproducibility](#-reproducibility)
- [A note on phenotype provenance](#-a-note-on-phenotype-provenance)
- [Limitations](#-limitations)
- [Repository map](#-repository-map)

</details>

## 🦠 Overview

Antimicrobial resistance often carries a fitness cost, and resistant bacterial lineages that acquire secondary, compensatory mutations offsetting that cost are expected to be favored by selection. Screening genome collections for these mutations is confounded by clonal population structure: because bacteria reproduce clonally, a mutation that arose once in a successful resistant lineage will be statistically associated with resistance across every descendant of that lineage — whether or not it does anything compensatory at all.

## ❓ Why CLADE?

A single structure-corrected association test reduces this problem but does not solve it.

> [!IMPORTANT]
> In development, CLADE was validated against a real 3,261-genome *bla*OXA-23 / *Acinetobacter baumannii* cohort, and the result made the point directly: the dataset's single strongest candidate by kinship-corrected regression alone (a DnaA synonymous variant, recurring across 47 lineages, q=5.1×10⁻³) failed on two further checks — its effect ran the wrong direction, and it mostly predated resistance acquisition rather than following it. Under a single-test pipeline, this candidate would have been reported as the strongest finding. **CLADE exists to catch exactly this.**

## 🧩 The six-stage framework

```mermaid
flowchart LR
    S1["1️⃣ Association<br/>structure-corrected"] --> S2["2️⃣ Recurrence<br/>clonal cross-tab"]
    S2 --> S3["3️⃣ Direction<br/>enriched in resistant?"]
    S3 --> S4["4️⃣ Temporal order<br/>Fitch parsimony"]
    S4 --> S5["5️⃣ Matched-neighbor<br/>McNemar"]
    S5 --> S6["6️⃣ Corroboration<br/>manual, literature"]

    classDef stage fill:#6366f1,stroke:#4338ca,color:#fff
    class S1,S2,S3,S4,S5,S6 stage
```

| Stage | Checks | Catches |
|---|---|---|
| 1. Structure-corrected association | Kinship-corrected regression (`pyseer`, supplied via `--stage1-results`), Firth regression as fallback | Associations that vanish once population structure is accounted for |
| 2. Clonal-recurrence | Single lineage vs. many (MLST cross-tab) | "Real recurring pattern" vs. "happened once, inherited ever since" |
| 3. Effect direction | Is the association enriched (positive) in resistant genomes? | Frequently-skipped check that can disqualify a strong-looking Stage 1 result outright |
| 4. Temporal ordering | Fitch-parsimony ancestral reconstruction | Candidates that predate resistance (lineage markers, not compensation) |
| 5. Matched-neighbor comparison | McNemar's test, real patristic distance | Disagreement between population-wide and local evidence |
| 6. External corroboration | Literature / homology search (manual — not automated by this package, see `clade.validation.stage6_corroboration`) | Statistically convergent candidates lacking a known mechanism |

## ✅ Evidence convergence

A candidate is reported as **convergent** only if it passes every applicable stage. Every candidate that enters the framework is classified, not just the survivors:

![convergent](https://img.shields.io/badge/convergent-brightgreen)
![unresolved](https://img.shields.io/badge/unresolved-yellow)
![rejected](https://img.shields.io/badge/rejected-red)
![clonal_artifact](https://img.shields.io/badge/clonal__artifact-orange)
![uninformative](https://img.shields.io/badge/uninformative-lightgrey)
![insufficient_evidence](https://img.shields.io/badge/insufficient__evidence-blue)
![not_individually_retested](https://img.shields.io/badge/not__individually__retested-inactive)

| Disposition | Meaning |
|---|---|
| `convergent` | **Explicit positive** result at every stage required for convergence (1, 3, 4, 5) |
| `unresolved` | Stages genuinely disagree (local matched-pair evidence says yes, population-wide or whole-tree evidence says no) |
| `rejected` | An explicit negative result at one or more stages |
| `clonal_artifact` | Signal concentrated in one lineage with no independent recurrence, plus a contradicting stage |
| `uninformative` | Near-fixed across the cohort regardless of any statistic |
| `insufficient_evidence` | No stage contradicts, but a required stage produced no usable result |
| `not_individually_retested` | Stage 1 was never run — a coverage limitation, not a scientific verdict |

> [!NOTE]
> **Missing evidence is never treated as support.** Every stage is tri-state (`pass` / `FAIL` / `not run`), and a candidate cannot reach `convergent` because a stage was skipped, a model failed to converge, or a genotype was uncalled. That distinction is what `insufficient_evidence` exists to make: it separates "tested and survived" from "never properly tested". The evidence table names the missing stages for every candidate.

See `src/clade/classification/disposition.py` for the exact rule, and its docstring for the honest caveat on what it can and can't fully capture for borderline candidates.

## 📦 Installation

```bash
git clone https://github.com/biocoderep/clade.git && cd clade
python -m venv .venv && source .venv/bin/activate
pip install -e ".[firth]"
```

> [!WARNING]
> `firthlogist` (Stage 1's non-kinship fallback) requires Python <3.11 and `scikit-learn<1.6` — a real dependency conflict hit and resolved during development (`firthlogist` depends on a private scikit-learn API removed in 1.6+). Use `environment.yml` for a conda setup with this already pinned.

## 🚀 Quick start

<details>
<summary>Run CLADE on a 60-genome synthetic cohort with known-correct answers</summary>

```bash
F=tests/fixtures/synthetic_cohort
clade validate \
  --genotypes      $F/genotypes.csv \
  --phenotype      $F/phenotype.tsv \
  --lineages       $F/lineages.tsv \
  --candidates     true_compensatory clonal_marker wrong_direction near_fixed \
  --tree           $F/tree.nwk \
  --distances      $F/distances.tsv \
  --stage1-results $F/stage1_results.tsv \
  --stage6         $F/stage6.json \
  --output         evidence_table.md
```

This runs all five automated stages on a 60-genome synthetic cohort whose
correct answers are known by construction, and should reproduce exactly:

| Candidate | Disposition | Why |
|---|---|---|
| `true_compensatory` | `convergent` | independent origins in several lineages, all post-resistance |
| `clonal_marker` | `rejected` | one ancestral origin; the association is absorbed by the lineage covariate |
| `wrong_direction` | `rejected` | depleted among resistant genomes (Stage 3 failure) |
| `near_fixed` | `uninformative` | present in almost every genome |

**That cohort is software-validation data, not biological evidence** — see
`tests/fixtures/synthetic_cohort/README.md`. Regenerate it deterministically
with `python tests/fixtures/make_synthetic_cohort.py`.

Omit `--tree` and `--distances` and the same run yields `insufficient_evidence`
rather than `convergent`, which is the point: Stages 4 and 5 are recorded as
*not run*, not as passed.

</details>

## 🧵 From raw reads: the Nextflow pipeline

The CLI above validates a candidate list you already have inputs for. To go
from raw sequencing reads to an evidence table in one command:

```bash
nextflow run biocoderep/clade -profile docker \
    --samplesheet     samples.csv \
    --reference       reference.fa \
    --candidates      candidates.csv \
    --resistance_gene blaOXA-23 \
    --outdir          results
```

The pipeline runs read QC, trimming, reference-based variant calling, assembly,
AMR and MLST typing, core-genome alignment, phylogeny, the patristic distance
matrix, and then CLADE Stages 1–5. Every tool is pinned to an exact version and
available as a container; `results/pipeline_info/software_versions.yml` records
what actually ran.

<details>
<summary>⚡ Verify the wiring in five seconds — no tools, data or network required</summary>

```bash
cd workflows/nextflow && nextflow run main.nf -profile test -stub-run
```

</details>

Full documentation: [`workflows/nextflow/README.md`](workflows/nextflow/README.md).

## 📄 Input/output format

- **Genotypes**: CSV, `Sample_ID` + one 0/1 column per candidate locus.
- **Phenotype**: TSV, sample ID + a single 0/1 column — derive this from a genotype source with `clade.provenance.phenotype.derive_phenotype_from_amr_matrix` rather than maintaining it by hand (see below).
- **Lineages**: TSV with `sample` and `ST` columns (MLST output).
- **Tree**: Newick, optional (enables Stage 4).
- **Distances**: TSV pairwise patristic distance matrix, optional (enables Stage 5).
- **Output**: a Markdown evidence table, one row per candidate, grouped by disposition — every candidate that entered CLADE appears, whatever the result.

## 🔁 Reproducibility

<details open>
<summary>How correctness is enforced, end to end</summary>

- `pyproject.toml` / `environment.yml` / `requirements.txt` pin exact version constraints.
- 108 tests (unit + integration), synthetic fixtures only — CI never depends on a real dataset. Run locally with `pytest tests/`.
- The tests that matter most are the **scientific controls**
  (`tests/unit/test_scientific_controls.py`) and the planted-truth end-to-end
  test (`tests/integration/test_cli_smoke.py`): small datasets with known
  correct answers, so a refactor that breaks the science fails the build even
  if the code still runs.
- `ruff check src/ tests/` for linting; both run in CI on every push (`.github/workflows/`).
- The Nextflow pipeline pins every external tool to an exact version, each
  available as a BioContainer whose tag is verified to exist. `-profile conda`
  and `-profile docker` install the same versions. Nothing resolves to `latest`.
- Every pipeline run writes `software_versions.yml`, collected from the
  processes themselves rather than from documentation, plus an execution
  report, trace and DAG under `pipeline_info/`.
- CI stub-runs the entire pipeline DAG on every push, against both the oldest
  supported Nextflow version and the current release, and checks `-resume`
  caches rather than recomputes.
- All inputs are validated before use (`clade.io.validation`). Malformed input
  fails with a named cause and a non-zero exit status rather than being
  silently reinterpreted — non-binary genotype codes, duplicated sample IDs,
  tree labels that do not match the sample IDs, non-square distance matrices
  and single-class phenotypes are all refused.
- Exit codes: `0` completed, `1` usage error, `2` invalid input, `3` missing
  optional dependency.

</details>

## 🕵️ A note on phenotype provenance

> [!IMPORTANT]
> `clade.provenance.phenotype` exists because of a real, caught error during development: a hand-maintained phenotype file was used through an earlier round of analysis before being discovered, by direct value comparison, to represent a *different* variable than its filename claimed. `derive_phenotype_from_amr_matrix` and `verify_phenotype_provenance` encode the fix: derive phenotypes programmatically from their generating source, and verify any externally-supplied phenotype file against that source before trusting it.

## ⚠️ Limitations

<details open>
<summary>Read this before trusting a convergent call</summary>

- Validated on one organism, one resistance determinant, one real dataset so far — generalization to other organisms or population structures is untested.
- No wet-lab validation has been performed on any candidate CLADE has reported as convergent; it establishes statistical/phylogenetic convergence, not causal proof.
- **Stage 4's dependency on tree topology is less bias-proof than might be assumed.** A direct comparison (SNP-only vs. full-alignment core-genome trees, built independently from the identical genome set) found the two disagree substantially in topology — not just branch length — even among high-confidence splits. This is exactly why CLADE requires convergence across multiple stages rather than trusting any one stage, including Stage 4, in isolation.
- **Stage 5 pairs are not independent.** Nothing prevents one susceptible genome
  from being the nearest neighbour of many resistant ones, but McNemar's test
  assumes independent pairs, so its p-value is anticonservative to the degree
  that reuse occurs. The default (reuse permitted) matches the published
  analysis; CLADE now *reports* the reuse (`Max_neighbor_reuse`,
  `N_unique_neighbors`) and offers `--unique-neighbors` for a one-to-one
  sensitivity check. This is a property of the method, not a bug that has been
  fixed — read Stage 5 p-values with the reported reuse in hand.
- **Stage 4 depends on how ambiguous ancestral states are resolved.** The
  default dates gains as late as possible (DELTRAN-like); `--ambiguity-resolution 1`
  dates them as early as possible (ACCTRAN-like), and
  `clade.validation.stage4_temporal_order.temporal_ordering_sensitivity` runs
  both and flags candidates whose verdict flips between them.
- Stage 6 (external corroboration) is a manual process by design — `clade.validation.stage6_corroboration.CorroborationRecord` gives it a structured place to be recorded, but does not automate literature or homology search.

</details>

## 📚 Citation

See [`CITATION.cff`](CITATION.cff).

## ⚖️ License

MIT — see [`LICENSE`](LICENSE).

---

## 🗺️ Repository map

| Directory | Contents |
|---|---|
| [`src/clade/`](src/clade/) | The framework: `io/`, `phylogeny/` (Fitch parsimony), `validation/` (Stages 1–6), `classification/` (disposition rule), `reporting/` (evidence-table generation), `provenance/` (phenotype derivation/verification), `cli.py` |
| [`tests/`](tests/) | `unit/` (per-module tests, hand-verified synthetic fixtures) and `integration/` (end-to-end CLI smoke test) |
| [`configs/`](configs/) | `example.yaml` (runnable, matches `tests/fixtures/`) and `case_study.yaml` (documents the real-dataset configuration; not runnable without external data) |
| [`workflows/nextflow/`](workflows/nextflow/) | End-to-end pipeline: `main.nf`, `modules/` (one per tool), `subworkflows/`, `bin/` (scripts each module invokes), `conf/` (resources, publishing, test profile), `nextflow_schema.json` |
| [`scripts/`](scripts/) | Standalone diagnostic/calibration tools, run by hand, not wired into any pipeline module — each documents a specific investigation (e.g. effective test count under candidate correlation, empirical FDR calibration) |
| [`containers/`](containers/) | `Dockerfile` for the `clade` image, built and published by CI |
| [`.github/`](.github/) | CI workflows (tests, lint, pipeline stub-run, container build), issue templates, PR template |

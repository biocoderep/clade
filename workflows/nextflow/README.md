# CLADE pipeline: raw reads → six-stage evidence table

A Nextflow DSL2 wrapper around the real genome-processing chain (fastp → Snippy/Shovill → Abricate/AMRFinderPlus/mlst → snippy-core → FastTree, reconstructed from `pipeline_orchestrator_linux_v2.py`) and CLADE's own six-stage validation CLI. Removes the manual-input-building step: raw reads in, `evidence_table.md` out.

**Stage coverage:** runs CLADE Stages 1–5 automatically. Stage 6 (external corroboration) is a deliberate manual step by design — see `clade.validation.stage6_corroboration` — this pipeline does not, and is not meant to, automate it.

## Quick start (no real tools needed)

```
cd workflows/nextflow
nextflow run main.nf -stub-run -profile test
```

This verifies the pipeline's wiring only (channels, dependencies, aggregation logic) using tiny synthetic fixtures — no fastp/Snippy/AMRFinderPlus/etc. installation, no network access, no real data. Same philosophy as CLADE's own synthetic-fixture-only Python test suite (`tests/`).

## Real run

```
nextflow run main.nf -profile conda \
  --samplesheet samplesheet.csv \
  --reference reference.fa \
  --candidates candidates.csv \
  --resistance_gene "blaOXA-23" \
  --outdir results
```

- `-profile conda` builds the genome-processing tool environment from `conda/genome_processing.yml` (fastp, FastQC, Snippy, Shovill, Abricate, AMRFinderPlus, `mlst`, FastTree, SRA Toolkit).
- `--samplesheet`: CSV with `sample_id,accession,fastq_1,fastq_2` — populate either `accession` (SRA) or `fastq_1`/`fastq_2` (local paths) per row.
- `--candidates`: CSV with `name,position` (e.g. `SecA_M21L,CP058289.1:3155360`) — genotypes are extracted from `core.vcf` at each exact position, matching how the real case study's candidates were extracted. This pipeline is agnostic to how candidates were originally proposed (see `CLADE_Framework_Specification.md` §6) — bring your own list.
- `--resistance_gene`: the AMRFinderPlus gene symbol defining your phenotype (e.g. `blaOXA-23`). The phenotype is derived programmatically from aggregated AMRFinderPlus calls via `clade.provenance.phenotype.derive_phenotype_from_amr_matrix` — not hand-maintained, for exactly the reason documented in that module's docstring (a real phenotype-identity error this project caught and fixed).

**Before running at real cohort scale**, read `docs/limitations/phylogenetic_alignment_sensitivity.md` (kept locally, not in this public repo — see the main manuscript repository). `FASTTREE` is hard-pinned to the SNP-only alignment (`core.aln`) by default: two independent attempts to build a full-alignment tree at 3,255-genome scale were killed by out-of-memory errors after 24–27 hours each (peak 497–521GB). This is a real, documented cost, not a hypothetical one.

## Pipeline stages

| Stage | What runs |
|---|---|
| Per-sample (parallel) | Download → FastQC → fastp → FastQC → Snippy (variant call) + Shovill (assembly) → Abricate/AMRFinderPlus/`mlst` |
| Aggregation | `snippy-core` → FastTree → patristic distance matrix → merge into CLADE's genotype/phenotype/lineage inputs |
| Validation | `clade validate` (Stages 1–5) |

## What's genuinely new here vs. reused unchanged

Every per-sample and aggregation command is the real, unchanged command from the original analysis (see the per-module code comments for exact provenance). Three scripts are genuinely new, written for this pipeline specifically:

- `bin/distance_matrix.py` — nothing previously computed a reusable patristic distance matrix; without it, CLADE's Stage 5 has no input and silently never runs.
- `bin/build_matrix.py` — bridges per-sample profiling tables + `core.vcf` to CLADE's expected input format. Candidate genotypes come from `core.vcf` (variant-level), not the profiling tables (gene-level) — a distinction that matters and is documented in the script's own docstring.
- `bin/run_clade.py` — the real `clade validate --candidates` CLI takes space-separated names, not a file path; this expands the candidates CSV correctly.

All three were tested against hand-verified synthetic data (not just assumed correct) before being wired into the pipeline.

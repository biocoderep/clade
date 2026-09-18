# CLADE Nextflow pipeline

Raw sequencing reads → six-stage CLADE evidence table, in one command.

This directory holds the workflow. For what CLADE *is* and why the six stages
exist, see the [repository README](../../README.md).

---

## Quick start

Check the wiring in five seconds — no tools, no data, no network:

```bash
nextflow run main.nf -profile test -stub-run
```

Run it for real:

```bash
nextflow run main.nf -profile docker \
    --samplesheet  samples.csv \
    --reference    reference.fa \
    --candidates   candidates.csv \
    --resistance_gene blaOXA-23 \
    --outdir       results
```

`nextflow run main.nf --help` prints every option.

---

## Inputs

### `--samplesheet`

A CSV. Rows may give an SRA accession, local FASTQ paths, or both kinds mixed —
a cohort of public data plus unreleased local isolates is the normal case:

```csv
sample_id,accession
SRR11589190,SRR11589190
```

```csv
sample_id,fastq_1,fastq_2
isolate_01,/data/iso01_R1.fastq.gz,/data/iso01_R2.fastq.gz
```

### `--reference`

Reference genome FASTA. **Candidate positions must be on this coordinate
system.** Mismatched coordinates fail silently — the variant is simply absent
from every sample, and the candidate is reported as untestable rather than as
an error.

### `--candidates` (optional)

```csv
name,position
htz92_2925_m21l,3155360
aminotransferase,3712044
```

CLADE *validates* a candidate list; it does not itself claim to solve
population-structure-safe candidate *discovery*. Where the list came from
does not matter — that is the point: a list from an inadequately-corrected
discovery method can still be assessed honestly here.

If `--candidates` is omitted, `DISCOVER_CANDIDATES` generates one
automatically from `core.vcf`: every biallelic, non-monomorphic,
low-missingness core-genome SNP, unranked. This is a brute-force fallback,
not a statistically validated discovery method — see
`bin/discover_candidates.py`.

### `--resistance_gene` (optional)

The AMRFinderPlus gene symbol that defines the binary phenotype (`blaOXA-23`).
The phenotype is derived programmatically from AMRFinderPlus output rather than
supplied as a curated file, because a hand-maintained phenotype file is a
well-known source of silent error.

If `--resistance_gene` is omitted, `DETECT_RESISTANCE_GENE` picks the gene
whose carriage prevalence across the cohort is closest to 50% (bounded by
`--auto_gene_min_prevalence`/`--auto_gene_max_prevalence`) — a documented
heuristic, not a validated phenotype assignment. Always check
`resistance_gene_candidates.tsv` in the output. See
`bin/detect_resistance_gene.py`.

### `--from_matrix` (alternate entry point)

Skips genome processing (alignment, variant calling, phylogeny, distances)
entirely, for when you already have those artifacts — the state this
project's own real case-study analysis was actually run from:

```bash
nextflow run main.nf --from_matrix \
  --matrix_genotypes genotypes.csv \
  --matrix_phenotype phenotype.tsv \
  --matrix_lineages lineages.tsv \
  --matrix_tree tree.nwk \
  --matrix_distances distances.tsv \
  --matrix_candidates candidates.csv   # omit for auto-discovery, as above
```

Bypasses `-profile test`'s resource ceiling, so set `--max_cpus`/`--max_memory`
explicitly if not running on a large machine.

---

## Outputs

```
results/
├── clade/
│   ├── FINAL_VERDICT.md             <- CONVERGENT / NOMINATED / NO CANDIDATE, per candidate
│   ├── clade_evidence_table.md      <- the full per-stage evidence
│   └── clade_evidence_table.tsv
├── detect_resistance_gene/          only when --resistance_gene was omitted
├── discover_candidates/             only when --candidates was omitted
├── clade_input/                     genotypes.csv, pheno.tsv, lineages.tsv
├── variants/core/                   core.aln, core.full.aln, core.vcf
├── phylogeny/                       tree.nwk, distances.tsv
├── assemblies/
├── mlst/  amrfinder/  abricate/
├── qc/                              fastqc_raw, fastqc_trimmed, fastp, multiqc
└── pipeline_info/
    ├── software_versions.yml        every tool version used in THIS run
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.html
```

`FINAL_VERDICT.md` folds every stage (including a conditional Stage 5 result,
when `--run_stage5_conditional` was set and the built-in Stage 5 was withheld
for insufficient power) into the same three-way verdict as
`clade.classification.disposition`. It is not a substitute for Stage 6: a
`NOMINATED FOR FOLLOW-UP` verdict still needs the manual corroboration step
below before it can be called `CONVERGENT`.

`software_versions.yml` is the file that makes a run reproducible. It records
the exact version of every tool that touched the data, collected from the
processes themselves rather than from documentation.

Both `core.aln` (SNP-only) and `core.full.aln` are kept. The tree is built from
the SNP-only alignment; the full alignment is retained so the ascertainment-bias
sensitivity check remains reproducible without re-running the cohort. That check
matters: on a 150-genome subsample the two alignments produced trees differing in
72% of internal edges, not merely in branch length.

---

## Profiles

| Profile | Use |
|---|---|
| `docker`, `singularity`, `apptainer`, `podman` | Containerised — **recommended** |
| `conda`, `mamba` | Where containers are not permitted |
| `test` | Tiny inputs for `-stub-run` |
| `local`, `slurm` | Executor presets |

Combine them: `-profile test,docker` or `-profile slurm,singularity`.

Containers are preferred over conda because an image is a fixed filesystem,
whereas a conda solve can drift as channels change even with versions pinned.

---

## Reproducibility

Every tool is pinned to an exact version, and each is available as a
BioContainer image whose tag was verified to exist in the registry. The conda
environments pin the same versions, so `-profile conda` and `-profile docker`
run the same software.

| Step | Tool | Version |
|---|---|---|
| Download | sra-tools | 3.4.1 |
| QC | FastQC | 0.12.1 |
| Trimming | fastp | 0.23.4 |
| Variant calling | Snippy | 4.6.0 |
| Assembly | Shovill | 1.4.2 |
| AMR (assembly) | Abricate | 1.4.0 |
| AMR (phenotype) | AMRFinderPlus | 4.0.23 |
| Typing | mlst | 2.35.0 |
| Phylogeny | FastTree | 2.2.0 |
| QC aggregation | MultiQC | 1.35 |

Nothing resolves to `latest`.

`-resume` is supported and exercised: a re-run after a change re-executes only
the affected tasks.

### What has and has not been tested

Honest scope, because this matters more than a green badge:

- **Tested:** the complete DAG runs end to end under `-stub-run` (34 tasks,
  every process, every channel join and every file-name expectation), on every
  profile, with `-resume`. CI runs this on each push.
- **Tested:** CLADE's own scientific logic, against a synthetic cohort with
  planted ground truth where the correct answer is fixed by construction
  (`pytest`, 108 tests).
- **Not tested here:** a full real-data run with real tools. The container tags
  are verified to exist and the commands are the ones this project used, but
  this pipeline has not been re-run end to end against the 3,261-genome cohort
  to confirm numerical identity with the published results. Treat the original
  analysis scripts as authoritative for those specific numbers until that
  cross-validation is done.

---

## Resources

Labels declare what a process wants; `--max_cpus`, `--max_memory`, `--max_time`
declare what the machine has. Requests are capped to the ceiling, so the same
pipeline runs on a laptop and a cluster without edits.

```bash
nextflow run main.nf -profile docker --max_cpus 40 --max_memory 256.GB ...
```

Retries are deliberately narrow: only out-of-memory and wall-clock kills
(exit 137, 140, 143 and friends) retry, with more resources on each attempt.
Any other failure stops the run, because retrying a genuine bug wastes hours
and hides the cause.

`SNIPPY_CORE` and `DISTANCE_MATRIX` carry the `process_high_memory` label. At
cohort scale these are the steps that will exhaust a small machine first — the
distance matrix is O(n²) and was 137 GB of text for 3,255 genomes.

---

## Stage 6

Stage 6 (external corroboration) is **not** automated. It is a literature and
homology judgement, and a pipeline that claimed to automate it would be
asserting evidence it had not gathered. Supply findings via `--stage6` as JSON:

```json
{ "htz92_2925_m21l": true, "aminotransferase": null }
```

`true` = support found, `false` = searched and none found, `null` = not
assessed. Absent the file, Stage 6 is recorded as missing — never as passed.

`bin/interactive_stage6.py` is a terminal walkthrough that builds this JSON
one candidate at a time (y/n/unknown, resumable via `--resume`). It refuses
to run outside an interactive terminal and is deliberately not wired into
any `.nf` module — run it by hand between a validation run and a rerun with
`--stage6`:

```bash
python3 bin/interactive_stage6.py results/clade/clade_evidence_table.tsv --output stage6.json
nextflow run main.nf -profile docker ... --stage6 stage6.json -resume
```

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | At least one candidate convergent across all required stages |
| 1 | Ran cleanly; nothing convergent |
| 2 | Insufficient evidence to reach a verdict |
| 3 | Error |

Exit 1 is a scientific result, not a failure. Most candidates should not survive.

---

## Troubleshooting

**"Process requirement exceeds available CPUs"** — lower `--max_cpus` to what
the machine actually has.

**AMRFinderPlus database errors** — the database is fetched once per run by
`AMRFINDER_UPDATE` and cached under `--outdir/databases`. Use
`--skip_amrfinder_update` only if the container already carries a database.

**A task failed** — the work directory in the error message contains
`.command.sh` (what ran), `.command.err` and `.command.log` (why it failed).
Fix, then re-run with `-resume`.

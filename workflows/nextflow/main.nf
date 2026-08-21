#!/usr/bin/env nextflow
/*
 * CLADE end-to-end pipeline: raw reads -> six-stage evidence table.
 *
 * Wraps existing, already-verified tools (real commands reconstructed from
 * pipeline_orchestrator_linux_v2.py, see CLADE_Complete_Manuscript.md §5.2)
 * and the already-tested `clade` CLI (Stages 1-5). Stage 6 (external
 * corroboration) is not run here by design — see clade.validation.stage6_corroboration.
 */
nextflow.enable.dsl = 2

include { PER_SAMPLE } from './subworkflows/per_sample'
include { AGGREGATE_AND_VALIDATE } from './subworkflows/aggregate_and_validate'
include { DOWNLOAD } from './modules/download'

workflow {
    if (!params.samplesheet) { error "Provide --samplesheet (CSV: sample_id,accession OR sample_id,fastq_1,fastq_2)" }
    if (!params.reference)   { error "Provide --reference (reference FASTA)" }
    if (!params.candidates)  { error "Provide --candidates (CSV: name,position)" }
    if (!params.resistance_gene) { error "Provide --resistance_gene (AMRFinderPlus gene symbol defining the phenotype)" }

    reference = file(params.reference)
    candidates_file = file(params.candidates)

    samplesheet_ch = Channel.fromPath(params.samplesheet)
        .splitCsv(header: true)

    // Branch: rows with an accession get downloaded; rows with local FASTQ
    // paths are used directly. Matches the samplesheet flexibility documented
    // in the plan (params.samplesheet comment).
    to_download = samplesheet_ch.filter { it.accession }.map { it.accession }
    local_reads = samplesheet_ch.filter { it.fastq_1 }
        .map { tuple(it.sample_id, file(it.fastq_1), file(it.fastq_2)) }

    DOWNLOAD(to_download)
    reads_ch = DOWNLOAD.out.mix(local_reads)

    PER_SAMPLE(reads_ch, reference)

    AGGREGATE_AND_VALIDATE(
        PER_SAMPLE.out.snippy_dirs,
        PER_SAMPLE.out.amrfinder.map { id, f -> f },
        PER_SAMPLE.out.mlst.map { id, f -> f },
        reference,
        candidates_file,
        params.resistance_gene
    )

    AGGREGATE_AND_VALIDATE.out.evidence_table.view { "CLADE evidence table: $it" }
}

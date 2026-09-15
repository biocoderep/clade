#!/usr/bin/env nextflow
/*
========================================================================================
    biocoderep/clade
========================================================================================
    Raw sequencing reads -> six-stage CLADE evidence table.

    CLADE validates candidate compensatory mutations in clonally-reproducing
    bacterial populations, where a variant that merely co-inherited with a
    successful resistant lineage is statistically indistinguishable from one
    that actually compensates. A single structure-corrected association test
    cannot separate them; CLADE requires a candidate to survive six sequential,
    independently-falsifiable checks.

    Documentation: https://github.com/biocoderep/clade
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { PER_SAMPLE             } from './subworkflows/per_sample'
include { AGGREGATE_AND_VALIDATE } from './subworkflows/aggregate_and_validate'
include { DOWNLOAD               } from './modules/download'
include { MULTIQC                } from './modules/multiqc'
include { COLLATE_VERSIONS       } from './modules/collate_versions'

/*
----------------------------------------------------------------------------------------
    HELP / VERSION
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    ============================================================================
     ${workflow.manifest.name} v${workflow.manifest.version}
     ${workflow.manifest.description}
    ============================================================================

    USAGE

      nextflow run biocoderep/clade -profile docker \\
          --samplesheet samples.csv \\
          --reference reference.fa \\
          --candidates candidates.csv \\
          --resistance_gene blaOXA-23 \\
          --outdir results

    REQUIRED

      --samplesheet       CSV. Either  sample_id,accession
                          or          sample_id,fastq_1,fastq_2
                          Mixed rows are allowed.
      --reference         Reference genome FASTA (the coordinate system your
                          candidate positions refer to).
      --candidates        CSV: name,position  -- 1-based positions on the
                          reference. CLADE validates this list; it does not
                          generate it.
      --resistance_gene   AMRFinderPlus gene symbol defining the binary
                          phenotype, e.g. blaOXA-23.

    CLADE STAGE OPTIONS

      --alpha                  Significance threshold                 [${params.alpha}]
      --min_carriers           Minimum carriers to test a candidate   [${params.min_carriers}]
      --max_freq               Above this cohort frequency a candidate
                               is near-fixed and uninformative        [${params.max_freq}]
      --ambiguity_resolution   Fitch tie-breaking: 0 DELTRAN-like,
                               1 ACCTRAN-like                         [${params.ambiguity_resolution}]
      --unique_neighbors       Enforce 1:1 matching in Stage 5        [${params.unique_neighbors}]
      --stage1_results         Use externally computed Stage 1 results
                               (e.g. pyseer output) instead of the
                               built-in Firth fit                     [${params.stage1_results}]
      --stage6                 JSON of manual corroboration findings  [${params.stage6}]

    OUTPUT

      --outdir                 Results directory                      [${params.outdir}]
      --publish_dir_mode       copy | symlink | link                  [${params.publish_dir_mode}]
      --save_trimmed           Keep trimmed FASTQs                    [${params.save_trimmed}]
      --save_assemblies        Keep assemblies                        [${params.save_assemblies}]
      --save_snippy_dirs       Keep per-sample Snippy dirs (large)    [${params.save_snippy_dirs}]

    RESOURCES  (ceilings for THIS machine; requests are capped to them)

      --max_cpus               [${params.max_cpus}]
      --max_memory             [${params.max_memory}]
      --max_time               [${params.max_time}]

    PROFILES

      docker, singularity, apptainer, podman   containerised (recommended)
      conda, mamba                             conda fallback
      test                                     tiny inputs, for -stub-run
      local, slurm                             executor presets

    A five-second wiring check, no tools or data required:

      nextflow run main.nf -profile test,docker -stub-run

    """.stripIndent()
}

/*
----------------------------------------------------------------------------------------
    INPUT VALIDATION
----------------------------------------------------------------------------------------
*/

def validateParams() {
    def errors = []

    if (!params.samplesheet)     errors << "--samplesheet is required"
    if (!params.reference)       errors << "--reference is required"
    if (!params.candidates)      errors << "--candidates is required"
    if (!params.resistance_gene) errors << "--resistance_gene is required"

    // Fail on a missing input file here, at second zero, rather than after the
    // cohort has already been downloaded and assembled.
    [ 'samplesheet': params.samplesheet,
      'reference'  : params.reference,
      'candidates' : params.candidates ].each { name, value ->
        if (value && !file(value).exists()) errors << "--${name} does not exist: ${value}"
    }

    if (params.alpha <= 0 || params.alpha >= 1)  errors << "--alpha must be between 0 and 1 (got ${params.alpha})"
    if (params.max_freq <= 0 || params.max_freq > 1) errors << "--max_freq must be in (0,1] (got ${params.max_freq})"
    if (params.min_carriers < 1)                 errors << "--min_carriers must be >= 1 (got ${params.min_carriers})"
    if (!(params.ambiguity_resolution in [0, 1])) errors << "--ambiguity_resolution must be 0 or 1 (got ${params.ambiguity_resolution})"
    if (!(params.publish_dir_mode in ['copy','symlink','link','rellink','move'])) {
        errors << "--publish_dir_mode '${params.publish_dir_mode}' is not valid"
    }

    if (errors) {
        log.error "Parameter validation failed:\n  - " + errors.join("\n  - ") +
                  "\n\nRun with --help for usage."
        System.exit(1)
    }
}

def summaryLog() {
    log.info """
    ----------------------------------------------------------------------------
     ${workflow.manifest.name} v${workflow.manifest.version}
    ----------------------------------------------------------------------------
     samplesheet      : ${params.samplesheet}
     reference        : ${params.reference}
     candidates       : ${params.candidates}
     resistance gene  : ${params.resistance_gene}
     outdir           : ${params.outdir}
     profile          : ${workflow.profile}
     container engine : ${workflow.containerEngine ?: 'none (conda or local)'}
     max cpus/mem/time: ${params.max_cpus} / ${params.max_memory} / ${params.max_time}
     Nextflow         : ${workflow.nextflow.version}
     run name         : ${workflow.runName}
    ----------------------------------------------------------------------------
    """.stripIndent()
}

/*
----------------------------------------------------------------------------------------
    MAIN WORKFLOW
----------------------------------------------------------------------------------------
*/

workflow {

    if (params.help)    { helpMessage(); System.exit(0) }
    if (params.version) { log.info "${workflow.manifest.name} v${workflow.manifest.version}"; System.exit(0) }

    validateParams()
    summaryLog()

    ch_reference  = Channel.value(file(params.reference,  checkIfExists: true))
    ch_candidates = Channel.value(file(params.candidates, checkIfExists: true))
    ch_versions   = Channel.empty()

    /*
     * A samplesheet row is either an SRA accession to fetch or a pair of local
     * FASTQs. Both are supported in the same sheet: a cohort assembled from
     * public data plus unreleased local isolates is the normal case, not an
     * edge case.
     */
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample_id?.trim() && !row.accession?.trim()) {
                error "Samplesheet row has neither sample_id nor accession: ${row}"
            }
            row
        }

    ch_to_download = ch_samplesheet.filter { it.accession?.trim() }.map { it.accession.trim() }
    ch_local_reads = ch_samplesheet
        .filter { it.fastq_1?.trim() }
        .map { row -> tuple(row.sample_id.trim(),
                            file(row.fastq_1.trim(), checkIfExists: true),
                            file(row.fastq_2.trim(), checkIfExists: true)) }

    DOWNLOAD(ch_to_download)
    ch_versions = ch_versions.mix(DOWNLOAD.out.versions.first().ifEmpty([]))

    ch_reads = DOWNLOAD.out.reads.mix(ch_local_reads)

    PER_SAMPLE(ch_reads, ch_reference)
    ch_versions = ch_versions.mix(PER_SAMPLE.out.versions)

    AGGREGATE_AND_VALIDATE(
        PER_SAMPLE.out.snippy_dirs,
        PER_SAMPLE.out.amrfinder,
        PER_SAMPLE.out.mlst,
        ch_reference,
        ch_candidates,
        params.resistance_gene
    )
    ch_versions = ch_versions.mix(AGGREGATE_AND_VALIDATE.out.versions)

    if (!params.skip_multiqc) {
        MULTIQC(PER_SAMPLE.out.qc.collect().ifEmpty([]))
        ch_versions = ch_versions.mix(MULTIQC.out.versions)
    }

    COLLATE_VERSIONS(ch_versions.unique().collectFile(name: 'collated_versions.yml'))
}

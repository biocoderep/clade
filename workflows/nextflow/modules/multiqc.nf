// Aggregate QC across the cohort.
//
// Worth its own step because the project it comes from lost every FastQC
// report to an aggressive per-sample cleanup, leaving read quality
// un-auditable after the fact. A single persisted MultiQC report prevents
// that failure mode recurring.
process MULTIQC {
    label 'process_single'

    conda "bioconda::multiqc=1.35"
    container "quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1"

    input:
    // stageAs '?/*' puts every input in its own numbered subdirectory.
    // Without it, identically-named reports collide: FastQC emits
    // <sample>_1_fastqc.zip for both the raw and the trimmed pass, and
    // Nextflow refuses to stage two files to one name.
    path multiqc_files, stageAs: '?/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data",        emit: data
    path "versions.yml",        emit: versions

    script:
    """
    multiqc --force .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed 's/multiqc, version //' )
    END_VERSIONS
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_data
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: 1.35
    END_VERSIONS
    """
}

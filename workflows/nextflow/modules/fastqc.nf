// Pre- and post-trim quality control — same module used twice via aliasing
// (FASTQC as FASTQC_RAW / FASTQC as FASTQC_TRIM), matching the orchestrator's
// two identical FastQC calls before and after fastp.
process FASTQC {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    """
    fastqc -o . -t 2 ${r1} ${r2}
    """

    stub:
    """
    touch ${r1.baseName}_fastqc.html ${r1.baseName}_fastqc.zip
    touch ${r2.baseName}_fastqc.html ${r2.baseName}_fastqc.zip
    """
}

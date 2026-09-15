// Read quality control. Included twice by alias (raw and post-trim) so the
// effect of trimming is visible rather than assumed.
process FASTQC {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::fastqc=0.12.1"
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("*.html"), emit: html
    path "*.zip",                         emit: zip
    path "versions.yml",                  emit: versions

    script:
    """
    fastqc --threads ${task.cpus} --outdir . ${r1} ${r2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed 's/FastQC v//' )
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}_1_fastqc.html ${sample_id}_1_fastqc.zip
    touch ${sample_id}_2_fastqc.html ${sample_id}_2_fastqc.zip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: 0.12.1
    END_VERSIONS
    """
}

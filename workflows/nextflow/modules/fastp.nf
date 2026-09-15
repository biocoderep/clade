// Adapter and quality trimming. The JSON report is retained because MultiQC
// needs it, and because per-sample retained-read counts are the first thing to
// check when a downstream sample behaves oddly.
process FASTP {
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "quay.io/biocontainers/fastp:0.23.4--h125f33a_5"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_1.trim.fastq.gz"), path("${sample_id}_2.trim.fastq.gz"), emit: reads
    path "${sample_id}.fastp.json", emit: json
    path "${sample_id}.fastp.html", emit: html
    path "versions.yml",            emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    fastp \\
        --in1 ${r1} --in2 ${r2} \\
        --out1 ${sample_id}_1.trim.fastq.gz \\
        --out2 ${sample_id}_2.trim.fastq.gz \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        --thread ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$( fastp --version 2>&1 | sed 's/fastp //' )
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > ${sample_id}_1.trim.fastq.gz
    echo "" | gzip > ${sample_id}_2.trim.fastq.gz
    echo '{}' > ${sample_id}.fastp.json
    touch ${sample_id}.fastp.html
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: 0.23.4
    END_VERSIONS
    """
}

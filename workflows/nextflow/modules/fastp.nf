// Adapter/quality trimming — real command, default parameters, unchanged
// from pipeline_orchestrator_linux_v2.py step 2.
process FASTP {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_fastp.html"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_trim_1.fastq.gz"), path("${sample_id}_trim_2.fastq.gz"), emit: reads
    path "${sample_id}_fastp.html", emit: report

    script:
    """
    fastp -i ${r1} -I ${r2} -o ${sample_id}_trim_1.fastq.gz -O ${sample_id}_trim_2.fastq.gz -h ${sample_id}_fastp.html --thread ${task.cpus}
    """

    stub:
    """
    touch ${sample_id}_trim_1.fastq.gz ${sample_id}_trim_2.fastq.gz ${sample_id}_fastp.html
    """
}

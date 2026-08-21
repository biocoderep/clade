// De novo assembly — real command, default parameters, unchanged from
// pipeline_orchestrator_linux_v2.py step 5.
process SHOVILL {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_assembly/contigs.fa")

    script:
    """
    shovill --outdir ${sample_id}_assembly --R1 ${r1} --R2 ${r2} --cpus ${task.cpus} --force
    """

    stub:
    """
    mkdir -p ${sample_id}_assembly
    touch ${sample_id}_assembly/contigs.fa
    """
}

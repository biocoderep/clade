// MLST typing on assembly contigs — real command, default scheme
// auto-detection, unchanged from pipeline_orchestrator_linux_v2.py step 6.
process MLST {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_mlst.tsv")

    script:
    """
    mlst ${contigs} > ${sample_id}_mlst.tsv
    """

    stub:
    """
    printf "${contigs}\\tabaumannii_2\\t1\\n" > ${sample_id}_mlst.tsv
    """
}

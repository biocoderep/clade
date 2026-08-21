// AMR gene profiling (ResFinder + NCBI databases) on assembly contigs —
// real commands, default parameters, unchanged from
// pipeline_orchestrator_linux_v2.py step 6.
process ABRICATE {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}_resfinder.tsv"), path("${sample_id}_ncbi.tsv")

    script:
    """
    abricate --db resfinder ${contigs} > ${sample_id}_resfinder.tsv
    abricate --db ncbi ${contigs} > ${sample_id}_ncbi.tsv
    """

    stub:
    """
    touch ${sample_id}_resfinder.tsv ${sample_id}_ncbi.tsv
    """
}

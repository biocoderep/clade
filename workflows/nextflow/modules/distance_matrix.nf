// NEW module (caught during plan recheck): without this, CLADE's --distances
// input is never populated and Stage 5 (matched-neighbor comparison) silently
// never runs. See bin/distance_matrix.py.
process DISTANCE_MATRIX {
    conda "${projectDir}/conda/genome_processing.yml"
    publishDir "${params.outdir}/core_genome", mode: 'copy'

    input:
    path tree

    output:
    path "distances.tsv"

    script:
    """
    python3 ${projectDir}/bin/distance_matrix.py --tree ${tree} --output distances.tsv
    """

    stub:
    """
    touch distances.tsv
    """
}

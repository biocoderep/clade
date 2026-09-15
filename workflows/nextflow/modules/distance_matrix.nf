// Pairwise patristic distances from the phylogeny, for CLADE Stage 5's
// matched-nearest-neighbour comparison.
//
// Separated from the Stage 5 step because this matrix is O(n^2) and expensive
// to recompute: at cohort scale it is the largest intermediate the pipeline
// produces, and Stage 5 may reasonably be re-run with different matching
// options against the same matrix.
process DISTANCE_MATRIX {
    label 'process_high_memory'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path tree

    output:
    path "distances.tsv", emit: distances
    path "versions.yml",  emit: versions

    script:
    """
    distance_matrix.py --tree ${tree} --output distances.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: \$( python3 -c "import dendropy; print(dendropy.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    printf '\\ts1\\ns1\\t0.0\\n' > distances.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: 5.0.11
    END_VERSIONS
    """
}

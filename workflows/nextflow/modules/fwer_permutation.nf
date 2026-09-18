// Family-wise error across the whole six-stage screen, by lineage-preserving
// permutation (see bin/fwer_permutation.py). Answers what the framework's
// stage design asserts but does not by itself measure: under no true
// candidate-phenotype relationship, how often does ANY candidate in the
// family reach convergence.
//
// Expensive -- each permutation re-runs Stages 1, 3, 4 and 5 for every
// candidate that survives the stage before it -- so it is opt-in
// (params.run_fwer_permutation) and not part of the default run.
process FWER_PERMUTATION {
    label 'process_high'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path genotypes
    path phenotype
    path lineages
    path tree
    path distances

    output:
    path "fwer_permutation.txt", emit: summary
    path "fwer_permutation.tsv", emit: raw
    path "versions.yml",         emit: versions

    script:
    """
    fwer_permutation.py \\
        --genotypes ${genotypes} \\
        --phenotype ${phenotype} \\
        --lineages ${lineages} \\
        --tree ${tree} \\
        --distances ${distances} \\
        --n-permutations ${params.fwer_n_permutations} \\
        --seed ${params.fwer_seed} \\
        --outdir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: \$( python3 -c "import dendropy; print(dendropy.__version__)" )
        lifelines: \$( python3 -c "import lifelines; print(lifelines.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    echo "FWER stub" > fwer_permutation.txt
    printf 'permutation\\tn_convergent\\n0\\t0\\n' > fwer_permutation.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dendropy: 5.0.11
    END_VERSIONS
    """
}

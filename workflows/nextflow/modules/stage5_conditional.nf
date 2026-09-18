// Stage 5 by conditional logistic regression, stratified on the control
// genome (see bin/stage5_conditional.py for why this replaced plain McNemar).
//
// Optional, and off by default (see params.run_stage5_conditional): it needs
// `lifelines`, and it is slower than the McNemar path already inside
// RUN_CLADE, since it fits one stratified Cox model per candidate rather than
// one contingency table.
process STAGE5_CONDITIONAL {
    label 'process_high'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path genotypes
    path phenotype
    path distances
    val  candidate_names

    output:
    path "stage5_conditional.tsv", emit: results
    path "versions.yml",           emit: versions

    script:
    def cand_args = candidate_names.join(' ')
    """
    stage5_conditional.py \\
        --genotypes ${genotypes} \\
        --phenotype ${phenotype} \\
        --distances ${distances} \\
        --candidates ${cand_args} \\
        --output stage5_conditional.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lifelines: \$( python3 -c "import lifelines; print(lifelines.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    echo -e "Candidate\\tconditional_p\\tconditional_coef" > stage5_conditional.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lifelines: 0.30.3
    END_VERSIONS
    """
}

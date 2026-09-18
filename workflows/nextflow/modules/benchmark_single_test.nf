// Compare CLADE's dispositions against a conventional single-test pipeline
// (classical MDS of the distance matrix + logistic regression + BH-FDR --
// pyseer's fixed-effects distance model, reimplemented; see
// bin/benchmark_single_test.py for why pyseer itself is not vendored).
//
// Optional (params.run_benchmark): answers the reviewer question a framework
// paper is judged on -- what does the extra machinery buy over the standard
// approach -- but is not needed for CLADE's own evidence table.
process BENCHMARK_SINGLE_TEST {
    label 'process_medium'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path genotypes
    path phenotype
    path distances
    path clade_evidence_tsv

    output:
    path "benchmark_single_test.tsv", emit: results
    path "versions.yml",              emit: versions

    script:
    """
    benchmark_single_test.py \\
        --genotypes ${genotypes} \\
        --phenotype ${phenotype} \\
        --distances ${distances} \\
        --clade-results ${clade_evidence_tsv} \\
        --max-dimensions ${params.benchmark_mds_dimensions} \\
        --output benchmark_single_test.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        statsmodels: \$( python3 -c "import statsmodels; print(statsmodels.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    echo -e "Candidate\\tsingletest_p\\tsingletest_hit" > benchmark_single_test.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        statsmodels: 0.14.0
    END_VERSIONS
    """
}

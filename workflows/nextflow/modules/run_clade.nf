// The six-stage validation itself.
//
// Stage 6 (external corroboration) is not automated here by design: it is a
// literature and homology judgement, and a pipeline that pretended to
// automate it would be asserting evidence it had not gathered. Supply findings
// via --stage6 as JSON if you have them; absent that, Stage 6 is recorded as
// missing, never as passed.
//
// Exit codes are meaningful and deliberately not swallowed:
//   0 = at least one convergent candidate
//   1 = ran cleanly, nothing convergent
//   2 = insufficient evidence for a verdict
//   3 = error
process RUN_CLADE {
    label 'process_medium'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path genotypes
    path phenotype
    path lineages
    path candidates_file
    path tree
    path distances
    path stage1_results
    path stage6

    output:
    path "clade_evidence_table.md",  emit: evidence_table
    path "clade_evidence_table.tsv", emit: evidence_tsv
    path "versions.yml",             emit: versions

    script:
    def stage1_arg = stage1_results.name != 'NO_STAGE1' ? "--stage1-results ${stage1_results}" : ''
    def stage6_arg = stage6.name         != 'NO_STAGE6' ? "--stage6 ${stage6}"                 : ''
    def unique_arg = params.unique_neighbors ? '--unique-neighbors' : ''
    """
    run_clade.py \\
        --genotypes ${genotypes} \\
        --phenotype ${phenotype} \\
        --lineages ${lineages} \\
        --candidates-file ${candidates_file} \\
        --tree ${tree} \\
        --distances ${distances} \\
        --output clade_evidence_table.md \\
        --output-tsv clade_evidence_table.tsv \\
        --alpha ${params.alpha} \\
        --min-carriers ${params.min_carriers} \\
        --max-freq ${params.max_freq} \\
        --temporal-margin ${params.temporal_margin} \\
        --ambiguity-resolution ${params.ambiguity_resolution} \\
        ${unique_arg} ${stage1_arg} ${stage6_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: \$( python3 -c "import clade; print(clade.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    echo "# CLADE evidence table (stub)" > clade_evidence_table.md
    printf 'Candidate\\tDisposition\\ncand_1\\tinsufficient_evidence\\n' > clade_evidence_table.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """
}

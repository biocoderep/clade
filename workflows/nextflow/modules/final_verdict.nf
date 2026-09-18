// The pipeline's last word: a candidate, or no candidate, stated once. See
// bin/final_verdict.py for the exact rule -- it mirrors CLADE's own
// disposition logic and loosens nothing.
process FINAL_VERDICT {
    label 'process_single'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path evidence_tsv
    path stage5_conditional  // may be the NO_STAGE5_CONDITIONAL placeholder

    output:
    path "FINAL_VERDICT.md", emit: verdict
    path "versions.yml",     emit: versions

    script:
    def cond_arg = stage5_conditional.name != 'NO_STAGE5_CONDITIONAL' \
        ? "--stage5-conditional ${stage5_conditional}" : ''
    """
    final_verdict.py \\
        --evidence-tsv ${evidence_tsv} \\
        ${cond_arg} \\
        --output FINAL_VERDICT.md

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """

    stub:
    """
    echo "# CLADE Final Verdict (stub)" > FINAL_VERDICT.md
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """
}

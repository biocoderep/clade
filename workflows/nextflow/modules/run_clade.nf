// Final step — calls the existing, tested `clade` CLI (Stages 1-5). Stage 6
// (external corroboration) is deliberately not run here; see the plan's
// "Stage coverage boundary" note and clade.validation.stage6_corroboration.
process RUN_CLADE {
    conda "${projectDir}/../../environment.yml"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path genotypes
    path phenotype
    path lineages
    path candidates_file
    path tree
    path distances

    output:
    path "evidence_table.md"

    script:
    """
    python3 ${projectDir}/bin/run_clade.py \
        --genotypes ${genotypes} \
        --phenotype ${phenotype} \
        --lineages ${lineages} \
        --candidates-file ${candidates_file} \
        --tree ${tree} \
        --distances ${distances} \
        --output evidence_table.md
    """

    stub:
    """
    touch evidence_table.md
    """
}

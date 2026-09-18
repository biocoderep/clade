// Generate a genome-wide, unranked candidates.csv directly from core.vcf
// (see bin/discover_candidates.py for why this makes no ranking claim, unlike
// the invalid MI/DCA discovery step this project's own case study found and
// discarded). Only runs when --candidates is not supplied.
process DISCOVER_CANDIDATES {
    label 'process_medium'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path core_vcf

    output:
    path "candidates.csv", emit: candidates
    path "versions.yml",   emit: versions

    script:
    def cap_arg = params.auto_candidates_max_sites ? "--max-sites ${params.auto_candidates_max_sites}" : ''
    """
    discover_candidates.py \\
        --vcf ${core_vcf} \\
        --output candidates.csv \\
        --max-missing-fraction ${params.auto_candidates_max_missing} \\
        ${cap_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """

    stub:
    """
    echo "name,position" > candidates.csv
    echo "stub_1,chr:1" >> candidates.csv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """
}

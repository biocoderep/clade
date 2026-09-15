// Assembly-based acquired-AMR gene screen. Reported alongside AMRFinderPlus
// rather than instead of it: the two use different databases and disagree at
// the margins, and that disagreement is informative.
process ABRICATE {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::abricate=1.4.0"
    container "quay.io/biocontainers/abricate:1.4.0--h05cac1d_0"

    input:
    tuple val(sample_id), path(contigs)
    val database

    output:
    tuple val(sample_id), path("${sample_id}.abricate.tsv"), emit: report
    path "versions.yml",                                     emit: versions

    script:
    """
    abricate --db ${database} --threads ${task.cpus} ${contigs} > ${sample_id}.abricate.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: \$( abricate --version 2>&1 | sed 's/abricate //' )
    END_VERSIONS
    """

    stub:
    """
    echo -e "#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE" > ${sample_id}.abricate.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abricate: 1.4.0
    END_VERSIONS
    """
}

// Per-sample AMR gene detection. The phenotype used by CLADE Stage 1 is
// derived from this output (params.resistance_gene), so it is the single most
// consequential step for the association result.
process AMRFINDER {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::ncbi-amrfinderplus=4.0.23"
    container "quay.io/biocontainers/ncbi-amrfinderplus:4.0.23--hf69ffd2_0"

    input:
    tuple val(sample_id), path(contigs)
    path db

    output:
    tuple val(sample_id), path("${sample_id}.amrfinder.tsv"), emit: report
    path "versions.yml",                                      emit: versions

    script:
    def db_arg = db.name != 'NO_DB' ? "--database ${db}" : ''
    """
    amrfinder --nucleotide ${contigs} ${db_arg} --threads ${task.cpus} > ${sample_id}.amrfinder.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$( amrfinder --version )
    END_VERSIONS
    """

    stub:
    """
    echo -e "Protein identifier\\tGene symbol\\tSequence name" > ${sample_id}.amrfinder.tsv
    echo -e "-\\t${params.resistance_gene}\\tstub" >> ${sample_id}.amrfinder.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: 4.0.23
    END_VERSIONS
    """
}

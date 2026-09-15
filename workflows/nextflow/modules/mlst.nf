// Sequence typing. MLST assignments are the lineage covariate for CLADE
// Stage 1 and the stratum for Stage 2's recurrence check, so this is load-
// bearing for the population-structure correction, not just descriptive.
process MLST {
    tag "$sample_id"
    label 'process_low'

    conda "bioconda::mlst=2.35.0"
    container "quay.io/biocontainers/mlst:2.35.0--hdfd78af_0"

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.mlst.tsv"), emit: report
    path "versions.yml",                                 emit: versions

    script:
    """
    mlst --threads ${task.cpus} ${contigs} > ${sample_id}.mlst.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: \$( mlst --version 2>&1 | sed 's/mlst //' )
    END_VERSIONS
    """

    stub:
    """
    echo -e "${sample_id}.fa\\tabaumannii\\t2" > ${sample_id}.mlst.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mlst: 2.35.0
    END_VERSIONS
    """
}

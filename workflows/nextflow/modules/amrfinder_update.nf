// One-time AMRFinderPlus database fetch.
//
// AMRFinderPlus hard-errors on a fresh install with no database. Running this
// once and gating every AMRFINDER task on its output makes that a real dataflow
// dependency rather than a documented prerequisite a user can miss — and means
// the database is fetched once per run, not once per sample.
process AMRFINDER_UPDATE {
    label 'process_network'
    storeDir "${params.outdir}/databases"

    conda "bioconda::ncbi-amrfinderplus=4.0.23"
    container "quay.io/biocontainers/ncbi-amrfinderplus:4.0.23--hf69ffd2_0"

    output:
    path "amrfinder_db", emit: db
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p amrfinder_db
    amrfinder_update --database amrfinder_db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$( amrfinder --version )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p amrfinder_db
    touch amrfinder_db/.stub
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: 4.0.23
    END_VERSIONS
    """
}

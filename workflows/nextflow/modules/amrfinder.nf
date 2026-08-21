// AMRFinderPlus profiling on assembly contigs, nucleotide mode — real command,
// unchanged from pipeline_orchestrator_linux_v2.py step 6. This is the same
// tool/mode CLADE's own phenotype provenance derivation
// (clade.provenance.phenotype.derive_phenotype_from_amr_matrix) expects its
// input Rtab to come from.
//
// `db_ready` is not used by the command itself — it exists purely to create
// a real Nextflow dataflow dependency on AMRFINDER_UPDATE, so every sample's
// AMRFINDER call waits for the one-time database download to finish first
// (see modules/amrfinder_update.nf for why this is required, not optional).
process AMRFINDER {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(contigs)
    val db_ready

    output:
    tuple val(sample_id), path("${sample_id}_amrfinder.tsv")

    script:
    """
    amrfinder -n ${contigs} > ${sample_id}_amrfinder.tsv
    """

    stub:
    """
    touch ${sample_id}_amrfinder.tsv
    """
}

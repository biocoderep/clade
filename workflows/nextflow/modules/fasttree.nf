// Phylogeny — real command, unchanged from the confirmed, byte-reproduced
// provenance of tree.nwk this session (FastTree -nt -gtr on the SNP-only
// alignment). Deliberately hard-pinned to core.aln, NOT core.full.aln:
// two independent full-cohort attempts against core.full.aln were killed by
// the OS out-of-memory handler after 24-27 hours each (peak 497-521GB) — see
// docs/limitations/phylogenetic_alignment_sensitivity.md. Switching inputs
// is an explicit, documented opt-in for users who have the infrastructure
// for it, never the default.
process FASTTREE {
    conda "${projectDir}/conda/genome_processing.yml"
    publishDir "${params.outdir}/core_genome", mode: 'copy'

    input:
    path core_aln

    output:
    path "tree.nwk"

    script:
    """
    FastTree -nt -gtr ${core_aln} > tree.nwk
    """

    stub:
    """
    echo '(A:0.1,B:0.1,(C:0.1,D:0.1):0.1);' > tree.nwk
    """
}

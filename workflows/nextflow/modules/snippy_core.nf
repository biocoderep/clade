// Core-genome alignment — real command, default masking (0 regions masked),
// unchanged from the confirmed real snippy-core 4.6.0 invocation
// (see CLADE_Complete_Manuscript.md §5.4 / snippy_core_full_run_v2.log).
// Produces both core.aln (SNP-only) and core.full.aln (every core-genome site) —
// only core.aln is used downstream by default; see FASTTREE for why.
process SNIPPY_CORE {
    conda "${projectDir}/conda/genome_processing.yml"
    publishDir "${params.outdir}/core_genome", mode: 'copy'

    input:
    path snippy_dirs
    path reference

    output:
    path "core.aln", emit: core_aln
    path "core.full.aln", emit: core_full_aln
    path "core.vcf", emit: core_vcf

    script:
    """
    snippy-core --ref ${reference} ${snippy_dirs}
    """

    stub:
    """
    touch core.aln core.full.aln core.vcf core.tab core.txt
    """
}

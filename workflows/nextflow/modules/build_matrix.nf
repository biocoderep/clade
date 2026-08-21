// Bridges per-sample profiling tables + core.vcf to CLADE's expected input
// format (genotype matrix, phenotype vector, lineage file). See
// bin/build_matrix.py for why candidate genotypes come from core.vcf, not
// the profiling tables, and why the phenotype is derived (not hand-written).
//
// amrfinder_files / mlst_files are collected lists — Nextflow stages every
// file in a `path` input list into the process's working directory, so
// build_matrix.py's `--profiling-dir .` finds them all via glob without any
// manual directory assembly.
process BUILD_MATRIX {
    conda "${projectDir}/conda/genome_processing.yml"
    publishDir "${params.outdir}/clade_input", mode: 'copy'

    input:
    path core_vcf
    path candidates_file
    path amrfinder_files
    path mlst_files
    val resistance_gene

    output:
    path "out/genotypes.csv", emit: genotypes
    path "out/pheno.tsv", emit: phenotype
    path "out/lineages.tsv", emit: lineages

    script:
    """
    python3 ${projectDir}/bin/build_matrix.py \
        --vcf ${core_vcf} \
        --candidates-file ${candidates_file} \
        --profiling-dir . \
        --resistance-gene "${resistance_gene}" \
        --outdir out
    """

    stub:
    """
    mkdir -p out
    touch out/genotypes.csv out/pheno.tsv out/lineages.tsv
    """
}

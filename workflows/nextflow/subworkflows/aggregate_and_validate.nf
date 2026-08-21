include { SNIPPY_CORE } from '../modules/snippy_core'
include { FASTTREE } from '../modules/fasttree'
include { DISTANCE_MATRIX } from '../modules/distance_matrix'
include { BUILD_MATRIX } from '../modules/build_matrix'
include { RUN_CLADE } from '../modules/run_clade'

workflow AGGREGATE_AND_VALIDATE {
    take:
    snippy_dirs      // collected list, one per sample
    amrfinder_files  // collected list, one per sample
    mlst_files       // collected list, one per sample
    reference
    candidates_file
    resistance_gene

    main:
    SNIPPY_CORE(snippy_dirs.collect(), reference)
    FASTTREE(SNIPPY_CORE.out.core_aln)
    DISTANCE_MATRIX(FASTTREE.out)

    BUILD_MATRIX(
        SNIPPY_CORE.out.core_vcf,
        candidates_file,
        amrfinder_files.collect(),
        mlst_files.collect(),
        resistance_gene
    )

    RUN_CLADE(
        BUILD_MATRIX.out.genotypes,
        BUILD_MATRIX.out.phenotype,
        BUILD_MATRIX.out.lineages,
        candidates_file,
        FASTTREE.out,
        DISTANCE_MATRIX.out
    )

    emit:
    evidence_table = RUN_CLADE.out
}

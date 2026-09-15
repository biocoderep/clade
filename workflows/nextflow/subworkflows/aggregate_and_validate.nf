/*
 * Everything that requires the whole cohort: core-genome alignment, phylogeny,
 * the distance matrix, matrix assembly, and the six-stage validation.
 */

include { SNIPPY_CORE     } from '../modules/snippy_core'
include { FASTTREE        } from '../modules/fasttree'
include { DISTANCE_MATRIX } from '../modules/distance_matrix'
include { BUILD_MATRIX    } from '../modules/build_matrix'
include { RUN_CLADE       } from '../modules/run_clade'

workflow AGGREGATE_AND_VALIDATE {

    take:
    snippy_dirs        // channel: collected per-sample snippy directories
    amrfinder_reports  // channel: collected AMRFinderPlus TSVs
    mlst_reports       // channel: collected MLST TSVs
    reference          // path
    candidates_file    // path
    resistance_gene    // val

    main:
    ch_versions = Channel.empty()

    SNIPPY_CORE(snippy_dirs.collect(), reference)
    ch_versions = ch_versions.mix(SNIPPY_CORE.out.versions)

    // Tree from the SNP-only alignment; core.full.aln is published alongside
    // so the ascertainment-bias sensitivity check remains reproducible.
    FASTTREE(SNIPPY_CORE.out.aln)
    ch_versions = ch_versions.mix(FASTTREE.out.versions)

    DISTANCE_MATRIX(FASTTREE.out.tree)
    ch_versions = ch_versions.mix(DISTANCE_MATRIX.out.versions)

    BUILD_MATRIX(
        SNIPPY_CORE.out.vcf,
        candidates_file,
        amrfinder_reports.mix(mlst_reports).collect(),
        resistance_gene
    )
    ch_versions = ch_versions.mix(BUILD_MATRIX.out.versions)

    // Optional inputs are passed as named placeholder files so the process
    // signature stays fixed; the module checks the name rather than branching
    // the DAG on presence.
    ch_stage1 = params.stage1_results ? Channel.fromPath(params.stage1_results, checkIfExists: true)
                                      : Channel.value(file("${projectDir}/assets/NO_STAGE1"))
    ch_stage6 = params.stage6         ? Channel.fromPath(params.stage6, checkIfExists: true)
                                      : Channel.value(file("${projectDir}/assets/NO_STAGE6"))

    RUN_CLADE(
        BUILD_MATRIX.out.genotypes,
        BUILD_MATRIX.out.phenotype,
        BUILD_MATRIX.out.lineages,
        candidates_file,
        FASTTREE.out.tree,
        DISTANCE_MATRIX.out.distances,
        ch_stage1,
        ch_stage6
    )
    ch_versions = ch_versions.mix(RUN_CLADE.out.versions)

    emit:
    evidence_table = RUN_CLADE.out.evidence_table
    evidence_tsv   = RUN_CLADE.out.evidence_tsv
    tree           = FASTTREE.out.tree
    core_vcf       = SNIPPY_CORE.out.vcf
    versions       = ch_versions
}

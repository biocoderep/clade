/*
 * Everything that requires the whole cohort: core-genome alignment, phylogeny,
 * the distance matrix, matrix assembly, and the six-stage validation.
 */

include { SNIPPY_CORE     } from '../modules/snippy_core'
include { FASTTREE        } from '../modules/fasttree'
include { DISTANCE_MATRIX } from '../modules/distance_matrix'
include { BUILD_MATRIX    } from '../modules/build_matrix'
include { RUN_CLADE             } from '../modules/run_clade'
include { STAGE5_CONDITIONAL    } from '../modules/stage5_conditional'
include { FWER_PERMUTATION      } from '../modules/fwer_permutation'
include { BENCHMARK_SINGLE_TEST } from '../modules/benchmark_single_test'

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

    // --- optional deeper validation, off by default ------------------------
    //
    // None of these three change CLADE's own evidence table. They are
    // additional scrutiny of it: a Stage 5 re-analysis that models matched-
    // pair reuse instead of assuming independence, a permutation estimate of
    // the whole screen's family-wise error, and a comparison against what a
    // conventional single-test pipeline would have returned on the same data.
    // Each re-fits many candidates and is materially slower than RUN_CLADE
    // itself, hence opt-in.
    ch_stage5_conditional = Channel.empty()
    if (params.run_stage5_conditional) {
        ch_candidate_names = candidates_file
            .splitCsv(header: true)
            .map { it.name }
            .collect()
        STAGE5_CONDITIONAL(
            BUILD_MATRIX.out.genotypes,
            BUILD_MATRIX.out.phenotype,
            DISTANCE_MATRIX.out.distances,
            ch_candidate_names
        )
        ch_versions = ch_versions.mix(STAGE5_CONDITIONAL.out.versions)
        ch_stage5_conditional = STAGE5_CONDITIONAL.out.results
    }

    ch_fwer = Channel.empty()
    if (params.run_fwer_permutation) {
        FWER_PERMUTATION(
            BUILD_MATRIX.out.genotypes,
            BUILD_MATRIX.out.phenotype,
            BUILD_MATRIX.out.lineages,
            FASTTREE.out.tree,
            DISTANCE_MATRIX.out.distances
        )
        ch_versions = ch_versions.mix(FWER_PERMUTATION.out.versions)
        ch_fwer = FWER_PERMUTATION.out.summary
    }

    ch_benchmark = Channel.empty()
    if (params.run_benchmark) {
        BENCHMARK_SINGLE_TEST(
            BUILD_MATRIX.out.genotypes,
            BUILD_MATRIX.out.phenotype,
            DISTANCE_MATRIX.out.distances,
            RUN_CLADE.out.evidence_tsv
        )
        ch_versions = ch_versions.mix(BENCHMARK_SINGLE_TEST.out.versions)
        ch_benchmark = BENCHMARK_SINGLE_TEST.out.results
    }

    emit:
    evidence_table     = RUN_CLADE.out.evidence_table
    evidence_tsv       = RUN_CLADE.out.evidence_tsv
    tree               = FASTTREE.out.tree
    core_vcf           = SNIPPY_CORE.out.vcf
    stage5_conditional = ch_stage5_conditional
    fwer_summary       = ch_fwer
    benchmark          = ch_benchmark
    versions           = ch_versions
}

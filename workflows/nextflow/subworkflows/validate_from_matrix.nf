/*
 * Entry point for cohorts that are already past genome processing --
 * genotype matrix, phenotype, lineages, phylogeny and distance matrix
 * already extracted, exactly the state this project's real case study was
 * actually run from (the raw-reads path in AGGREGATE_AND_VALIDATE has never
 * been exercised end to end on real data; this one has, by hand, in the
 * session that produced these files).
 *
 * Runs the identical RUN_CLADE / STAGE5_CONDITIONAL / FWER_PERMUTATION /
 * BENCHMARK_SINGLE_TEST modules as the raw-reads path -- nothing about
 * validation itself is reimplemented here, only genome processing is
 * skipped. A result produced through this entry point and one produced by
 * running AGGREGATE_AND_VALIDATE to completion on the same underlying reads
 * should be identical, module for module.
 */

include { RUN_CLADE             } from '../modules/run_clade'
include { STAGE5_CONDITIONAL    } from '../modules/stage5_conditional'
include { FWER_PERMUTATION      } from '../modules/fwer_permutation'
include { BENCHMARK_SINGLE_TEST } from '../modules/benchmark_single_test'

workflow VALIDATE_FROM_MATRIX {

    take:
    genotypes        // path: Sample_ID-indexed CSV, one 0/1 column per candidate
    phenotype        // path: sample<TAB>phenotype TSV
    lineages         // path: sample<TAB>ST TSV (extra columns ignored)
    candidates_file  // path: name,position CSV (position unused in this mode)
    tree             // path: Newick phylogeny
    distances        // path: sample x sample patristic distance TSV

    main:
    ch_versions = Channel.empty()

    ch_stage1 = params.stage1_results ? Channel.fromPath(params.stage1_results, checkIfExists: true)
                                      : Channel.value(file("${projectDir}/assets/NO_STAGE1"))
    ch_stage6 = params.stage6         ? Channel.fromPath(params.stage6, checkIfExists: true)
                                      : Channel.value(file("${projectDir}/assets/NO_STAGE6"))

    RUN_CLADE(
        genotypes, phenotype, lineages, candidates_file, tree, distances,
        ch_stage1, ch_stage6
    )
    ch_versions = ch_versions.mix(RUN_CLADE.out.versions)

    ch_stage5_conditional = Channel.empty()
    if (params.run_stage5_conditional) {
        ch_candidate_names = candidates_file.splitCsv(header: true).map { it.name }.collect()
        STAGE5_CONDITIONAL(genotypes, phenotype, distances, ch_candidate_names)
        ch_versions = ch_versions.mix(STAGE5_CONDITIONAL.out.versions)
        ch_stage5_conditional = STAGE5_CONDITIONAL.out.results
    }

    ch_fwer = Channel.empty()
    if (params.run_fwer_permutation) {
        FWER_PERMUTATION(genotypes, phenotype, lineages, tree, distances)
        ch_versions = ch_versions.mix(FWER_PERMUTATION.out.versions)
        ch_fwer = FWER_PERMUTATION.out.summary
    }

    ch_benchmark = Channel.empty()
    if (params.run_benchmark) {
        BENCHMARK_SINGLE_TEST(genotypes, phenotype, distances, RUN_CLADE.out.evidence_tsv)
        ch_versions = ch_versions.mix(BENCHMARK_SINGLE_TEST.out.versions)
        ch_benchmark = BENCHMARK_SINGLE_TEST.out.results
    }

    emit:
    evidence_table     = RUN_CLADE.out.evidence_table
    evidence_tsv       = RUN_CLADE.out.evidence_tsv
    stage5_conditional = ch_stage5_conditional
    fwer_summary       = ch_fwer
    benchmark          = ch_benchmark
    versions           = ch_versions
}

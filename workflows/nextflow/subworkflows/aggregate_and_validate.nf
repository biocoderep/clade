/*
 * Everything that requires the whole cohort: core-genome alignment, phylogeny,
 * the distance matrix, matrix assembly, and the six-stage validation.
 *
 * Candidates and the resistance-gene phenotype are each EITHER supplied by
 * the user OR generated automatically -- see DISCOVER_CANDIDATES and
 * DETECT_RESISTANCE_GENE. Both are deferred to inside this subworkflow
 * rather than resolved in main.nf, because automatic candidate discovery
 * needs core.vcf, which does not exist until SNIPPY_CORE has run.
 */

include { SNIPPY_CORE            } from '../modules/snippy_core'
include { FASTTREE               } from '../modules/fasttree'
include { DISTANCE_MATRIX        } from '../modules/distance_matrix'
include { DETECT_RESISTANCE_GENE } from '../modules/detect_resistance_gene'
include { DISCOVER_CANDIDATES    } from '../modules/discover_candidates'
include { BUILD_MATRIX           } from '../modules/build_matrix'
include { RUN_CLADE              } from '../modules/run_clade'
include { STAGE5_CONDITIONAL     } from '../modules/stage5_conditional'
include { FWER_PERMUTATION       } from '../modules/fwer_permutation'
include { BENCHMARK_SINGLE_TEST  } from '../modules/benchmark_single_test'
include { FINAL_VERDICT          } from '../modules/final_verdict'

workflow AGGREGATE_AND_VALIDATE {

    take:
    snippy_dirs         // channel: collected per-sample snippy directories
    amrfinder_reports   // channel: per-sample AMRFinderPlus TSVs, uncollected
    mlst_reports        // channel: collected MLST TSVs
    reference            // path
    candidates_param     // val: path string to a name,position CSV, or null/false
    resistance_gene_param // val: gene symbol string, or null/false

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

    // --- phenotype: supplied, or detected from the cohort's own AMR calls --
    ch_amrfinder_collected = amrfinder_reports.collect()
    if (resistance_gene_param) {
        ch_resistance_gene = Channel.value(resistance_gene_param)
    } else {
        log.info "No --resistance_gene supplied: detecting one automatically " +
                 "from the cohort's own AMRFinderPlus calls (see " +
                 "results/detect_resistance_gene/resistance_gene_candidates.tsv " +
                 "once this run completes -- this is a heuristic starting point, " +
                 "not a validated phenotype assignment)."
        DETECT_RESISTANCE_GENE(ch_amrfinder_collected, ch_amrfinder_collected.map { it.size() })
        ch_versions = ch_versions.mix(DETECT_RESISTANCE_GENE.out.versions)
        // A plain file output, not `env`, to sidestep a static-analysis issue
        // that Nextflow 26's parser hits on `env` outputs in this module.
        ch_resistance_gene = DETECT_RESISTANCE_GENE.out.gene_file.map { it.text.trim() }
    }

    // --- candidates: supplied, or discovered genome-wide from core.vcf -----
    if (candidates_param) {
        ch_candidates = Channel.value(file(candidates_param, checkIfExists: true))
    } else {
        log.info "No --candidates supplied: generating an unranked, genome-wide " +
                 "candidate list directly from core.vcf (see DISCOVER_CANDIDATES " +
                 "-- every core-genome SNP becomes one candidate; CLADE's own six " +
                 "stages, not this step, separate signal from clonal noise)."
        DISCOVER_CANDIDATES(SNIPPY_CORE.out.vcf)
        ch_versions = ch_versions.mix(DISCOVER_CANDIDATES.out.versions)
        ch_candidates = DISCOVER_CANDIDATES.out.candidates
    }

    BUILD_MATRIX(
        SNIPPY_CORE.out.vcf,
        ch_candidates,
        amrfinder_reports.mix(mlst_reports).collect(),
        ch_resistance_gene
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
        ch_candidates,
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
    ch_stage5_conditional = Channel.value(file("${projectDir}/assets/NO_STAGE5_CONDITIONAL"))
    if (params.run_stage5_conditional) {
        ch_candidate_names = ch_candidates
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

    // --- the pipeline's last word: a candidate, or no candidate ------------
    FINAL_VERDICT(RUN_CLADE.out.evidence_tsv, ch_stage5_conditional)
    ch_versions = ch_versions.mix(FINAL_VERDICT.out.versions)

    emit:
    evidence_table     = RUN_CLADE.out.evidence_table
    evidence_tsv       = RUN_CLADE.out.evidence_tsv
    verdict            = FINAL_VERDICT.out.verdict
    tree               = FASTTREE.out.tree
    core_vcf           = SNIPPY_CORE.out.vcf
    stage5_conditional = ch_stage5_conditional
    fwer_summary       = ch_fwer
    benchmark          = ch_benchmark
    versions           = ch_versions
}

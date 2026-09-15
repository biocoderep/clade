/*
 * Everything that happens once per genome: QC, trimming, variant calling,
 * assembly, and assembly-based typing.
 *
 * Two independent paths run from the same trimmed reads. Variants come from
 * the reference-based Snippy path; typing and AMR calls come from the assembly
 * path. They are kept separate so that assembly error cannot reach the
 * genotype matrix that Stage 1 tests.
 */

include { FASTQC as FASTQC_RAW  } from '../modules/fastqc'
include { FASTQC as FASTQC_TRIM } from '../modules/fastqc'
include { FASTP                 } from '../modules/fastp'
include { SNIPPY                } from '../modules/snippy'
include { SHOVILL               } from '../modules/shovill'
include { ABRICATE              } from '../modules/abricate'
include { AMRFINDER_UPDATE      } from '../modules/amrfinder_update'
include { AMRFINDER             } from '../modules/amrfinder'
include { MLST                  } from '../modules/mlst'

workflow PER_SAMPLE {

    take:
    reads_ch      // channel: tuple(sample_id, r1, r2)
    reference     // path

    main:
    ch_versions = Channel.empty()
    ch_qc       = Channel.empty()

    if (!params.skip_fastqc) {
        FASTQC_RAW(reads_ch)
        ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())
        ch_qc       = ch_qc.mix(FASTQC_RAW.out.zip)
    }

    FASTP(reads_ch)
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_qc       = ch_qc.mix(FASTP.out.json)

    if (!params.skip_fastqc) {
        FASTQC_TRIM(FASTP.out.reads)
        ch_versions = ch_versions.mix(FASTQC_TRIM.out.versions.first())
        ch_qc       = ch_qc.mix(FASTQC_TRIM.out.zip)
    }

    // --- reference-based path: the genotypes Stage 1 will test ---
    SNIPPY(FASTP.out.reads, reference)
    ch_versions = ch_versions.mix(SNIPPY.out.versions.first())

    // --- assembly path: typing and AMR only ---
    SHOVILL(FASTP.out.reads)
    ch_versions = ch_versions.mix(SHOVILL.out.versions.first())

    ABRICATE(SHOVILL.out.contigs, params.abricate_db)
    ch_versions = ch_versions.mix(ABRICATE.out.versions.first())

    MLST(SHOVILL.out.contigs)
    ch_versions = ch_versions.mix(MLST.out.versions.first())

    /*
     * The AMRFinderPlus database is fetched once and every per-sample task is
     * gated on it through the channel, making the dependency real rather than
     * a README instruction. When skipped, a placeholder file keeps the channel
     * cardinality identical so the DAG shape does not change between modes.
     */
    if (params.skip_amrfinder_update) {
        ch_db = Channel.fromPath("${projectDir}/assets/NO_DB", checkIfExists: false).ifEmpty(file('NO_DB'))
    } else {
        AMRFINDER_UPDATE()
        ch_db       = AMRFINDER_UPDATE.out.db
        ch_versions = ch_versions.mix(AMRFINDER_UPDATE.out.versions)
    }

    AMRFINDER(SHOVILL.out.contigs, ch_db.first())
    ch_versions = ch_versions.mix(AMRFINDER.out.versions.first())

    emit:
    snippy_dirs = SNIPPY.out.snippy_dir.map { id, dir -> dir }
    amrfinder   = AMRFINDER.out.report.map  { id, f   -> f   }
    mlst        = MLST.out.report.map       { id, f   -> f   }
    abricate    = ABRICATE.out.report
    qc          = ch_qc
    versions    = ch_versions
}

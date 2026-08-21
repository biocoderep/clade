include { DOWNLOAD } from '../modules/download'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_TRIM } from '../modules/fastqc'
include { FASTP } from '../modules/fastp'
include { SNIPPY } from '../modules/snippy'
include { SHOVILL } from '../modules/shovill'
include { ABRICATE } from '../modules/abricate'
include { AMRFINDER_UPDATE } from '../modules/amrfinder_update'
include { AMRFINDER } from '../modules/amrfinder'
include { MLST } from '../modules/mlst'

workflow PER_SAMPLE {
    take:
    reads_ch   // tuple(sample_id, r1, r2) -- either from DOWNLOAD or supplied directly
    reference

    main:
    FASTQC_RAW(reads_ch)
    FASTP(reads_ch)
    FASTQC_TRIM(FASTP.out.reads)

    SNIPPY(FASTP.out.reads, reference)
    SHOVILL(FASTP.out.reads)

    ABRICATE(SHOVILL.out)
    MLST(SHOVILL.out)

    // One-time database fetch (found missing during code review: AMRFinderPlus
    // hard-errors with no database on a fresh install). Runs once; every
    // sample's AMRFINDER call is gated on it via .combine(), a real dataflow
    // dependency, not just documentation.
    AMRFINDER_UPDATE()
    AMRFINDER(SHOVILL.out, AMRFINDER_UPDATE.out)

    emit:
    snippy_dirs = SNIPPY.out.map { id, dir -> dir }
    amrfinder   = AMRFINDER.out
    mlst        = MLST.out
    abricate    = ABRICATE.out
}

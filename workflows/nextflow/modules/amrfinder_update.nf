// NEW module (gap found during code review, not part of the original plan):
// AMRFinderPlus refuses to run at all without its curated database first —
// confirmed directly: `amrfinder --database_version` hard-errors with
// "No valid AMRFinder database is found" on a fresh install. The per-sample
// AMRFINDER module never handled this; without this process, the pipeline
// would fail on anyone's first real run. Runs once, not per-sample; every
// AMRFINDER call depends on its output to guarantee ordering.
process AMRFINDER_UPDATE {
    conda "${projectDir}/conda/genome_processing.yml"

    output:
    val true

    script:
    """
    amrfinder -u
    """

    stub:
    """
    true
    """
}

// Approximate-maximum-likelihood phylogeny from the core SNP alignment.
//
// FastTree is used because a full-alignment ML tree at cohort scale was not
// tractable here: two attempts on 3,255 genomes were killed by the OOM handler
// after 24-27 hours at ~500 GB. FastTree has no ascertainment-bias correction
// (unlike RAxML --asc-corr or IQ-TREE +ASC); that limitation is real and is
// why CLADE never relies on Stage 4 alone.
process FASTTREE {
    label 'process_high'

    conda "bioconda::fasttree=2.2.0"
    container "quay.io/biocontainers/fasttree:2.2.0--h7b50bb2_1"

    input:
    path alignment

    output:
    path "tree.nwk",     emit: tree
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: '-nt -gtr'
    """
    export OMP_NUM_THREADS=${task.cpus}
    FastTree ${args} ${alignment} > tree.nwk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$( FastTree -help 2>&1 | head -1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 )
    END_VERSIONS
    """

    stub:
    """
    echo "(s1:0.1,s2:0.1);" > tree.nwk
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: 2.2.0
    END_VERSIONS
    """
}

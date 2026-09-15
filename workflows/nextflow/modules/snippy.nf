// Reference-based variant calling. The whole output directory is emitted, not
// just the VCF: snippy-core requires the per-sample directory structure
// downstream.
process SNIPPY {
    tag "$sample_id"
    label 'process_medium'

    conda "bioconda::snippy=4.6.0"
    container "quay.io/biocontainers/snippy:4.6.0--hdfd78af_6"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}"), emit: snippy_dir
    path "versions.yml",                        emit: versions

    script:
    """
    snippy \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        --outdir ${sample_id} \\
        --ref ${reference} \\
        --R1 ${r1} --R2 ${r2} \\
        --force

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: \$( snippy --version 2>&1 | sed 's/snippy //' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${sample_id}
    touch ${sample_id}/snps.vcf ${sample_id}/snps.aligned.fa ${sample_id}/snps.tab
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy: 4.6.0
    END_VERSIONS
    """
}

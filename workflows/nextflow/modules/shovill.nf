// De novo assembly (SPAdes under the hood), used for the assembly-based
// typing and AMR steps that follow. Not used for variant calling — those come
// from the reference-based Snippy path, so assembly error cannot propagate
// into the genotype matrix.
process SHOVILL {
    tag "$sample_id"
    label 'process_high'

    conda "bioconda::shovill=1.4.2"
    container "quay.io/biocontainers/shovill:1.4.2--hdfd78af_1"

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}.fa"), emit: contigs
    path "versions.yml",                           emit: versions

    script:
    """
    shovill \\
        --cpus ${task.cpus} \\
        --ram ${task.memory.toGiga()} \\
        --outdir shovill_out \\
        --R1 ${r1} --R2 ${r2} \\
        --force
    mv shovill_out/contigs.fa ${sample_id}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$( shovill --version 2>&1 | sed 's/shovill //' )
    END_VERSIONS
    """

    stub:
    """
    echo ">contig1" > ${sample_id}.fa
    echo "ACGTACGTACGT" >> ${sample_id}.fa
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: 1.4.2
    END_VERSIONS
    """
}

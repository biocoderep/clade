// Fetch paired-end reads from the SRA by run accession.
// Network-bound and the most common source of transient failure in a large
// cohort run, hence the process_network label's extended retry policy.
process DOWNLOAD {
    tag "$accession"
    label 'process_network'

    conda "bioconda::sra-tools=3.4.1"
    container "quay.io/biocontainers/sra-tools:3.4.1--h4304569_1"

    input:
    val accession

    output:
    tuple val(accession), path("${accession}_1.fastq.gz"), path("${accession}_2.fastq.gz"), emit: reads
    path "versions.yml", emit: versions

    script:
    """
    prefetch --max-size u --output-directory . ${accession}
    fasterq-dump --split-files --threads ${task.cpus} --outdir . ${accession}
    gzip -f ${accession}_1.fastq ${accession}_2.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-tools: \$( fasterq-dump --version 2>&1 | grep -oP '\\d+\\.\\d+\\.\\d+' | head -1 )
    END_VERSIONS
    """

    stub:
    """
    echo -e "@r1\\nACGT\\n+\\nIIII" | gzip > ${accession}_1.fastq.gz
    echo -e "@r1\\nACGT\\n+\\nIIII" | gzip > ${accession}_2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sra-tools: 3.4.1
    END_VERSIONS
    """
}

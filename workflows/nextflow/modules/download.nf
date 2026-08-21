// Raw read retrieval — matches pipeline_orchestrator_linux_v2.py step 0 exactly
// (prefetch + fasterq-dump --split-files + gzip), reconstructed from the real
// server chain (see CLADE_Complete_Manuscript.md §5.2).
process DOWNLOAD {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    val sample_id

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    prefetch ${sample_id} --output-directory raw --max-size 100GB
    fasterq-dump raw/${sample_id}/${sample_id}.sra --split-files --outdir . --threads ${task.cpus} --force
    gzip ${sample_id}_1.fastq ${sample_id}_2.fastq
    """

    stub:
    """
    touch ${sample_id}_1.fastq.gz ${sample_id}_2.fastq.gz
    """
}

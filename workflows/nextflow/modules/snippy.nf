// Reference-based variant calling — real command, default parameters, unchanged
// from pipeline_orchestrator_linux_v2.py step 4. Output directory is needed
// whole (not just the VCF) by SNIPPY_CORE downstream.
process SNIPPY {
    tag "$sample_id"
    conda "${projectDir}/conda/genome_processing.yml"

    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}_snippy")

    script:
    """
    snippy --cpus ${task.cpus} --outdir ${sample_id}_snippy --ref ${reference} --R1 ${r1} --R2 ${r2} --force
    """

    stub:
    """
    mkdir -p ${sample_id}_snippy
    touch ${sample_id}_snippy/snps.vcf ${sample_id}_snippy/snps.aligned.fa
    """
}

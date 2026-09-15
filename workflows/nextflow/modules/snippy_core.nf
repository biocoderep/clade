// Core-genome alignment across all samples.
//
// Emits BOTH alignments deliberately. The SNP-only alignment is what the tree
// is built from; the full alignment is retained because building a tree from
// SNP-only data without ascertainment-bias correction inflates branch lengths,
// and in this project a 150-genome comparison found it also changes topology
// substantially (72% Robinson-Foulds). Keeping core.full.aln means that
// sensitivity check is reproducible rather than requiring a re-run.
process SNIPPY_CORE {
    label 'process_high_memory'

    conda "bioconda::snippy=4.6.0"
    container "quay.io/biocontainers/snippy:4.6.0--hdfd78af_6"

    input:
    path snippy_dirs
    path reference

    output:
    path "core.aln",      emit: aln
    path "core.full.aln", emit: full_aln
    path "core.vcf",      emit: vcf
    path "core.txt",      emit: stats
    path "versions.yml",  emit: versions

    script:
    """
    snippy-core --ref ${reference} ${snippy_dirs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy-core: \$( snippy-core --version 2>&1 | sed 's/snippy-core //' )
    END_VERSIONS
    """

    stub:
    """
    printf '>ref\\nACGT\\n>s1\\nACGT\\n' > core.aln
    printf '>ref\\nACGTACGT\\n>s1\\nACGTACGT\\n' > core.full.aln
    printf '##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n' > core.vcf
    touch core.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snippy-core: 4.6.0
    END_VERSIONS
    """
}

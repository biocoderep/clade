// Assemble CLADE's three inputs from pipeline outputs: the candidate genotype
// matrix (from core.vcf by exact chromosome:position match), the binary
// phenotype (from AMRFinderPlus calls for params.resistance_gene), and the
// lineage assignments (from MLST).
//
// Genotypes come from the VCF rather than from assemblies so that assembly
// error cannot enter the association test. No-calls (./., .) are preserved as
// missing rather than silently recoded as reference — conflating the two turns
// absent data into evidence of absence and biases every downstream stage the
// same way.
process BUILD_MATRIX {
    label 'process_medium'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path core_vcf
    path candidates_file
    path profiling_files
    val  resistance_gene

    output:
    path "genotypes.csv", emit: genotypes
    path "pheno.tsv",     emit: phenotype
    path "lineages.tsv",  emit: lineages
    path "versions.yml",  emit: versions

    script:
    """
    mkdir -p profiling
    cp ${profiling_files} profiling/ 2>/dev/null || true

    build_matrix.py \\
        --vcf ${core_vcf} \\
        --candidates-file ${candidates_file} \\
        --profiling-dir profiling \\
        --resistance-gene '${resistance_gene}' \\
        --outdir out
    mv out/genotypes.csv out/pheno.tsv out/lineages.tsv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python3 --version | sed 's/Python //' )
        clade: \$( python3 -c "import clade; print(clade.__version__)" )
    END_VERSIONS
    """

    stub:
    """
    echo "Sample_ID,cand_1" > genotypes.csv
    echo "s1,1"            >> genotypes.csv
    printf 'sample\\tphenotype\\ns1\\t1\\n' > pheno.tsv
    printf 'sample\\tST\\ns1\\t2\\n'        > lineages.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """
}

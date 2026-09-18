// Auto-detect the phenotype-defining resistance gene from the cohort's own
// AMRFinderPlus calls (see bin/detect_resistance_gene.py for the heuristic
// and why it is reported and overridable rather than hidden). Only runs
// when --resistance_gene is not supplied.
process DETECT_RESISTANCE_GENE {
    label 'process_low'

    conda "${projectDir}/../../environment.yml"
    container "${params.clade_container}"

    input:
    path amrfinder_reports
    val  n_samples

    output:
    path "resistance_gene.txt",            emit: gene_file
    path "resistance_gene_candidates.tsv", emit: report
    path "versions.yml",                   emit: versions

    script:
    """
    mkdir -p amr
    for f in ${amrfinder_reports}; do
        cp "\$f" amr/
    done

    detect_resistance_gene.py \\
        --amrfinder-dir amr \\
        --n-samples ${n_samples} \\
        --min-prevalence ${params.auto_gene_min_prevalence} \\
        --max-prevalence ${params.auto_gene_max_prevalence} \\
        --output resistance_gene.txt \\
        --report resistance_gene_candidates.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """

    stub:
    """
    echo "stub_gene" > resistance_gene.txt
    printf "gene\\tn_carriers\\tprevalence\\n" > resistance_gene_candidates.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clade: ${workflow.manifest.version}
    END_VERSIONS
    """
}

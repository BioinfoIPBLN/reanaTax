// SparCC, as implemented by FastSpar.
//
// Three tools in sequence, which is how FastSpar is meant to be used and why
// this is one process rather than three: the bootstrap tables are large, purely
// intermediate, and exist only to feed the p-value step in the same directory.
//
//   fastspar            the point estimate: correlations of the underlying
//                       basis, inferred from log-ratio variances rather than
//                       from the closed proportions.
//   fastspar_bootstrap  resampled count tables under the null.
//   fastspar_pvalues    where the real correlation falls in that null.
//
// The bootstrap loop is serial in the script but each fastspar call is
// threaded, so the task's CPUs are used throughout. `--iterations 5` on the
// bootstraps is FastSpar's own recommendation: the null does not need the
// convergence the point estimate does, and running it to 50 would multiply the
// cost of the whole step by ten for no change in the p-values.
process SPARCC_CORRELATION {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastspar:1.0.0--h1b620e3_6'
        : 'quay.io/biocontainers/fastspar:1.0.0--h1b620e3_6'}"

    input:
    tuple val(meta), path(otu)
    val iterations
    val exclusion_iterations
    val permutations

    output:
    tuple val(meta), path("*.sparcc_correlation.tsv"), emit: correlation
    tuple val(meta), path("*.sparcc_covariance.tsv"), emit: covariance
    tuple val(meta), path("*.sparcc_pvalues.tsv"), emit: pvalues
    tuple val("${task.process}"), val('fastspar'), eval("fastspar --version 2>&1 | sed -n 's/^.*fastspar[ v]*//p' | head -n1"), emit: versions_fastspar, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastspar \\
        ${args} \\
        --otu_table ${otu} \\
        --correlation ${prefix}.sparcc_correlation.tsv \\
        --covariance ${prefix}.sparcc_covariance.tsv \\
        --iterations ${iterations} \\
        --exclusion_iterations ${exclusion_iterations} \\
        --threads ${task.cpus}

    mkdir -p bootstrap_counts bootstrap_correlation

    fastspar_bootstrap \\
        --otu_table ${otu} \\
        --number ${permutations} \\
        --prefix bootstrap_counts/boot \\
        --threads ${task.cpus}

    for table in bootstrap_counts/boot_*.tsv; do
        name=\$(basename "\${table}" .tsv)
        fastspar \\
            --otu_table "\${table}" \\
            --correlation "bootstrap_correlation/cor_\${name}.tsv" \\
            --covariance "bootstrap_correlation/cov_\${name}.tsv" \\
            --iterations 5 \\
            --threads ${task.cpus} \\
            > /dev/null
    done

    fastspar_pvalues \\
        --otu_table ${otu} \\
        --correlation ${prefix}.sparcc_correlation.tsv \\
        --prefix bootstrap_correlation/cor_boot_ \\
        --permutations ${permutations} \\
        --outfile ${prefix}.sparcc_pvalues.tsv \\
        --threads ${task.cpus}

    # The resampled tables are the size of the input times --permutations and
    # have no use once the p-values exist.
    rm -rf bootstrap_counts bootstrap_correlation
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sparcc_correlation.tsv ${prefix}.sparcc_covariance.tsv ${prefix}.sparcc_pvalues.tsv
    """
}

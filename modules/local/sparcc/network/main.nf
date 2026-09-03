// FastSpar's two square matrices as one edge table, with the p-values adjusted.
//
// The adjustment is why this is a step rather than a formatting detail. A
// 200-taxon network is 19,900 pairs, so at an unadjusted 0.05 about a thousand
// edges are expected from nothing at all - which is how co-occurrence networks
// acquire their reputation. Benjamini-Hochberg across every pair at once.
process SPARCC_NETWORK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(correlation), path(pvalues), path(taxa)
    val permutations
    val min_correlation
    val p_threshold

    output:
    tuple val(meta), path("*.sparcc_edges_all.tsv"), emit: edges
    tuple val(meta), path("*.sparcc_edges_significant.tsv"), emit: significant
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sparcc_network.py \\
        ${args} \\
        --correlation ${correlation} \\
        --pvalues ${pvalues} \\
        --taxa ${taxa} \\
        --permutations ${permutations} \\
        --min-correlation ${min_correlation} \\
        --p-threshold ${p_threshold} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sparcc_edges_all.tsv ${prefix}.sparcc_edges_significant.tsv ${prefix}_sparcc_mqc.tsv
    """
}

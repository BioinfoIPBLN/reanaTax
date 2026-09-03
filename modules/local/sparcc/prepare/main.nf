// The combined Bracken table as an OTU table, and the feasibility check.
//
// Both jobs in one step on purpose. The conversion is trivial; the question of
// whether this cohort can support a co-occurrence network at all is not, and
// asking it anywhere later would mean the run has already been paid for. See
// bin/sparcc_prepare.py for what is refused and why.
process SPARCC_PREPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(counts)
    val min_prevalence
    val min_reads
    val top_taxa
    val permutations
    val p_threshold
    val force

    output:
    tuple val(meta), path("*.sparcc_otu.tsv"), emit: otu
    tuple val(meta), path("*.sparcc_taxa.tsv"), emit: taxa
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def force_arg = force ? '--force' : ''
    """
    sparcc_prepare.py \\
        ${args} \\
        --counts ${counts} \\
        --min-prevalence ${min_prevalence} \\
        --min-reads ${min_reads} \\
        --top-taxa ${top_taxa} \\
        --permutations ${permutations} \\
        --p-threshold ${p_threshold} \\
        ${force_arg} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sparcc_otu.tsv ${prefix}.sparcc_taxa.tsv
    """
}

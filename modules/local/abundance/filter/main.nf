// Drop taxa that never reach a meaningful relative abundance.
//
// Not a contamination caller - see bin/filter_abundance.py for why that
// distinction matters. The rows it removes are written out alongside the
// filtered table so the decision stays auditable.
process ABUNDANCE_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(table, stageAs: 'input/*')
    path drop_list
    val min_rel_abundance
    val min_samples
    val min_reads

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: filtered
    tuple val(meta), path("*.removed.tsv"), emit: removed
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A collection, because the minimizer and host-k-mer filters each write
    // their own list and a taxon condemned by either is removed.
    def drop_files = (drop_list instanceof List ? drop_list : [drop_list]).findAll { entry -> entry }
    def drop_arg = drop_files ? "--drop-taxids ${drop_files.join(' ')}" : ''
    """
    filter_abundance.py \\
        ${args} \\
        ${drop_arg} \\
        --input ${table} \\
        --output ${prefix}.filtered.tsv \\
        --removed ${prefix}.removed.tsv \\
        --min-rel-abundance ${min_rel_abundance} \\
        --min-samples ${min_samples} \\
        --min-reads ${min_reads}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.tsv ${prefix}.removed.tsv
    """
}

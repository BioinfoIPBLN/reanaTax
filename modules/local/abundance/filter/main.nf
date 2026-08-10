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
    val min_rel_abundance
    val min_samples

    output:
    tuple val(meta), path("*.filtered.tsv"), emit: filtered
    tuple val(meta), path("*.removed.tsv"), emit: removed
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_abundance.py \\
        ${args} \\
        --input ${table} \\
        --output ${prefix}.filtered.tsv \\
        --removed ${prefix}.removed.tsv \\
        --min-rel-abundance ${min_rel_abundance} \\
        --min-samples ${min_samples}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.tsv ${prefix}.removed.tsv
    """
}

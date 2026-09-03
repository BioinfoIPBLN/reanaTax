//
// Per-sample featureCounts tables into one gene x sample matrix.
//
// featureCounts names its count column after the BAM it was handed, so the
// tables cannot simply be pasted together - see bin/merge_featurecounts.py.
//
process HOSTCOUNTS_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(tables)

    output:
    tuple val(meta), path("*.host_counts.tsv") , emit: counts
    tuple val(meta), path("*.gene_lengths.tsv"), emit: lengths
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_featurecounts.py ${tables} \\
        --output ${prefix}.host_counts.tsv \\
        --lengths ${prefix}.gene_lengths.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.host_counts.tsv ${prefix}.gene_lengths.tsv
    """
}

//
// Host gene expression against microbial abundance, across samples.
//
// Monteleone et al.'s third layer, and the one that justifies --quantify_host:
// both sides come from the SAME library, so an association between them is not
// confounded by sample handling the way two separate assays would be.
//
// The script refuses to run when the cohort is too small for any pair to reach
// significance - see the feasibility guard in bin/host_microbe_correlation.py,
// which is the substance of the module rather than a safety rail.
//
process HOSTMICROBE_CORRELATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(host_counts), path(microbial)

    output:
    tuple val(meta), path("*.host_microbe_correlation.tsv"), emit: results
    tuple val(meta), path("*_mqc.tsv")                     , emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    host_microbe_correlation.py \\
        --host ${host_counts} \\
        --microbial ${microbial} \\
        --output ${prefix}.host_microbe_correlation.tsv \\
        --mqc ${prefix}.host_microbe_correlation_mqc.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.host_microbe_correlation.tsv ${prefix}.host_microbe_correlation_mqc.tsv
    """
}

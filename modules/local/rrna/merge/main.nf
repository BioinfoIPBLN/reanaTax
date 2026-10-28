// The cohort's rRNA share per taxon, with the minimizer filter's verdict,
// coverage and duplication beside it when that filter ran.
process RRNA_SHARE_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/../../read/accounting/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(tables, stageAs: 'samples/*')
    path evidence

    output:
    tuple val(meta), path("*.rrna_share.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def evidence_arg = evidence ? "--evidence ${evidence}" : ''
    """
    rrna_share.py merge \\
        ${args} \\
        --tables ${tables} \\
        ${evidence_arg} \\
        --output ${prefix}.rrna_share.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rrna_share.tsv
    """
}

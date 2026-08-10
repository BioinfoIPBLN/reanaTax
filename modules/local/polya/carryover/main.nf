// Ask whether the microbial reads in a poly(A)-selected library got there by
// oligo-dT internal mispriming or by non-specific carry-over, by comparing
// internal poly-A/T runs against a composition-matched permutation null.
process POLYA_CARRYOVER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.polya.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    polya_carryover.py \\
        ${args} \\
        --fastq ${reads} \\
        --sample ${prefix} \\
        --output ${prefix}.polya.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.polya.tsv
    """
}

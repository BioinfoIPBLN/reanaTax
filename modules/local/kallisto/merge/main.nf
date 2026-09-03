//
// Per-sample kallisto abundances into one counts matrix, in the same shape
// HOSTCOUNTS_MERGE produces - so one differential-expression path serves both
// quantifiers and there is no second copy of it to keep correct.
//
process KALLISTO_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'quay.io/biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(abundances, stageAs: 'abundance/*')
    path tx2gene

    output:
    tuple val(meta), path("*.counts.tsv") , emit: counts
    tuple val(meta), path("*.tpm.tsv")    , emit: tpm
    tuple val(meta), path("*.lengths.tsv"), emit: lengths
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def map_arg = tx2gene ? "--tx2gene ${tx2gene}" : ''
    """
    merge_kallisto.py abundance/* \\
        --counts ${prefix}.counts.tsv \\
        --tpm ${prefix}.tpm.tsv \\
        --lengths ${prefix}.lengths.tsv \\
        ${map_arg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv ${prefix}.tpm.tsv ${prefix}.lengths.tsv
    """
}

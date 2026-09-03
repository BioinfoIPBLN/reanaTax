// Join fastp, SortMeRNA, HISAT2 and Kraken2 into a single per-sample read-fate
// table, and derive the host carry-over metric from it.
//
// One task for the whole cohort rather than one per sample: the inputs are a
// few kilobytes each, and the table is only useful when every sample is in it.
process READ_ACCOUNTING {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(fastp_json, stageAs: 'fastp/*'), path(sortmerna_logs, stageAs: 'sortmerna/*'), path(hisat2_logs, stageAs: 'hisat2/*'), path(kraken2_reports, stageAs: 'kraken2/*')
    val host_taxid

    output:
    tuple val(meta), path("*.read_accounting.tsv"), emit: tsv
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastp_arg = fastp_json ? "--fastp ${fastp_json}" : ''
    def sortmerna_arg = sortmerna_logs ? "--sortmerna ${sortmerna_logs}" : ''
    def hisat2_arg = hisat2_logs ? "--hisat2 ${hisat2_logs}" : ''
    def kraken2_arg = kraken2_reports ? "--kraken2 ${kraken2_reports}" : ''
    """
    read_accounting.py \\
        ${args} \\
        ${fastp_arg} \\
        ${sortmerna_arg} \\
        ${hisat2_arg} \\
        ${kraken2_arg} \\
        --host-taxid ${host_taxid} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.read_accounting.tsv
    """
}

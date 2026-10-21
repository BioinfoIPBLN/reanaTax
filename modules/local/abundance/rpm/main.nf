// A combined table restated as reads per million of the whole library.
//
// Every fraction the pipeline reports is a share of the classified or
// surviving reads, which in host RNA-seq hides the scale: one read can be 1% of
// what is left. This divides each count by the sample's trimmed read total
// from the read accounting instead. A unit, not a filter - see
// bin/rpm_total.py.
process ABUNDANCE_RPM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(table)
    path accounting

    output:
    tuple val(meta), path("*.rpm_total.tsv"), emit: rpm
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: table.baseName
    """
    rpm_total.py \\
        --table ${table} \\
        --accounting ${accounting} \\
        --output ${prefix}.rpm_total.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: table.baseName
    """
    touch ${prefix}.rpm_total.tsv
    """
}

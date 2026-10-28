// One library's rRNA share per taxon: the rRNA read IDs joined to Kraken2's
// per-read calls and rolled up the library's own report. See bin/rrna_share.py.
process RRNA_SHARE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/../../read/accounting/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(ids), path(assignments), path(report)

    output:
    tuple val(meta), path("*.rrna_share.tsv"), emit: tsv
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rrna_share.py sample \\
        --sample ${meta.id} \\
        --ids ${ids} \\
        --assignments ${assignments} \\
        --report ${report} \\
        --output ${prefix}.rrna_share.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'sample\\ttaxid\\trank\\tname\\treads\\trrna_reads\\n' > ${prefix}.rrna_share.tsv
    """
}

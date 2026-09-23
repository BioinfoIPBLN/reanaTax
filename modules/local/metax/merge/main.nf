// Cohort tables from the per-sample Metax profiles: every sample's rows in one
// long table, and taxon-by-sample tables of abundance and reads.
//
// Metax writes its profile without a header row, so bin/metax_merge.py is also
// where the columns get their names. A sample in which nothing passed Metax's
// coverage filters stays in the wide tables as a column of zeros rather than
// dropping out of them.
process METAX_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(profiles, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.long.tsv")                        , emit: combined
    tuple val(meta), path("*.{abundance,reads}.tsv")           , emit: merged
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    metax_merge.py \\
        --profiles profiles/*.profile.txt \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.long.tsv ${prefix}.abundance.tsv ${prefix}.reads.tsv
    """
}

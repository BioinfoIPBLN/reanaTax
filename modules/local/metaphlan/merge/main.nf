// One table for the cohort, from the per-sample MetaPhlAn profiles.
//
// merge_metaphlan_tables.py ships with MetaPhlAn and keys the columns off the
// input file names, so the profiles are staged under `<sample>_metaphlan.txt`
// upstream and the column headers come out as sample ids.
process METAPHLAN_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metaphlan:4.2.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/metaphlan:4.2.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(profiles, stageAs: 'profiles/*')

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version | sed -n 's/.*version //p'"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_metaphlan_tables.py profiles/*_metaphlan.txt > ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    """
}

// Each species' core gene families, the MetaPhlAn marker candidates, and the
// genomes the uniqueness test maps them against. See bin/mpa_candidates.py.
process DB_MPA_CANDIDATES {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(clusters), path(genes, stageAs: 'genes/*'), path(taxonomy)

    output:
    tuple val(meta), path("candidates.fna"), path("candidates.tsv"), path("targets.tsv"), path("species.tsv"), emit: candidates
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mpa_candidates.py \\
        ${args} \\
        --clusters ${clusters} \\
        --genes genes/*.genes.tsv \\
        --genomes genes/*.genome.tsv \\
        --sequences genes/*.genes.fna \\
        --taxonomy ${taxonomy}
    """

    stub:
    """
    touch candidates.fna candidates.tsv targets.tsv species.tsv
    """
}

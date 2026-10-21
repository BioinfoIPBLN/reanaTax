// One genome of the --db_input sheet, with every header carrying its taxid
// (`kraken:taxid|<n>`) and a seqid-to-taxid map beside it: the form both
// kraken2-build and krakenuniq-build accept without NCBI's accession tables.
// See bin/db_prepare_genome.py.
process DB_PREPARE_GENOME {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(fasta, stageAs: 'input/*')

    output:
    tuple val(meta), path("${meta.id}.fna"), path("${meta.id}.fna.map"), emit: genome
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    db_prepare_genome.py \\
        --fasta ${fasta} \\
        --id ${meta.id} \\
        --taxid ${meta.taxid}
    """

    stub:
    """
    touch ${meta.id}.fna ${meta.id}.fna.map
    """
}

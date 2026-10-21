// Every genome's taxid checked against the taxonomy before anything is built,
// and the lineage of each written as the metadata sylph-tax attaches to sylph's
// genome-level calls. See bin/db_lineage.py.
process DB_LINEAGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(genomes)
    path taxonomy, stageAs: 'taxonomy'

    output:
    tuple val(meta), path("*.sylph_taxonomy.tsv"), emit: sylph_taxonomy
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    db_lineage.py \\
        --taxonomy taxonomy \\
        --genomes ${genomes} \\
        --output ${prefix}.sylph_taxonomy.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sylph_taxonomy.tsv
    """
}

// One genome's protein-coding genes of marker length, for the MetaPhlAn build:
// from the sheet's cds file when there is one, from pyrodigal otherwise. See
// bin/mpa_genes.py.
process DB_MPA_GENES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pyrodigal:3.7.1--py312h247cb63_1'
        : 'quay.io/biocontainers/pyrodigal:3.7.1--py312h247cb63_1'}"

    input:
    tuple val(meta), path(fasta), path(cds, stageAs: 'cds/*')

    output:
    tuple val(meta), path("${meta.id}.genes.fna"), path("${meta.id}.genes.faa"), path("${meta.id}.genes.tsv"), path("${meta.id}.genome.tsv"), emit: genes
    tuple val("${task.process}"), val('pyrodigal'), eval('python3 -c "import pyrodigal; print(pyrodigal.__version__)"'), emit: versions_pyrodigal, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cds_arg = cds ? "--cds ${cds}" : ''
    """
    mpa_genes.py \\
        ${args} \\
        --fasta ${fasta} \\
        --id ${meta.id} \\
        --taxid ${meta.taxid} \\
        ${cds_arg}
    """

    stub:
    """
    touch ${meta.id}.genes.fna ${meta.id}.genes.faa ${meta.id}.genes.tsv ${meta.id}.genome.tsv
    """
}

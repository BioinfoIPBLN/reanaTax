// The MetaPhlAn database: each species' markers chosen, the pickle MetaPhlAn
// reads written, and the markers indexed with the Bowtie2 of the MetaPhlAn
// container the analysis runs (4.2.6), as a large index - MetaPhlAn looks for
// six *.bt2l files beside <name>.pkl. See bin/mpa_select.py.
process DB_METAPHLAN_BUILD {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metaphlan:4.2.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/metaphlan:4.2.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(candidates_fna), path(candidates_tsv), path(targets), path(species), path(hits), path(taxonomy)

    output:
    // The database's files and folders, not the `db` directory itself:
    // publishDir skips everything under a directory output, so the
    // `db/...` patterns in conf/modules.config would publish nothing.
    tuple val(meta), path("db/*"), emit: db
    tuple val(meta), path("*.mpa_species.tsv"), emit: summary
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version | sed -n 's/.*version //p'"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def name = task.ext.prefix ?: "${meta.id}"
    """
    mpa_select.py \\
        ${args} \\
        --candidates ${candidates_tsv} \\
        --sequences ${candidates_fna} \\
        --targets ${targets} \\
        --species ${species} \\
        --hits ${hits} \\
        --taxonomy ${taxonomy} \\
        --name ${name}

    mkdir db
    mv ${name}.pkl ${name}.fna db/
    bowtie2-build --large-index --threads ${task.cpus} db/${name}.fna db/${name}
    """

    stub:
    def name = task.ext.prefix ?: "${meta.id}"
    """
    mkdir db
    touch db/${name}.pkl db/${name}.fna ${name}.mpa_species.tsv
    for part in 1 2 3 4 rev.1 rev.2; do touch db/${name}.\${part}.bt2l; done
    """
}

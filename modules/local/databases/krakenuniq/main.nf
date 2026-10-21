// A KrakenUniq database from the genomes of the --db_input sheet.
//
// krakenuniq-build reads every FASTA in the library directory and the `.map`
// files beside them (sequence id, taxid), which DB_PREPARE_GENOME writes for
// each genome, so no accession lookup is needed. The container is the one the
// KrakenUniq classification module uses, which carries krakenuniq-build and
// jellyfish.
process DB_KRAKENUNIQ_BUILD {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakenuniq:1.0.4--pl5321h6dccd9a_2'
        : 'quay.io/biocontainers/krakenuniq:1.0.4--pl5321h6dccd9a_2'}"

    input:
    tuple val(meta), path(genomes, stageAs: 'genomes/*')
    path taxonomy, stageAs: 'taxonomy'

    output:
    tuple val(meta), path("db"), emit: db
    tuple val("${task.process}"), val('krakenuniq'), eval("krakenuniq --version 2>&1 | head -1 | sed 's/KrakenUniq version //'"), emit: versions_krakenuniq, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p db/library db/taxonomy
    cp taxonomy/names.dmp taxonomy/nodes.dmp db/taxonomy/
    for file in genomes/*; do
        cp -L "\$file" db/library/
    done

    krakenuniq-build --db db --threads ${task.cpus} ${args}
    """

    stub:
    """
    mkdir -p db
    touch db/database.kdb db/database.idx db/taxDB
    """
}

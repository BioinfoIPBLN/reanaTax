// A Kraken2 database from the genomes of the --db_input sheet.
//
// The same Kraken2 the classification modules run (2.17.2), so the database is
// read by the version that built it. Every genome is added to the library with
// its taxid in the headers (see DB_PREPARE_GENOME), masked for low complexity
// as kraken2-build does by default, and the index is built against the NCBI
// taxonomy given. `kraken2-inspect` then writes inspect.txt, the per-clade
// minimizer counts the minimizer filter divides by for its coverage figure: a
// database without it gets no coverage evidence.
//
// The library stays in the task directory, not in the published database:
// DB_BRACKEN_BUILD classifies it to estimate Bracken's distributions, and
// nothing reads it afterwards.
process DB_KRAKEN2_BUILD {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0852ebcb1e7341e3bda5be3cde7fbc92d0f69af04448e877a36d55bc331013ae/data'
        : 'community.wave.seqera.io/library/kraken2_bracken_coreutils_findutils_pruned:64b27b48dbc9e832'}"

    input:
    tuple val(meta), path(genomes, stageAs: 'genomes/*')
    path taxonomy, stageAs: 'taxonomy'

    output:
    // The database's files and folders, not the `db` directory itself:
    // publishDir skips everything under a directory output, so the
    // `db/...` patterns in conf/modules.config would publish nothing.
    tuple val(meta), path("db/*"), emit: db
    tuple val("${task.process}"), val('kraken2'), eval("kraken2 --version | head -1 | sed 's/Kraken version //'"), emit: versions_kraken2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def library_args = task.ext.args2 ?: ''
    """
    mkdir -p db/taxonomy
    cp taxonomy/names.dmp taxonomy/nodes.dmp db/taxonomy/

    for genome in genomes/*.fna; do
        kraken2-build --add-to-library "\$genome" --db db --threads ${task.cpus} ${library_args}
    done

    kraken2-build --build --db db --threads ${task.cpus} ${args}
    kraken2-inspect --db db --threads ${task.cpus} > db/inspect.txt
    """

    stub:
    """
    mkdir -p db
    touch db/hash.k2d db/opts.k2d db/taxo.k2d db/inspect.txt
    """
}

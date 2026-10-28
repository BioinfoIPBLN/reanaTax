// Bracken's k-mer distributions for a database DB_KRAKEN2_BUILD made, one per
// read length in --db_bracken_read_lengths.
//
// One task for every length: the expensive part of bracken-build is classifying
// the whole library with Kraken2, which it does once and reuses for the other
// lengths. It writes into the database directory, so the database is linked
// file by file into a directory of this task's own rather than written through
// the staged input. The k-mer length is read from the database's inspect.txt,
// because Bracken's -k has to match the one the database was built with.
process DB_BRACKEN_BUILD {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/08/0852ebcb1e7341e3bda5be3cde7fbc92d0f69af04448e877a36d55bc331013ae/data'
        : 'community.wave.seqera.io/library/kraken2_bracken_coreutils_findutils_pruned:64b27b48dbc9e832'}"

    input:
    tuple val(meta), path(db, stageAs: 'input_db/*')
    val read_lengths

    output:
    tuple val(meta), path("bracken/database*mers.kmer_distrib"), emit: distributions
    tuple val("${task.process}"), val('bracken'), eval("bracken -v 2>&1 | head -1 | sed 's/Bracken v//'"), emit: versions_bracken, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir bracken
    for entry in input_db/*; do
        ln -s "\$(readlink -f "\$entry")" bracken/
    done
    kmer_len=\$(sed -n 's/.*k = \\([0-9]*\\).*/\\1/p' bracken/inspect.txt | head -1)

    for length in ${read_lengths.join(' ')}; do
        bracken-build -d bracken -t ${task.cpus} -k "\${kmer_len:-35}" -l "\$length" ${args}
    done
    """

    stub:
    """
    mkdir bracken
    for length in ${read_lengths.join(' ')}; do
        touch bracken/database\${length}mers.kmer_distrib
    done
    """
}

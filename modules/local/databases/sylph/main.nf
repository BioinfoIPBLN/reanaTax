// A sylph database: the genomes of the --db_input sheet, sketched.
//
// sylph records each genome by the path it was given, and sylph-tax attaches a
// lineage to it by that file name (DB_LINEAGE writes the matching metadata), so
// the genomes are sketched under the names `<id>.fna`. The default -c 200
// matches sylph's prebuilt databases; pass `--db_sylph_args '-c 100'` for a
// denser one, and the same -c to `--sylph_args` when profiling against it.
process DB_SYLPH_SKETCH {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:1.0.0--hb42e459_0'
        : 'quay.io/biocontainers/sylph:1.0.0--hb42e459_0'}"

    input:
    tuple val(meta), path(genomes)

    output:
    tuple val(meta), path("*.syldb"), emit: db
    tuple val("${task.process}"), val('sylph'), eval("sylph --version 2>&1 | sed 's/sylph //'"), emit: versions_sylph, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sylph sketch \\
        -g ${genomes} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.syldb
    """
}

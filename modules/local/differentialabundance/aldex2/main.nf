// Differential abundance with ALDEx2, on the combined Bracken count table.
//
// The method is fixed per process because each one needs its own container:
// there is no biocontainer carrying both ALDEx2 and ANCOMBC. See
// bin/differential_abundance.R for why this takes counts and not fractions,
// and ezRun's EzAppDiffShot for the filter and contrast conventions it follows.
process DIFFERENTIAL_ABUNDANCE_ALDEX2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-aldex2:1.38.0--r44hdfd78af_0'
        : 'quay.io/biocontainers/bioconductor-aldex2:1.38.0--r44hdfd78af_0'}"

    input:
    tuple val(meta), path(counts)
    path metadata

    output:
    tuple val(meta), path("*.results.tsv"), emit: results
    tuple val(meta), path("*.volcano.png"), emit: volcano, optional: true
    tuple val(meta), path("*_mqc.tsv")    , emit: mqc, optional: true
    tuple val("${task.process}"), val('ALDEx2'), eval("Rscript -e 'cat(as.character(packageVersion(\"ALDEx2\")))'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    differential_abundance.R \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --method aldex2 \\
        --prefix ${prefix} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.aldex2.results.tsv
    """
}

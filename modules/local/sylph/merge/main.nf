// Cohort tables from the per-sample sylph and sylph-tax outputs: every
// sample's genome calls in one file, and one clade-by-sample table per
// abundance when there is a taxonomy.
//
// bin/sylph_merge.py rather than `sylph-tax merge`, for one reason: a sample in
// which sylph detected nothing has to stay in the table as a column of zeros.
// Dropping it would make an empty library look like one that was never
// profiled. The layout is the one `sylph-tax merge` writes.
process SYLPH_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/sylph-tax:1.9.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(profiles, stageAs: 'profiles/*'), path(taxprofs, stageAs: 'taxprof/*')

    output:
    tuple val(meta), path("*.genomes.tsv")  , emit: genomes
    tuple val(meta), path("*_abundance.tsv"), emit: merged, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def taxprof_arg = taxprofs ? '--taxprof taxprof/*.sylphmpa' : ''
    """
    sylph_merge.py \\
        --profiles profiles/*.sylph.tsv \\
        ${taxprof_arg} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.genomes.tsv
    """
}

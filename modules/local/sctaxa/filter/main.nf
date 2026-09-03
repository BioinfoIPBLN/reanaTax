// The evidence filters' drop lists, applied to the cell-by-taxon matrix.
//
// Until this existed the drop lists reached only the cohort tables - the
// combined Kraken2 report and the combined Bracken table - so on a single-cell
// run --minimizer_filter and --shuffle_control changed the abundance tables and
// left the matrix the enrichment test reads untouched.
//
// The enrichment test is where it matters most. Its correction is levied across
// every taxon x cell type pair, so every background taxon that survives costs
// power for the real ones. See bin/sc_filter_taxa.py.
//
// One task for the cohort: a drop list is a cohort-level verdict, and the
// matrix is one table.
process SCTAXA_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(counts, stageAs: 'counts/*')
    path drop_lists, stageAs: 'drop/*'

    output:
    tuple val(meta), path("*.cell_taxa.tsv"), emit: counts
    tuple val(meta), path("*.cell_taxa_removed.tsv"), emit: removed
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def drop = drop_lists ? "--drop drop/*" : ''
    """
    sc_filter_taxa.py \\
        ${args} \\
        --counts counts/* \\
        ${drop} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cell_taxa.tsv ${prefix}.cell_taxa_removed.tsv ${prefix}_sc_filter_mqc.tsv
    """
}

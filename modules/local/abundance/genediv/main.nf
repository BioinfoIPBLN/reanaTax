// Does a taxon's evidence spread over its genome, or pile onto one locus?
//
// PRISM's argument, applied to a table this pipeline can already produce. Every
// filter upstream counts reads, or the k-mers behind reads; none of them can
// see WHERE those reads landed, and that is the whole difference between a
// species that is present and a species whose name has been attached to one
// conserved stretch of sequence. HUMAnN's species-stratified gene-family table
// answers it directly, so no new reference and no new alignment are needed.
//
// Run twice: once on the gene-family table and once on the same table regrouped
// onto KEGG orthologs, which is the difference between "many gene families" and
// "many products". A taxon can reach dozens of UniRef90 families that all
// encode the same thing.
//
// Unlike the other evidence filters this one cannot feed ABUNDANCE_FILTER
// inside the taxonomy subworkflow - HUMAnN runs downstream of it, and a taxid
// list cannot travel backwards - so its drop list is applied by a further
// filter pass at the workflow level.
process ABUNDANCE_GENEDIV {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(genefamilies, stageAs: 'tables/*')
    path bracken
    val label
    val min_families
    val max_top_fraction
    val min_abundance

    output:
    tuple val(meta), path("*_diversity.tsv"), emit: evidence
    tuple val(meta), path("*_diversity_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bracken_arg = bracken ? "--bracken ${bracken}" : ''
    """
    gene_diversity_filter.py \\
        ${args} \\
        --genefamilies tables/* \\
        ${bracken_arg} \\
        --label ${label} \\
        --min-genes ${min_families} \\
        --max-top-fraction ${max_top_fraction} \\
        --min-abundance ${min_abundance} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${label}_diversity.tsv ${prefix}.${label}_diversity_drop.txt \\
        ${prefix}_${label}_diversity_mqc.tsv
    """
}

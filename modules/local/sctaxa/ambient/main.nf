// Is a taxon associated with cells, or is it in the medium?
//
// An ANNOTATION, not a filter, and the distinction is the point. An empty
// droplet holds reagent contaminants, ambient soup, microbes released by lysed
// cells AND genuine extracellular organisms from the tissue - a luminal or
// mucosal organism is in that pool by construction, not by artefact. So a
// taxon as common in empty droplets as in cells is reported as UNDECIDED, not
// as contamination. Only an external measurement of the kit settles that, and
// DECONTAM_FILTER's blanks are the only thing here that supplies one.
//
// Needs both halves of the Solo.out tree: filtered/ names the cells, raw/
// names every barcode including the empty droplets. The raw matrix is only
// ever streamed for its column sums, never held.
//
// One task per sample, because a barcode means nothing outside its own library.
process SCTAXA_AMBIENT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(solo), path(cell_taxa)
    val feature_type
    val min_umis
    val min_empty_umis
    val max_empty_umis
    val min_droplets
    val enriched_ratio
    val depleted_ratio
    val p_threshold
    val drop_ambient

    output:
    tuple val(meta), path("*.sc_ambient.tsv"), emit: evidence
    tuple val(meta), path("*.sc_ambient_drop.txt"), emit: drop_list, optional: true
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def drop_arg = drop_ambient ? '--drop-ambient' : ''
    """
    sc_ambient.py \\
        ${args} \\
        --matrix ${solo} \\
        --feature-type ${feature_type} \\
        --cell-taxa ${cell_taxa} \\
        --sample ${prefix} \\
        --min-umis ${min_umis} \\
        --min-empty-umis ${min_empty_umis} \\
        --max-empty-umis ${max_empty_umis} \\
        --min-droplets ${min_droplets} \\
        --enriched-ratio ${enriched_ratio} \\
        --depleted-ratio ${depleted_ratio} \\
        --p-threshold ${p_threshold} \\
        ${drop_arg} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sc_ambient.tsv ${prefix}_sc_ambient_mqc.tsv
    """
}

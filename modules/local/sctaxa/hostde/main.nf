// The host half of SAHMI: what a cell's transcriptome does when it is carrying
// something.
//
// SAHMI's cell-by-taxon matrix is not the result, it is the instrument. The
// result is the comparison it makes possible - the host transcriptome of cells
// carrying a taxon against bystander cells OF THE SAME TYPE - and this module
// is the step that had been missing here.
//
// Both halves of that comparison come out of the same STARsolo pass: the
// cell-by-gene matrix is STARSOLO's Solo.out, and the cell-by-taxon matrix was
// built from the reads STARsolo could not place, keyed on the same corrected
// barcode. So the two views cannot disagree about which cell is which.
//
// One task per sample, because a barcode means nothing outside its own library.
process SCTAXA_HOSTDE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(solo), path(cell_taxa)
    path cells
    val feature_type
    val min_umis
    val min_cells
    val min_pct
    val logfc_threshold
    val top_taxa
    val p_threshold
    val force_pooled

    output:
    tuple val(meta), path("*.sc_host_de.tsv"), emit: de
    tuple val(meta), path("*.sc_host_de_groups.tsv"), emit: groups
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def cells_arg = cells ? "--cells ${cells}" : ''
    def pooled_arg = force_pooled ? '--force-pooled' : ''
    """
    sc_host_de.py \\
        ${args} \\
        --matrix ${solo} \\
        --feature-type ${feature_type} \\
        --cell-taxa ${cell_taxa} \\
        --sample ${prefix} \\
        ${cells_arg} \\
        ${pooled_arg} \\
        --min-umis ${min_umis} \\
        --min-cells ${min_cells} \\
        --min-pct ${min_pct} \\
        --logfc-threshold ${logfc_threshold} \\
        --top-taxa ${top_taxa} \\
        --p-threshold ${p_threshold} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sc_host_de.tsv ${prefix}.sc_host_de_groups.tsv ${prefix}_sc_host_de_mqc.tsv
    """
}

// The cell-by-taxon matrix for plate-based single-cell data.
//
// In Smart-seq2 and its relatives every well is its own library, so there is no
// barcode to recover and no barcode-aware aligner to run: each cell has already
// been classified by the ordinary bulk route, and the matrix is those reports
// stacked. CSI-Microbes treats this layout as a separate route for exactly this
// reason, and so does this pipeline.
//
// One task for the whole cohort, unlike SCTAXA_COUNTS - the input is one small
// report per cell rather than one huge read-assignment file per library.
//
// Emits the same long format SCTAXA_COUNTS does, so SCTAXA_ENRICHMENT cannot
// tell the two routes apart.
process SCTAXA_PLATE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')
    path metadata
    val ranks
    val min_reads

    output:
    tuple val(meta), path("*.cell_taxa.tsv"), emit: counts
    tuple val(meta), path("*.cell_taxa_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sc_plate_matrix.py \\
        ${args} \\
        --reports reports/* \\
        --metadata ${metadata} \\
        --ranks '${ranks}' \\
        --min-reads ${min_reads} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cell_taxa.tsv ${prefix}.cell_taxa_summary.tsv
    """
}

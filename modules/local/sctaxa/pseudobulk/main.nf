// The host response to a taxon, tested across donors rather than across cells.
//
// SCTAXA_HOSTDE compares infected with bystander cells inside one library, with
// cells as the units: a statement about that library. Here the per-library sums
// it writes under --sc_host_de_pseudobulk are set side by side, one matrix per
// taxon x cell type with a column per donor and state, for a paired DESeq2 /
// edgeR test in which the donor is the unit of replication. See
// bin/sc_pseudobulk.py.
//
// One task for the cohort: a donor-level test needs every donor at once.
process SCTAXA_PSEUDOBULK {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(tables)
    val min_donors
    val min_count
    val min_samples

    output:
    tuple val(meta), path("*.pseudobulk_counts.tsv"), emit: counts, optional: true
    tuple val(meta), path("*.pseudobulk_metadata.tsv"), emit: metadata
    tuple val(meta), path("*.pseudobulk_groups.tsv"), emit: groups
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sc_pseudobulk.py \\
        --tables ${tables} \\
        --min-donors ${min_donors} \\
        --min-count ${min_count} \\
        --min-samples ${min_samples} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pseudobulk_metadata.tsv ${prefix}.pseudobulk_groups.tsv
    """
}

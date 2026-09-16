//
// Everything under the host's own clade, listed for removal.
//
// --drop_host_taxon removes one taxid, and that is not where host leakage
// lands. A host read that misses the host genome is assigned to the nearest
// relative the database does hold, so the leakage arrives wearing a sibling's
// name: on the CSI-Microbes 10x cohort, with Homo sapiens already dropped,
// Pan, Pongo, Gorilla, Macaca and Tupaia were still 84% of everything the
// pipeline called microbial in the uninfected library.
//
// Like the other evidence filters it removes nothing itself - the taxids go to
// ABUNDANCE_FILTER, so one step owns every removal. Unlike them it reaches the
// single-cell matrix as well, because it travels on the shared drop list rather
// than on FILTER_*'s --drop-taxid; --drop_host_taxon never did, which is why
// Homo sapiens rows outlived it in cell_taxa.
//
// The lineage is read from the run's OWN reports rather than from a taxonomy
// dump, so the ranks are the ones this database uses and there is nothing to
// keep in step with it.
//
process HOSTCLADE_EXPAND {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')
    val taxid
    val rank

    output:
    tuple val(meta), path("*.host_clade_drop.txt")    , emit: drop_list
    tuple val(meta), path("*.host_clade_evidence.tsv"), emit: evidence
    tuple val(meta), path("*_mqc.tsv")                , emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    host_clade.py \\
        --reports reports/* \\
        --taxid ${taxid} \\
        --rank ${rank} \\
        --drop-list ${prefix}.host_clade_drop.txt \\
        --evidence ${prefix}.host_clade_evidence.tsv \\
        --mqc ${prefix}.host_clade_mqc.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.host_clade_drop.txt ${prefix}.host_clade_evidence.tsv
    """
}

// Sum the per-sample host-k-mer tables and decide which taxa are host leakage.
//
// Cohort-level for the same reason the minimizer filter is: whether a taxon's
// reads are mostly host is a property of the taxon, and a single shallow
// library is a poor place to settle it. The taxids it condemns go to
// ABUNDANCE_FILTER alongside the minimizer list rather than being removed here,
// so every removal from the combined tables happens in one place.
process HOSTKMER_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(tables, stageAs: 'tables/*')
    val host_taxid
    val max_host_fraction
    val min_reads

    output:
    tuple val(meta), path("*.host_kmer_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.host_kmer_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    host_kmer_filter.py \\
        ${args} \\
        --merge \\
        --reads tables/* \\
        --host-taxid ${host_taxid} \\
        --max-host-fraction ${max_host_fraction} \\
        --min-reads ${min_reads} \\
        --output ${prefix}.host_kmer_evidence.tsv \\
        --drop-list ${prefix}.host_kmer_drop.txt \\
        --mqc ${prefix}_host_kmer_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.host_kmer_evidence.tsv ${prefix}.host_kmer_drop.txt ${prefix}_host_kmer_mqc.tsv
    """
}

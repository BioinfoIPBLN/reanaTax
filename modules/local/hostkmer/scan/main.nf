// Per-sample pass over Kraken2's read-level output, counting reads that carry
// host k-mers whatever they were classified as.
//
// One task per sample because the input is one line per read - gigabytes for a
// deep library - and streaming it is the expensive part. What comes out is a
// few hundred rows, which HOSTKMER_MERGE then sums across the cohort.
//
// See bin/host_kmer_filter.py for why the test is at k-mer rather than read
// assignment level, and for the requirement that the host be in the database.
process HOSTKMER_SCAN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(classifiedreads)
    val host_taxid

    output:
    tuple val(meta), path("*.host_kmer.tsv"), emit: table
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    host_kmer_filter.py \\
        ${args} \\
        --reads ${classifiedreads} \\
        --sample ${prefix} \\
        --host-taxid ${host_taxid} \\
        --output ${prefix}.host_kmer.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.host_kmer.tsv
    """
}

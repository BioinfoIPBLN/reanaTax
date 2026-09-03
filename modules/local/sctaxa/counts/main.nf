// The cell-by-taxon matrix, from Kraken2's read-level output.
//
// One task per sample: the input is one line per read. The barcode rides in the
// read name, put there by STARSOLO_UNMAPPED, because Kraken2 truncates a read
// name at the first whitespace and would otherwise drop it.
//
// The read-level filters are SAHMI's (host k-mers, min_frac lineage coherence,
// homopolymers); the UMI deduplication and the use of STARsolo's corrected
// barcode are not - SAHMI extracts the UMI and never uses it, and reads the
// barcode positionally with no whitelist. See bin/sc_taxa_counts.py.
process SCTAXA_COUNTS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(classifiedreads), path(report), path(fastq)
    val host_taxid
    val min_frac
    val max_homopolymer
    val ranks
    val umi_dedup
    val min_umis

    output:
    tuple val(meta), path("*.cell_taxa.tsv"), emit: counts
    tuple val(meta), path("*.cell_taxa_summary.tsv"), emit: summary
    tuple val(meta), path("*.cell_taxa_sweep.tsv"), emit: sweep
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def host_arg = host_taxid ? "--host-taxid ${host_taxid}" : ''
    // The FASTQ is only read when a homopolymer limit is set: the sequences are
    // not in Kraken2's output, and streaming them for nothing would double the
    // I/O of the most expensive step in this branch.
    def fastq_arg = max_homopolymer.toInteger() > 0 && fastq ? "--fastq ${fastq} --max-homopolymer ${max_homopolymer}" : ''
    def dedup_arg = umi_dedup ? '' : '--no-umi-dedup'
    """
    sc_taxa_counts.py \\
        ${args} \\
        --reads ${classifiedreads} \\
        --report ${report} \\
        --sample ${prefix} \\
        ${host_arg} \\
        --min-frac ${min_frac} \\
        --ranks '${ranks}' \\
        ${fastq_arg} \\
        ${dedup_arg} \\
        --min-umis ${min_umis} \\
        --sweep ${prefix}.cell_taxa_sweep.tsv \\
        --output ${prefix}.cell_taxa.tsv \\
        --summary ${prefix}.cell_taxa_summary.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cell_taxa.tsv ${prefix}.cell_taxa_summary.tsv ${prefix}.cell_taxa_sweep.tsv
    """
}

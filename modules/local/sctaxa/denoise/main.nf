// Does a taxon's evidence scale across CELLS?
//
// The same question --minimizer_correlation asks across samples, at a much
// finer grain: a study has a handful of libraries but thousands of barcodes, so
// this correlation has power where the sample-level one needs at least five
// samples to reach significance at all. It is SAHMI's single-cell contribution.
//
// Per sample, because barcodes only mean anything within the library that
// produced them. See bin/sc_kmer_denoise.py for why the distinct k-mers have to
// be recovered from the reads rather than read out of the report, and for the
// saturated-distinct case the correlation alone cannot express.
process SCTAXA_DENOISE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(classifiedreads), path(fastq)
    val host_taxid
    val kmer_len
    val min_barcodes
    val correlation_p
    val adjust
    val max_barcodes_per_taxon

    output:
    tuple val(meta), path("*.sc_kmer_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.sc_kmer_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def host_arg = host_taxid ? "--host-taxid ${host_taxid}" : ''
    """
    sc_kmer_denoise.py \\
        ${args} \\
        --reads ${classifiedreads} \\
        --fastq ${fastq} \\
        --sample ${prefix} \\
        ${host_arg} \\
        --kmer-len ${kmer_len} \\
        --min-barcodes ${min_barcodes} \\
        --correlation-p ${correlation_p} \\
        --adjust ${adjust} \\
        --max-barcodes-per-taxon ${max_barcodes_per_taxon} \\
        --evidence ${prefix}.sc_kmer_evidence.tsv \\
        --drop-list ${prefix}.sc_kmer_drop.txt \\
        --mqc ${prefix}_sc_kmer_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sc_kmer_evidence.tsv ${prefix}.sc_kmer_drop.txt ${prefix}_sc_kmer_mqc.tsv
    """
}

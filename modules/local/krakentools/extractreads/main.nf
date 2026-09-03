//
// Pull the reads Kraken2 assigned to one taxon back out of the non-host FASTQ.
//
// This is what turns a classification into an experiment. A Bracken row says
// "n reads looked like Plasmodium"; these are those reads, and once they are
// aligned to that organism's own genome they can be counted per gene and tested
// between groups like any other RNA-seq.
//
// Two flags are not optional. --fastq-output, because KrakenTools writes FASTA
// by default and everything downstream expects FASTQ. And --max, because it
// silently caps at 100 million reads if left alone; it is exposed as
// --target_max_reads so a truncation is a decision rather than a surprise.
//
process KRAKENTOOLS_EXTRACTREADS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakentools:1.2.1--pyh7e72e81_0'
        : 'quay.io/biocontainers/krakentools:1.2.1--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(reads), path(classifiedreads), path(report)
    val taxid

    output:
    tuple val(meta), path("*.taxon_*.fastq.gz"), emit: reads
    tuple val(meta), path("*.extract.log")     , emit: log
    tuple val("${task.process}"), val('krakentools'), eval("extract_kraken_reads.py --version 2>&1 | sed -n '1s/.*[[:space:]]//p'"), emit: versions_krakentools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --compression_level, for the gzip stream below.
    def level = task.ext.compression != null ? "-${task.ext.compression} " : ''
    def label = taxid.toString().replaceAll(/[^0-9]+/, '_')
    def seq_in = meta.single_end ? "-s ${reads}" : "-s1 ${reads[0]} -s2 ${reads[1]}"
    def seq_out = meta.single_end
        ? "-o ${prefix}.taxon_${label}.fastq"
        : "-o ${prefix}.taxon_${label}_1.fastq -o2 ${prefix}.taxon_${label}_2.fastq"
    """
    extract_kraken_reads.py \\
        -k ${classifiedreads} \\
        -r ${report} \\
        -t ${taxid.toString().replaceAll(',', ' ')} \\
        ${seq_in} \\
        ${seq_out} \\
        --fastq-output \\
        ${args} 2>&1 | tee ${prefix}.extract.log

    gzip ${level}-f ${prefix}.taxon_${label}*.fastq
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def label = taxid.toString().replaceAll(/[^0-9]+/, '_')
    def touches = meta.single_end
        ? "echo '' | gzip > ${prefix}.taxon_${label}.fastq.gz"
        : "echo '' | gzip > ${prefix}.taxon_${label}_1.fastq.gz\n    echo '' | gzip > ${prefix}.taxon_${label}_2.fastq.gz"
    """
    ${touches}
    touch ${prefix}.extract.log
    """
}

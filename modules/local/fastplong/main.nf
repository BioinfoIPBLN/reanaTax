//
// Quality and adapter trimming for long-read libraries: fastp's long-read
// counterpart, from the same authors, with the same report layout. The JSON
// keeps fastp's before_filtering/after_filtering summary, so read accounting
// reads it unchanged; only the mean length is `read_mean_length` rather than
// fastp's per-mate `read1_mean_length`. Adapters are detected at both ends
// unless --adapter_fasta lists them.
//
process FASTPLONG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/fastplong:0.7.0--h43da1c4_0'
        : 'quay.io/biocontainers/fastplong:0.7.0--h43da1c4_0'}"

    input:
    tuple val(meta), path(reads), path(adapter_fasta)
    val save_trimmed_fail

    output:
    tuple val(meta), path("*.fastplong.fastq.gz"), emit: reads
    tuple val(meta), path("*.fastplong.json"), emit: json
    tuple val(meta), path("*.fastplong.html"), emit: html
    tuple val(meta), path("*.fastplong.log"), emit: log
    tuple val(meta), path("*.fail.fastq.gz"), optional: true, emit: reads_fail
    tuple val("${task.process}"), val('fastplong'), eval("fastplong --version 2>&1 | sed -n 's/^fastplong //p'"), emit: versions_fastplong, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def adapters = adapter_fasta ? "--adapter_fasta ${adapter_fasta}" : ''
    def failed = save_trimmed_fail ? "--failed_out ${prefix}.fail.fastq.gz" : ''
    """
    fastplong \\
        --in ${reads} \\
        --out ${prefix}.fastplong.fastq.gz \\
        --json ${prefix}.fastplong.json \\
        --html ${prefix}.fastplong.html \\
        --thread ${task.cpus} \\
        ${adapters} \\
        ${failed} \\
        ${args} \\
        2> >(tee ${prefix}.fastplong.log >&2)
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.fastplong.fastq.gz
    touch ${prefix}.fastplong.json ${prefix}.fastplong.html ${prefix}.fastplong.log
    """
}

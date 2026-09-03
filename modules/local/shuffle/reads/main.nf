// A copy of every library with the sequence destroyed and the composition kept.
//
// The point of the copy is that it is classified by the same database with the
// same settings as the real reads, so whatever the classifier reports from it
// is what that database produces from base composition alone. See
// bin/shuffle_reads.py for what each shuffle preserves and why the
// dinucleotide-preserving one is the default.
//
// Runs on the NON-HOST reads, not the trimmed ones: the control has to answer
// the question the classification actually asked, and the classification was
// asked of the reads that survived host depletion.
process SHUFFLE_READS {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reads)
    val method
    val seed
    val max_reads

    output:
    tuple val(meta), path("*_shuffled{,_1,_2}.fastq.gz"), emit: reads
    tuple val(meta), path("*.shuffle_stats.tsv"), emit: stats
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shuffle_reads.py \\
        ${args} \\
        --reads ${reads} \\
        --prefix ${prefix} \\
        --method ${method} \\
        --seed ${seed} \\
        --max-reads ${max_reads} \\
        --stats ${prefix}.shuffle_stats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}_shuffled.fastq.gz
    touch ${prefix}.shuffle_stats.tsv
    """
}

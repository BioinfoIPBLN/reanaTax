//
// The same assembly as MEGAHIT, traded the other way: better contiguity for
// considerably more time and memory.
//
// Two constraints come from metaSPAdes itself rather than from this pipeline.
// It accepts exactly one paired-end library, so a pool has to be concatenated
// first - which MEGAHIT does not need, and which is the reason MEGAHIT is the
// default. And it refuses single-end input outright, so that is caught here
// with a sentence instead of a SPAdes traceback several minutes in.
//
process SPADES_META {
    tag "${meta.id}"
    label 'process_high'
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/spades:4.3.0--hde4eca7_1'
        : 'quay.io/biocontainers/spades:4.3.0--hde4eca7_1'}"

    input:
    tuple val(meta), path(reads1, stageAs: 'r1/*'), path(reads2, stageAs: 'r2/*')

    output:
    tuple val(meta), path("*.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("*.spades.log")   , emit: log
    tuple val("${task.process}"), val('spades'), eval("spades.py --version 2>&1 | sed 's/^SPAdes genome assembler v//'"), emit: versions_spades, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def level = task.ext.compression != null ? "-${task.ext.compression} " : ''
    if (!reads2) {
        error("SPADES_META: metaSPAdes needs paired-end reads. Use --assembly_assembler megahit for this cohort, which assembles single-end input.")
    }
    def left = (reads1 instanceof List ? reads1 : [reads1]).join(' ')
    def right = (reads2 instanceof List ? reads2 : [reads2]).join(' ')
    """
    cat ${left} > pooled_1.fastq.gz
    cat ${right} > pooled_2.fastq.gz

    spades.py \\
        --meta \\
        -1 pooled_1.fastq.gz \\
        -2 pooled_2.fastq.gz \\
        -o assembly \\
        -t ${task.cpus} \\
        -m ${task.memory.toGiga()} \\
        ${args} 2>&1 | tee ${prefix}.spades.log

    gzip ${level}-c assembly/contigs.fasta > ${prefix}.contigs.fa.gz
    rm -f pooled_1.fastq.gz pooled_2.fastq.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.contigs.fa.gz
    touch ${prefix}.spades.log
    """
}

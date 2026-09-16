//
// De novo assembly of the non-host fraction.
//
// Not an RNA assembler, on purpose, even though the input is usually RNA-seq.
// What comes out of host depletion is overwhelmingly bacterial and archaeal:
// unspliced, operonic, with nothing for an isoform model to reconstruct. What
// RNA-seq does impose is savage coverage skew - rRNA and a handful of
// transcripts swamp the rest - and that is a coverage problem, which is exactly
// what a multi-k metagenome assembler is built for. Trinity, which
// nf-rnaSeqMetagen uses here, is solving the wrong half of the problem.
//
// The exception is a eukaryotic microbe or a spliced virus, where transcript
// structure is real. That case is --target_taxid: you already know the organism,
// and TARGETED_TAXON aligns to its actual genome rather than guessing at one.
//
// Libraries arrive already grouped. MEGAHIT takes comma-separated lists, so a
// pool is assembled without concatenating anything to disk first.
//
process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/megahit:1.2.9--haf24da9_8'
        : 'quay.io/biocontainers/megahit:1.2.9--haf24da9_8'}"

    input:
    tuple val(meta), path(reads1, stageAs: 'r1/*'), path(reads2, stageAs: 'r2/*')

    output:
    tuple val(meta), path("*.contigs.fa.gz"), emit: contigs
    tuple val(meta), path("*.megahit.log")  , emit: log
    tuple val("${task.process}"), val('megahit'), eval("megahit --version | sed 's/^MEGAHIT v//'"), emit: versions_megahit, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def level = task.ext.compression != null ? "-${task.ext.compression} " : ''
    def left = (reads1 instanceof List ? reads1 : [reads1]).join(',')
    def right = reads2 ? (reads2 instanceof List ? reads2 : [reads2]).join(',') : ''
    def input_arg = right ? "-1 ${left} -2 ${right}" : "-r ${left}"
    """
    megahit \\
        ${input_arg} \\
        -t ${task.cpus} \\
        -m ${task.memory.toBytes()} \\
        --out-dir assembly \\
        --out-prefix ${prefix} \\
        ${args} 2>&1 | tee ${prefix}.megahit.log

    gzip ${level}-c assembly/${prefix}.contigs.fa > ${prefix}.contigs.fa.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | gzip > ${prefix}.contigs.fa.gz
    touch ${prefix}.megahit.log
    """
}

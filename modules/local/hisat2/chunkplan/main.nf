//
// How many read pairs a library holds, so it can be aligned in chunks.
//
// HISAT2 takes `-s/--skip` and `-u/--upto` in units of reads-or-pairs, which
// means a library can be split into slices without ever writing the slices to
// disk - each chunk task opens the same FASTQ and aligns its own window. What
// that needs, and the only thing it needs, is the total up front.
//
process HISAT2_CHUNKPLAN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.chunkplan.txt"), emit: plan

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def first = reads instanceof List ? reads[0] : reads
    """
    # Without pipefail a failed decompression still leaves wc reporting 0, and a
    # count of zero would be read downstream as "one empty chunk" rather than as
    # the error it is.
    set -o pipefail
    LINES=\$(zcat -f ${first} | wc -l)
    echo \$(( LINES / 4 )) > ${prefix}.chunkplan.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo 4 > ${prefix}.chunkplan.txt
    """
}

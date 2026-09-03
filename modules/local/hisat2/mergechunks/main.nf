//
// Join the per-chunk outputs of a chunked HISAT2 run back into the three files
// a single alignment would have produced, under the same names, so nothing
// downstream can tell the difference.
//
// The chunk files are staged into subdirectories because the outputs here are
// declared as globs: a chunk BAM sitting beside the merged one would match
// `*.bam` and be emitted as a second host BAM.
//
// Chunk names are zero-padded (`part001`, `part002`, ...), so shell glob
// expansion is already in chunk order and the concatenations are deterministic
// without an explicit sort.
//
process HISAT2_MERGECHUNKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(bams, stageAs: 'chunk_bam/*'), path(summaries, stageAs: 'chunk_log/*'), path(fastqs, stageAs: 'chunk_fastq/*')
    val save_unaligned

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.hisat2.summary.log")    , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_fastq = ''
    if (save_unaligned) {
        // gzip members concatenate: `cat a.gz b.gz` is itself a valid gzip
        // stream, so there is no decompression pass here.
        merge_fastq = meta.single_end
            ? "cat chunk_fastq/*.unmapped.fastq.gz > ${prefix}.unmapped.fastq.gz"
            : """cat chunk_fastq/*.unmapped_1.fastq.gz > ${prefix}.unmapped_1.fastq.gz
    cat chunk_fastq/*.unmapped_2.fastq.gz > ${prefix}.unmapped_2.fastq.gz"""
    }
    """
    set -o pipefail

    # `samtools cat` takes the header of the first file for the whole output, so
    # every chunk has to carry the same @RG - which is why the chunk tasks are
    # given an explicit --rg-id in conf/modules.config instead of letting the
    # module derive one from the (per-chunk) prefix. It concatenates without
    # decompressing or re-sorting; the sort happens downstream as it always did.
    samtools cat -o ${prefix}.bam chunk_bam/*.bam

    merge_hisat2_summary.py chunk_log/*.hisat2.summary.log --output ${prefix}.hisat2.summary.log

    ${merge_fastq}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_fastq = save_unaligned
        ? (meta.single_end
            ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz"
            : "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz\n    echo '' | gzip > ${prefix}.unmapped_2.fastq.gz")
        : ''
    """
    touch ${prefix}.bam
    touch ${prefix}.hisat2.summary.log
    ${merge_fastq}
    """
}

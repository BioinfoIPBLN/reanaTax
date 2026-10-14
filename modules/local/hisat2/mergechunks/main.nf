//
// Join the per-chunk outputs of a chunked HISAT2 run back into the three files
// a single alignment would have produced, under the same names, so nothing
// downstream can tell the difference.
//
// The chunk files are staged into subdirectories because the outputs here are
// declared as globs: a chunk BAM sitting beside the merged one would match
// `*.bam` and be emitted as a second host BAM.
//
// Chunk names are zero-padded (`part0001`, `part0002`, ...), so shell glob
// expansion is already in chunk order and the concatenations are deterministic
// without an explicit sort.
//
// With `save_bam` false there are no chunk BAMs to join - the split aligner
// writes none then - and the sample gets no host BAM, as an unchunked run would
// not. The nf-core aligner writes one regardless; those are handed in anyway,
// so --cleanup_intermediates can still delete them.
//
process HISAT2_MERGECHUNKS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(bams, stageAs: 'chunk_bam/*'), path(summaries, stageAs: 'chunk_log/*'), path(fastqs, stageAs: 'chunk_fastq/*'), path(chunk_reads, stageAs: 'chunk_reads/*')
    val save_bam
    val save_unaligned
    val cleanup_intermediates
    val work_dir

    output:
    tuple val(meta), path("*.bam")    , optional:true, emit: bam
    tuple val(meta), path("*.summary.log")           , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --cleanup_intermediates. The chunk files exist only for this sample's
    // alignment - the slices of the library the chunk tasks read, and the chunk
    // BAMs concatenated here - and for a chunked run they are the largest things
    // on disk: a second copy of the library and a full copy of every host
    // alignment, in pieces, kept for the lifetime of the run.
    //
    // Deleted from this task rather than from a cleanup process of its own,
    // because the ordering that matters - "every chunk aligned and the merge
    // succeeded" - is then a property of the script rather than of a channel
    // join. `set -e` is in effect, so a failed merge never reaches these lines.
    //
    // Two guards. The staged names are symlinks, so it is the link TARGET that
    // has to go; and that target is only removed when it is inside the work
    // directory, so a staging mode that hard-links or copies real inputs -
    // or a future caller that hands this module something it does not own -
    // cannot lose data.
    //
    // The work directory arrives as a String input rather than as
    // `workflow.workDir` in the script. Nextflow hashes every variable a script
    // mentions, whether or not the branch that uses it runs, and a Path is
    // hashed by its size and modification time. The work directory's change
    // whenever a run creates task folders, so this task - and everything
    // downstream of it - missed the cache on every resume.
    def cleanup = cleanup_intermediates
        ? """
    shopt -s nullglob
    work_root=\$(readlink -f "${work_dir}")
    for chunk in chunk_bam/*.bam chunk_reads/*.fastq.gz; do
        target=\$(readlink -f "\$chunk" || true)
        case "\$target" in
            "\$work_root"/*) rm -f "\$target" ;;
        esac
    done
    """
        : ''
    def merge_fastq = ''
    if (save_unaligned) {
        // gzip members concatenate: `cat a.gz b.gz` is itself a valid gzip
        // stream, so there is no decompression pass here.
        merge_fastq = meta.single_end
            ? "cat chunk_fastq/*.unmapped.fastq.gz > ${prefix}.unmapped.fastq.gz"
            : """cat chunk_fastq/*.unmapped_1.fastq.gz > ${prefix}.unmapped_1.fastq.gz
    cat chunk_fastq/*.unmapped_2.fastq.gz > ${prefix}.unmapped_2.fastq.gz"""
    }
    // `samtools cat` takes the header of the first file for the whole output, so
    // every chunk has to carry the same @RG - which is why the chunk tasks are
    // given an explicit --rg-id in conf/modules.config instead of letting the
    // module derive one from the (per-chunk) prefix. It concatenates without
    // decompressing or re-sorting; the sort happens downstream as it always did.
    def merge_bam = save_bam ? "samtools cat -o ${prefix}.bam chunk_bam/*.bam" : ''
    """
    set -o pipefail

    ${merge_bam}

    # The chunks were aligned with HISAT2 or Bowtie2, which write the same
    # summary format; the merged log keeps the aligner's name.
    aligner=hisat2
    if compgen -G 'chunk_log/*.bowtie2.summary.log' > /dev/null; then
        aligner=bowtie2
    fi
    merge_hisat2_summary.py chunk_log/*.summary.log --output ${prefix}.\${aligner}.summary.log

    ${merge_fastq}
    ${cleanup}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def merge_fastq = save_unaligned
        ? (meta.single_end
            ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz"
            : "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz\n    echo '' | gzip > ${prefix}.unmapped_2.fastq.gz")
        : ''
    def merge_bam = save_bam ? "touch ${prefix}.bam" : ''
    """
    ${merge_bam}
    touch ${prefix}.hisat2.summary.log
    ${merge_fastq}
    """
}

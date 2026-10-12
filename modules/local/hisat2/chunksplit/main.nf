//
// Cut a library into files of --hisat2_chunk_size reads-or-pairs, so it can be
// aligned in chunks.
//
// One streaming pass per mate: decompress, cut every 4 x size lines, and
// recompress each slice as it is written. The chunks used to be offsets into
// the whole FASTQ instead (HISAT2's `-s/--skip` and `-u/--upto`), which writes
// nothing - but every chunk then decompressed and parsed the library up to its
// own window, so the reading grew with the square of the chunk count. On four
// libraries of ~540 M pairs, 27 chunks each, the last chunks of a library took
// about an hour longer than the first. One pass here costs minutes per library,
// and a second copy of it on disk until the merge, which --cleanup_intermediates
// then removes.
//
// A FASTQ record is four lines, as fastp and the archives write it. Both mates
// are cut at the same record counts, so chunk N of R1 and chunk N of R2 hold the
// same pairs. The task stops if the mates come out in different numbers of
// chunks, and HISAT2 itself refuses a chunk whose mates hold different numbers
// of reads.
//
process HISAT2_CHUNKSPLIT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads)
    val chunk_size

    output:
    tuple val(meta), path("chunks/*.fastq.gz"), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mates = reads instanceof List ? reads : [reads]
    def lines = (chunk_size as long) * 4
    // The mates are cut at the same time and share the task's cores: zcat and
    // split take one each per mate, and bgzip - the one compressor in the image
    // that can use more than one - the rest.
    def threads = Math.max(1, (task.cpus as int).intdiv(mates.size()) - 2)
    // --compression_level, as everywhere else. Unset leaves bgzip at its own.
    def level = task.ext.compression != null ? "-l ${task.ext.compression} " : ''
    def cut = mates.withIndex().collect { fq, i -> "cut_mate ${fq} ${i + 1} & pids+=(\$!)" }.join('\n    ')
    def empty = mates.withIndex().collect { _fq, i -> "bgzip -c < /dev/null > chunks/${prefix}.part0001_${i + 1}.fastq.gz" }.join('\n        ')
    def mate_check = mates.size() > 1
        ? """
    if [ "\$(find chunks -name '*_1.fastq.gz' | wc -l)" -ne "\$(find chunks -name '*_2.fastq.gz' | wc -l)" ]; then
        echo "${meta.id}: R1 and R2 were cut into different numbers of chunks, so the mates do not hold the same number of reads." >&2
        exit 1
    fi"""
        : ''
    """
    mkdir -p chunks

    # split runs the filter once per slice, with \$FILE set to that slice's name.
    cut_mate() {
        zcat -f "\$1" \\
            | split -l ${lines} -a 4 --numeric-suffixes=1 \\
                --additional-suffix="_\$2.fastq.gz" \\
                --filter='bgzip -c -@ ${threads} ${level}> "\$FILE"' \\
                - chunks/${prefix}.part
    }

    pids=()
    ${cut}
    for pid in "\${pids[@]}"; do
        wait "\$pid"
    done

    # split writes nothing for an empty library. One empty chunk instead, so the
    # sample still produces the summary and the empty outputs downstream expects.
    if [ -z "\$(ls chunks)" ]; then
        ${empty}
    fi
    ${mate_check}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mates = reads instanceof List ? reads : [reads]
    def empty = mates.withIndex().collect { _fq, i -> "echo '' | gzip > chunks/${prefix}.part0001_${i + 1}.fastq.gz" }.join('\n    ')
    """
    mkdir -p chunks
    ${empty}
    """
}

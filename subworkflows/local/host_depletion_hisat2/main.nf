//
// Split reads into a host fraction (sorted, indexed BAM) and a non-host
// fraction (FASTQ) by aligning against the host genome with HISAT2.
//

include { HISAT2_ALIGN             } from '../../../modules/nf-core/hisat2/align/main'
include { HISAT2_ALIGN_SPLIT       } from '../../../modules/local/hisat2/alignsplit/main'
include { HISAT2_ALIGN as HISAT2_ALIGN_CHUNK             } from '../../../modules/nf-core/hisat2/align/main'
include { HISAT2_ALIGN_SPLIT as HISAT2_ALIGN_SPLIT_CHUNK } from '../../../modules/local/hisat2/alignsplit/main'
include { HISAT2_CHUNKPLAN         } from '../../../modules/local/hisat2/chunkplan/main'
include { HISAT2_MERGECHUNKS       } from '../../../modules/local/hisat2/mergechunks/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { QUALIMAP_BAMQC           } from '../../../modules/nf-core/qualimap/bamqc/main'
include { SUBREAD_FEATURECOUNTS    } from '../../../modules/nf-core/subread/featurecounts/main'

workflow HOST_DEPLETION_HISAT2 {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    ch_index // channel: [ val(meta), path(hisat2_index_dir) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    save_host_bam // boolean: sort, index and stat the aligned (host) reads
    skip_qualimap // boolean
    qualimap_gff // string: optional GFF/GTF to add feature-level stats, or null
    quantify_gtf // string: GTF to count host reads against, or null to skip
    require_both_mates_unmapped // boolean: release a pair only if NEITHER mate aligned
    chunk_size // integer: align this many read pairs per task, or 0 to align the library in one

    main:

    def ch_multiqc_files = channel.empty()

    //
    // MODULE: Align against the host. `save_unaligned` is what produces the
    // non-host FASTQs this pipeline is really after.
    //
    // The default module takes its non-host reads from `--un-conc-gz`, i.e.
    // every pair that failed to align CONCORDANTLY, so a pair with one host
    // mate leaves whole. HISAT2_ALIGN_SPLIT drops `--no-mixed --no-discordant`
    // and splits on explicit SAM flags instead, releasing a pair only when both
    // mates are unmapped. Same outputs, same names.
    def ch_aligned = channel.empty()
    def ch_summary = channel.empty()
    def ch_fastq = channel.empty()

    //
    // Chunked alignment. HISAT2's memory grows with how much it has aligned in
    // one process - on libraries of a few hundred million pairs with
    // --very-sensitive it can reach the hundreds of GB - and the only thing
    // that resets it is a fresh process. Chunking gives one process per slice.
    //
    // The slices are offsets, not files: `-s/--skip` and `-u/--upto` count
    // reads-or-pairs, so each chunk task opens the same FASTQ and aligns its
    // own window. Nothing is split, recompressed or written to scratch, and
    // because the chunks are independent tasks they also run in parallel -
    // which turns the longest step in the pipeline from one very long task into
    // as many short ones as the executor will schedule.
    //
    if (chunk_size && chunk_size > 0) {
        HISAT2_CHUNKPLAN(ch_reads)

        def ch_chunks = HISAT2_CHUNKPLAN.out.plan
            .map { meta, plan -> [meta, plan.text.trim() as long] }
            .join(ch_reads)
            .flatMap { meta, total, reads ->
                // Ceiling division, and at least one chunk even for an empty
                // library, so a sample with no surviving reads still produces the
                // summary and the empty outputs the rest of the pipeline expects.
                def size = chunk_size as long
                def n_chunks = total > 0 ? (total + size - 1).intdiv(size) : 1
                (0..<n_chunks).collect { index ->
                    [
                        meta + [
                            chunk: String.format('part%03d', index + 1),
                            chunk_skip: index * size,
                            chunk_upto: size,
                        ],
                        reads,
                    ]
                }
            }

        def ch_chunk_bam = channel.empty()
        def ch_chunk_summary = channel.empty()
        def ch_chunk_fastq = channel.empty()

        if (require_both_mates_unmapped) {
            HISAT2_ALIGN_SPLIT_CHUNK(ch_chunks, ch_index.first(), channel.value([[:], []]), true)
            ch_chunk_bam = HISAT2_ALIGN_SPLIT_CHUNK.out.bam
            ch_chunk_summary = HISAT2_ALIGN_SPLIT_CHUNK.out.summary
            ch_chunk_fastq = HISAT2_ALIGN_SPLIT_CHUNK.out.fastq
        }
        else {
            HISAT2_ALIGN_CHUNK(ch_chunks, ch_index.first(), channel.value([[:], []]), true)
            ch_chunk_bam = HISAT2_ALIGN_CHUNK.out.bam
            ch_chunk_summary = HISAT2_ALIGN_CHUNK.out.summary
            ch_chunk_fastq = HISAT2_ALIGN_CHUNK.out.fastq
        }

        // Strip the chunk keys to recover the sample's own meta, then regroup on
        // it. Each list is sorted by name so a resumed run hands the merge task
        // its inputs in the same order as the original one: groupTuple makes no
        // promise about order, and an unsorted list would change the task hash
        // from run to run and defeat -resume.
        def chunk_keys = ['chunk', 'chunk_skip', 'chunk_upto']

        def ch_merged_bam = ch_chunk_bam
            .map { meta, files -> [meta.findAll { key, _value -> !(key in chunk_keys) }, files] }
            .groupTuple()
            .map { meta, files -> [meta, files.flatten().sort { entry -> entry.name }] }

        def ch_merged_summary = ch_chunk_summary
            .map { meta, files -> [meta.findAll { key, _value -> !(key in chunk_keys) }, files] }
            .groupTuple()
            .map { meta, files -> [meta, files.flatten().sort { entry -> entry.name }] }

        def ch_merged_fastq = ch_chunk_fastq
            .map { meta, files -> [meta.findAll { key, _value -> !(key in chunk_keys) }, files] }
            .groupTuple()
            .map { meta, files -> [meta, files.flatten().sort { entry -> entry.name }] }

        // `true`, like the unchunked calls above: the non-host FASTQs are the
        // pipeline's product, not an optional extra. params.save_unaligned only
        // decides whether they are also published.
        HISAT2_MERGECHUNKS(
            ch_merged_bam.join(ch_merged_summary).join(ch_merged_fastq),
            true,
        )

        ch_aligned = HISAT2_MERGECHUNKS.out.bam
        ch_summary = HISAT2_MERGECHUNKS.out.summary
        ch_fastq = HISAT2_MERGECHUNKS.out.fastq
    }
    else if (require_both_mates_unmapped) {
        HISAT2_ALIGN_SPLIT(ch_reads, ch_index.first(), channel.value([[:], []]), true)
        ch_aligned = HISAT2_ALIGN_SPLIT.out.bam
        ch_summary = HISAT2_ALIGN_SPLIT.out.summary
        ch_fastq = HISAT2_ALIGN_SPLIT.out.fastq
    }
    else {
        HISAT2_ALIGN(ch_reads, ch_index.first(), channel.value([[:], []]), true)
        ch_aligned = HISAT2_ALIGN.out.bam
        ch_summary = HISAT2_ALIGN.out.summary
        ch_fastq = HISAT2_ALIGN.out.fastq
    }

    ch_multiqc_files = ch_multiqc_files.mix(ch_summary.map { _meta, log -> log })

    def ch_bam = channel.empty()
    def ch_bai = channel.empty()
    def ch_qualimap = channel.empty()

    if (save_host_bam) {
        //
        // SUBWORKFLOW: Sort + index + stats, so the host BAM is directly usable
        // (and its mapping rate lands in MultiQC).
        //
        BAM_SORT_STATS_SAMTOOLS(
            ch_aligned,
            ch_fasta.map { meta, fasta -> [meta, fasta, []] }.first(),
        )
        ch_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        ch_bai = BAM_SORT_STATS_SAMTOOLS.out.index
        ch_multiqc_files = ch_multiqc_files
            .mix(BAM_SORT_STATS_SAMTOOLS.out.stats.map { _meta, stats -> stats })
            .mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat.map { _meta, flagstat -> flagstat })
            .mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.map { _meta, idxstats -> idxstats })

        //
        // MODULE: Qualimap BamQC on the sorted host BAM. samtools stats above
        // already give the read-level counts; Qualimap adds what they cannot -
        // coverage depth and its uniformity across the reference, duplication
        // rate, GC of the mapped reads, mapping-quality and insert-size
        // distributions - which is how you tell "little host in the library"
        // apart from "the host reference is wrong".
        //
        if (!skip_qualimap) {
            QUALIMAP_BAMQC(
                ch_bam,
                qualimap_gff ? file(qualimap_gff, checkIfExists: true) : [],
            )
            ch_qualimap = QUALIMAP_BAMQC.out.results
            ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQC.out.results.map { _meta, results -> results })
        }
    }

    //
    // MODULE: Host gene counts.
    //
    // The BAM is already exactly what featureCounts needs: the HISAT2 module
    // pipes through `samtools view -F 256`, so secondary alignments never
    // reached it and every aligned read appears once. Without that filter,
    // `--very-sensitive` (-k 50) would have counted each read up to 50 times.
    //
    // This is what turns the host half of the split from a QC by-product into
    // data - the expression matrix that the microbial profile can be correlated
    // against, from the very same library.
    //
    def ch_host_counts = channel.empty()
    def ch_host_counts_summary = channel.empty()

    if (save_host_bam && quantify_gtf) {
        SUBREAD_FEATURECOUNTS(
            ch_bam.map { meta, bam -> [meta, bam, file(quantify_gtf, checkIfExists: true)] }
        )
        ch_host_counts = SUBREAD_FEATURECOUNTS.out.counts
        ch_host_counts_summary = SUBREAD_FEATURECOUNTS.out.summary
        ch_multiqc_files = ch_multiqc_files.mix(SUBREAD_FEATURECOUNTS.out.summary.map { _meta, summary -> summary })
    }

    //
    // For paired-end input HISAT2 writes `<id>.unmapped_1.fastq.gz` and
    // `<id>.unmapped_2.fastq.gz`; sort so mate 1 always comes first.
    //
    def ch_unaligned = ch_fastq.map { meta, fastq ->
        def files = (fastq instanceof List ? fastq : [fastq]).sort { fq -> fq.name }
        [meta, files]
    }

    emit:
    reads = ch_unaligned // channel: [ val(meta), [ path(fastq) ] ]
    bam = ch_bam // channel: [ val(meta), path(bam) ]
    bai = ch_bai // channel: [ val(meta), path(bai) ]
    qualimap = ch_qualimap // channel: [ val(meta), path(results_dir) ]
    host_counts = ch_host_counts // channel: [ val(meta), path(featureCounts.tsv) ]
    host_counts_summary = ch_host_counts_summary // channel: [ val(meta), path(summary) ]
    summary = ch_summary // channel: [ val(meta), path(log) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

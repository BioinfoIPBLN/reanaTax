//
// Split reads into a host fraction (sorted, indexed BAM) and a non-host
// fraction (FASTQ) by aligning against the host genome: HISAT2 for RNA,
// Bowtie2 for DNA, minimap2 for long reads (--host_aligner; see hostAligner
// below).
//

include { HISAT2_ALIGN             } from '../../../modules/nf-core/hisat2/align/main'
include { HISAT2_ALIGN_SPLIT       } from '../../../modules/local/hisat2/alignsplit/main'
include { HISAT2_ALIGN as HISAT2_ALIGN_CHUNK             } from '../../../modules/nf-core/hisat2/align/main'
include { HISAT2_ALIGN_SPLIT as HISAT2_ALIGN_SPLIT_CHUNK } from '../../../modules/local/hisat2/alignsplit/main'
include { BOWTIE2_ALIGN_SPLIT      } from '../../../modules/local/bowtie2/alignsplit/main'
include { BOWTIE2_ALIGN_SPLIT as BOWTIE2_ALIGN_SPLIT_CHUNK } from '../../../modules/local/bowtie2/alignsplit/main'
include { MINIMAP2_ALIGN_SPLIT     } from '../../../modules/local/minimap2/alignsplit/main'
include { HISAT2_CHUNKSPLIT        } from '../../../modules/local/hisat2/chunksplit/main'
include { HISAT2_MERGECHUNKS       } from '../../../modules/local/hisat2/mergechunks/main'
include { SAMTOOLS_FAIDX           } from '../../../modules/local/samtools/faidx/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/local/bam_sort_stats_samtools/main'
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
    output_format // string: 'bam' or 'cram' for the sorted host alignment
    cleanup_intermediates // boolean: delete a sample's chunk files once its merge has succeeded
    ch_bowtie2_index // channel: [ val(meta), path(bowtie2_index_dir) ] - empty unless a library needs it
    aligner_mode // string: 'auto' (DNA to Bowtie2, the rest to HISAT2), 'hisat2' or 'bowtie2'
    ch_minimap2_index // channel: [ val(meta + preset), path(mmi) ] - empty unless a long-read library needs it
    minimap2_preset // string: the minimap2 preset for every long-read library, or null to choose per library

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
    // Long-read libraries go to minimap2 whole, each with the index built for
    // its preset. Chunking is for HISAT2's memory, which grows with what one
    // process has aligned; minimap2's does not.
    //
    def ch_by_aligner = ch_reads.branch { meta, _reads ->
        minimap2: hostAligner(meta, aligner_mode) == 'minimap2'
        short_reads: true
    }
    def ch_long_input = ch_by_aligner.minimap2
        .map { meta, reads -> [minimap2Preset(meta, minimap2_preset), meta, reads] }
        .combine(ch_minimap2_index.map { meta, index -> [meta.preset, meta, index] }, by: 0)
        .multiMap { _preset, meta, reads, meta_index, index ->
            reads: [meta, reads]
            index: [meta_index, index]
        }
    MINIMAP2_ALIGN_SPLIT(ch_long_input.reads, ch_long_input.index, true, save_host_bam)
    def ch_short_reads = ch_by_aligner.short_reads

    //
    // Chunked alignment. HISAT2's memory grows with how much it has aligned in
    // one process - on libraries of a few hundred million pairs with
    // --very-sensitive it can reach the hundreds of GB - and the only thing
    // that resets it is a fresh process. Chunking gives one process per slice.
    //
    // The slices are files, cut by HISAT2_CHUNKSPLIT in one pass over the
    // library. They used to be offsets into the whole FASTQ, which made every
    // chunk read the library up to its own window (see that module). The
    // chunks are independent tasks, so they also run in parallel - which turns
    // the longest step in the pipeline from one very long task into as many
    // short ones as the executor will schedule.
    //
    if (chunk_size && chunk_size > 0) {
        HISAT2_CHUNKSPLIT(ch_short_reads, chunk_size)

        // One list per sample, whether the glob matched one file or several.
        def ch_split = HISAT2_CHUNKSPLIT.out.reads
            .map { meta, files -> [meta, files instanceof List ? files : [files]] }

        // One task per chunk, named from its files (`<id>.part0001_1.fastq.gz`
        // is part0001), mate 1 first. A single-end chunk is handed over as the
        // file itself, the way the unchunked path hands over the library.
        // toSorted rather than sort, which would reorder the task's own output
        // list in place (see ch_unaligned below).
        def ch_chunks = ch_split
            .flatMap { meta, files ->
                files
                    .groupBy { fastq -> chunkName(fastq) }
                    .collect { chunk, mates ->
                        def sorted = mates.toSorted { fastq -> fastq.name }
                        [meta + [chunk: chunk], meta.single_end ? sorted[0] : sorted]
                    }
            }

        def ch_chunk_bam = channel.empty()
        def ch_chunk_summary = channel.empty()
        def ch_chunk_fastq = channel.empty()

        def ch_chunks_routed = ch_chunks.branch { meta, _reads ->
            bowtie2: hostAligner(meta, aligner_mode) == 'bowtie2'
            hisat2: true
        }

        if (require_both_mates_unmapped) {
            HISAT2_ALIGN_SPLIT_CHUNK(ch_chunks_routed.hisat2, ch_index.first(), channel.value([[:], []]), true, save_host_bam)
            ch_chunk_bam = HISAT2_ALIGN_SPLIT_CHUNK.out.bam
            ch_chunk_summary = HISAT2_ALIGN_SPLIT_CHUNK.out.summary
            ch_chunk_fastq = HISAT2_ALIGN_SPLIT_CHUNK.out.fastq
        }
        else {
            HISAT2_ALIGN_CHUNK(ch_chunks_routed.hisat2, ch_index.first(), channel.value([[:], []]), true)
            ch_chunk_bam = HISAT2_ALIGN_CHUNK.out.bam
            ch_chunk_summary = HISAT2_ALIGN_CHUNK.out.summary
            ch_chunk_fastq = HISAT2_ALIGN_CHUNK.out.fastq
        }

        // Bowtie2 always releases a pair only when neither mate aligned: there is
        // no --un-conc variant of it here.
        BOWTIE2_ALIGN_SPLIT_CHUNK(ch_chunks_routed.bowtie2, ch_bowtie2_index.first(), true, save_host_bam)
        ch_chunk_bam = ch_chunk_bam.mix(BOWTIE2_ALIGN_SPLIT_CHUNK.out.bam)
        ch_chunk_summary = ch_chunk_summary.mix(BOWTIE2_ALIGN_SPLIT_CHUNK.out.summary)
        ch_chunk_fastq = ch_chunk_fastq.mix(BOWTIE2_ALIGN_SPLIT_CHUNK.out.fastq)

        // Strip the chunk key to recover the sample's own meta, then regroup on
        // it. Each list is sorted by name so a resumed run hands the merge task
        // its inputs in the same order as the original one: groupTuple makes no
        // promise about order, and an unsorted list would change the task hash
        // from run to run and defeat -resume.
        //
        // Each sample is released as soon as its own chunks are aligned. A plain
        // groupTuple waits for the channel to close, i.e. for every chunk of
        // every sample, so no merge, no second host pass and no classification
        // could start until the whole cohort had aligned. groupKey carries the
        // sample's chunk count, from the same files that made the chunks; an
        // incomplete group is still emitted at the end (remainder), as before.
        def chunk_keys = ['chunk']
        def ch_chunk_count = ch_split.map { meta, files -> [meta, files.collect { fastq -> chunkName(fastq) }.unique().size()] }

        def ch_merged_summary = regroupChunks(ch_chunk_summary, ch_chunk_count, chunk_keys)
        def ch_merged_fastq = regroupChunks(ch_chunk_fastq, ch_chunk_count, chunk_keys)
        // The split aligners write chunk BAMs only when the host BAM is kept;
        // the nf-core HISAT2 one always does. Decided per library, since a run
        // can mix both aligners: a library with none coming is handed an empty
        // list, or the join below would drop it.
        def ch_bam_routes = ch_merged_summary.branch { meta, _summaries ->
            with_bam: save_host_bam || (!require_both_mates_unmapped && hostAligner(meta, aligner_mode) == 'hisat2')
            without: true
        }
        def ch_merged_bam = ch_bam_routes.with_bam
            .join(regroupChunks(ch_chunk_bam, ch_chunk_count, chunk_keys))
            .map { meta, _summaries, bams -> [meta, bams] }
            .mix(ch_bam_routes.without.map { meta, _summaries -> [meta, []] })

        // `true`, like the unchunked calls above: the non-host FASTQs are the
        // pipeline's product, not an optional extra. params.save_unaligned only
        // decides whether they are also published. The library's chunk files go
        // in as well, for --cleanup_intermediates to delete after the merge.
        HISAT2_MERGECHUNKS(
            ch_merged_bam.join(ch_merged_summary).join(ch_merged_fastq).join(ch_split),
            save_host_bam,
            true,
            cleanup_intermediates,
            workflow.workDir.toString(),
        )

        ch_aligned = HISAT2_MERGECHUNKS.out.bam
        ch_summary = HISAT2_MERGECHUNKS.out.summary
        ch_fastq = HISAT2_MERGECHUNKS.out.fastq
    }
    else {
        def ch_routed = ch_short_reads.branch { meta, _reads ->
            bowtie2: hostAligner(meta, aligner_mode) == 'bowtie2'
            hisat2: true
        }

        if (require_both_mates_unmapped) {
            HISAT2_ALIGN_SPLIT(ch_routed.hisat2, ch_index.first(), channel.value([[:], []]), true, save_host_bam)
            ch_aligned = HISAT2_ALIGN_SPLIT.out.bam
            ch_summary = HISAT2_ALIGN_SPLIT.out.summary
            ch_fastq = HISAT2_ALIGN_SPLIT.out.fastq
        }
        else {
            HISAT2_ALIGN(ch_routed.hisat2, ch_index.first(), channel.value([[:], []]), true)
            ch_aligned = HISAT2_ALIGN.out.bam
            ch_summary = HISAT2_ALIGN.out.summary
            ch_fastq = HISAT2_ALIGN.out.fastq
        }

        BOWTIE2_ALIGN_SPLIT(ch_routed.bowtie2, ch_bowtie2_index.first(), true, save_host_bam)
        ch_aligned = ch_aligned.mix(BOWTIE2_ALIGN_SPLIT.out.bam)
        ch_summary = ch_summary.mix(BOWTIE2_ALIGN_SPLIT.out.summary)
        ch_fastq = ch_fastq.mix(BOWTIE2_ALIGN_SPLIT.out.fastq)
    }

    ch_aligned = ch_aligned.mix(MINIMAP2_ALIGN_SPLIT.out.bam)
    ch_summary = ch_summary.mix(MINIMAP2_ALIGN_SPLIT.out.summary)
    ch_fastq = ch_fastq.mix(MINIMAP2_ALIGN_SPLIT.out.fastq)

    ch_multiqc_files = ch_multiqc_files.mix(ch_summary.map { _meta, log -> log })

    def ch_bam = channel.empty()
    def ch_bai = channel.empty()
    def ch_qualimap = channel.empty()

    if (save_host_bam) {
        //
        // SUBWORKFLOW: Sort + index + stats, so the host BAM is directly usable
        // (and its mapping rate lands in MultiQC).
        //
        // CRAM stores differences from a reference, so the writer, the reader
        // and samtools stats all need it indexed. Built once here rather than
        // left to htslib, which would otherwise rebuild it inside every task.
        // A BAM run needs none of it and pays for none of it.
        def ch_fasta_fai = ch_fasta.map { meta, fasta -> [meta, fasta, []] }.first()
        if (output_format == 'cram') {
            SAMTOOLS_FAIDX(ch_fasta)
            ch_fasta_fai = ch_fasta.join(SAMTOOLS_FAIDX.out.fai).first()
        }

        BAM_SORT_STATS_SAMTOOLS(
            ch_aligned,
            ch_fasta_fai,
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
            // Qualimap stops on a BAM with no mapped reads ("The BAM file is
            // empty or corrupt"), and a library with no host reads - a pure
            // culture, reads already depleted, the test profile's SARS-CoV-2
            // against a human genome - produces exactly that. There is nothing
            // in it to describe, so such BAMs are left out.
            def ch_bam_with_reads = ch_bam
                .join(BAM_SORT_STATS_SAMTOOLS.out.flagstat)
                .filter { _meta, _bam, flagstat ->
                    def mapped = flagstat.text.readLines().find { line -> line.contains(' mapped (') }
                    mapped && (mapped.tokenize(' ')[0] as long) > 0
                }
                .map { meta, bam, _flagstat -> [meta, bam] }
            QUALIMAP_BAMQC(
                ch_bam_with_reads,
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
    // toSorted, not sort: Groovy's List.sort sorts IN PLACE, and this list is
    // the task's own output object. Nextflow binds a task's outputs to its
    // channels before it publishes them, so the publisher can still be
    // iterating the list while this operator reorders it - a
    // ConcurrentModificationException that aborts the run ("Unexpected error
    // while finalizing task ... cause: null").
    //
    def ch_unaligned = ch_fastq.map { meta, fastq ->
        def files = (fastq instanceof List ? fastq : [fastq]).toSorted { fq -> fq.name }
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

// Which aligner depletes this library. A long-read library (meta.long_reads,
// set from --long_reads) always goes to minimap2: neither short-read aligner
// takes reads of tens of kilobases. HISAT2 is splice-aware, which RNA needs
// and DNA does not: it would split a genomic read across an intron that is not
// one. Under 'auto' a library known to be DNA (meta.molecule, from the archive's
// library strategy and source or the samplesheet's `molecule` column) goes to
// Bowtie2, and everything else - RNA, and libraries of unknown kind, since this
// is an RNA-seq pipeline first - to HISAT2. 'hisat2' or 'bowtie2' forces one.
def hostAligner(meta, mode) {
    if (meta.long_reads) {
        return 'minimap2'
    }
    if (mode && mode != 'auto') {
        return mode
    }
    return meta.molecule == 'dna' ? 'bowtie2' : 'hisat2'
}

// The minimap2 preset a long-read library is aligned with, unless
// --minimap2_preset names one for every library. Spliced for RNA and for
// libraries of unknown kind, as hostAligner treats them; genomic for DNA. Each
// in its platform's flavour: HiFi accuracy for PacBio, the ONT error model
// otherwise.
def minimap2Preset(meta, override) {
    if (override) {
        return override
    }
    def rna = meta.molecule != 'dna'
    if (meta.long_reads == 'pacbio') {
        return rna ? 'splice:hq' : 'map-hifi'
    }
    return rna ? 'splice' : 'map-ont'
}

// The chunk a file belongs to, from the name HISAT2_CHUNKSPLIT gives it:
// `<id>.part0001_1.fastq.gz` is part0001. Anchored at the end, so dots or
// underscores in the sample id cannot shift it.
def chunkName(fastq) {
    def match = fastq.name =~ /\.(part\d+)_[12]\.fastq\.gz$/
    if (!match) {
        error("Unexpected chunk file name '${fastq.name}': HISAT2_CHUNKSPLIT names them <id>.partNNNN_<mate>.fastq.gz.")
    }
    return match[0][1]
}

// Collect one sample's chunk outputs back under the sample's own meta, released
// as soon as all of that sample's chunks have arrived (see the call site).
def regroupChunks(ch_chunk_outputs, ch_chunk_count, chunk_keys) {
    return ch_chunk_outputs
        .map { meta, files -> [meta.findAll { key, _value -> !(key in chunk_keys) }, files] }
        .combine(ch_chunk_count, by: 0)
        .map { meta, files, n_chunks -> [groupKey(meta, n_chunks as int), files] }
        .groupTuple(remainder: true)
        .map { group, files -> [group.getGroupTarget(), files.flatten().sort { entry -> entry.name }] }
}

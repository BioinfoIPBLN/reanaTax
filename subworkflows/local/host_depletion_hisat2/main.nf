//
// Split reads into a host fraction (sorted, indexed BAM) and a non-host
// fraction (FASTQ) by aligning against the host genome with HISAT2.
//

include { HISAT2_ALIGN             } from '../../../modules/nf-core/hisat2/align/main'
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

    main:

    def ch_multiqc_files = channel.empty()

    //
    // MODULE: Align against the host. `save_unaligned` is what produces the
    // non-host FASTQs this pipeline is really after.
    //
    HISAT2_ALIGN(
        ch_reads,
        ch_index.first(),
        channel.value([[:], []]),
        true,
    )
    ch_multiqc_files = ch_multiqc_files.mix(HISAT2_ALIGN.out.summary.map { _meta, log -> log })

    def ch_bam = channel.empty()
    def ch_bai = channel.empty()
    def ch_qualimap = channel.empty()

    if (save_host_bam) {
        //
        // SUBWORKFLOW: Sort + index + stats, so the host BAM is directly usable
        // (and its mapping rate lands in MultiQC).
        //
        BAM_SORT_STATS_SAMTOOLS(
            HISAT2_ALIGN.out.bam,
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
    def ch_unaligned = HISAT2_ALIGN.out.fastq.map { meta, fastq ->
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
    summary = HISAT2_ALIGN.out.summary // channel: [ val(meta), path(log) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

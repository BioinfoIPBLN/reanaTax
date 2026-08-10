//
// Split reads into a host fraction (sorted, indexed BAM) and a non-host
// fraction (FASTQ) by aligning against the host genome with HISAT2.
//

include { HISAT2_ALIGN             } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS  } from '../../../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { QUALIMAP_BAMQC           } from '../../../modules/nf-core/qualimap/bamqc/main'

workflow HOST_DEPLETION_HISAT2 {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    ch_index // channel: [ val(meta), path(hisat2_index_dir) ]
    ch_fasta // channel: [ val(meta), path(fasta) ]
    save_host_bam // boolean: sort, index and stat the aligned (host) reads
    skip_qualimap // boolean
    qualimap_gff // string: optional GFF/GTF to add feature-level stats, or null

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
    summary = HISAT2_ALIGN.out.summary // channel: [ val(meta), path(log) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

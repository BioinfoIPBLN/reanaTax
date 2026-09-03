//
// Sort, index and collect stats for an alignment, in BAM or in CRAM.
//
// Deliberately the same NAME as the nf-core subworkflow it stands in for, and
// the same emit contract, so every `withName: '.*:BAM_SORT_STATS_SAMTOOLS:...'`
// selector in conf/modules.config keeps meaning what it did. The one difference
// is the two lines below: nf-core's version reads SAMTOOLS_SORT.out.bam and
// nothing else, so with `--output-fmt cram` in ext.args the sort succeeds, the
// CRAM is written, and the subworkflow emits an empty channel - the run then
// silently loses its host alignments. Everything downstream of the sort
// (index, stats, flagstat, idxstats) already handles CRAM given the reference.
//
// Only the host-depletion passes use this. Targeted alignment stays on the
// nf-core subworkflow: its output is one taxon's reads, so the saving would be
// negligible, and its reference is the target genome rather than --host.
//
include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ch_fasta_fai // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    SAMTOOLS_SORT(ch_bam, ch_fasta_fai, '')

    // Exactly one of these carries anything, decided by --output-fmt in
    // ext.args. `mix` rather than a conditional so the format stays a property
    // of the configuration, which is where a user can override it per process.
    def ch_sorted = SAMTOOLS_SORT.out.bam.mix(SAMTOOLS_SORT.out.cram)

    SAMTOOLS_INDEX(ch_sorted)

    def ch_bam_bai = ch_sorted.join(SAMTOOLS_INDEX.out.index, by: [0])

    BAM_STATS_SAMTOOLS(ch_bam_bai, ch_fasta_fai)

    emit:
    bam = ch_sorted // channel: [ val(meta), [ bam|cram ] ]
    index = SAMTOOLS_INDEX.out.index // channel: [ val(meta), [ bai|csi|crai ] ]
    stats = BAM_STATS_SAMTOOLS.out.stats // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}

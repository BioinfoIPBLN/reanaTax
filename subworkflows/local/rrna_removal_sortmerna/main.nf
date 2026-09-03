//
// Strip ribosomal RNA with SortMeRNA before anything else looks at the reads.
//
// Placed BETWEEN trimming and host depletion, which is not where Sola-Leyva et
// al. put it (they run it on the non-host fraction). Two reasons:
//
//   1. Cost. rRNA is 94-99% of a total-RNA/RiboZero library, so depleting it
//      first shrinks HISAT2's input by up to 20x instead of making HISAT2 align
//      a library that is almost entirely rRNA.
//   2. Accounting. read_accounting.py derives `classified_pct` from the LAST
//      HISAT2 pass's unaligned count, which is exactly Kraken2's input. Insert
//      a filter after HISAT2 and that denominator is silently too large;
//      insert it before, and every number downstream stays true by
//      construction.
//
// The trade is that host rRNA is removed before host depletion can count it,
// so `host_removed_reads` gets smaller. That is the honest number: an rRNA read
// is neither host transcriptome nor microbial signal.
//
// Pairs are held together with `--paired_in`, which flags BOTH mates as rRNA
// when EITHER aligns - so the surviving output holds only pairs in which
// neither mate is rRNA. That is the same rule, and the same reasoning, as
// --require_both_mates_unmapped in host depletion: a uniform, slightly
// aggressive loss beats a sample-specific leak. It also means the output stays
// synchronised, so the repair.sh pass a hand-rolled version needs is unnecessary.
//

include { commaList                   } from '../utils_nfcore_reanatax_pipeline'
include { SORTMERNA as SORTMERNA_INDEX } from '../../../modules/nf-core/sortmerna/main'
include { SORTMERNA as SORTMERNA_READS } from '../../../modules/nf-core/sortmerna/main'

workflow RRNA_REMOVAL_SORTMERNA {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    sortmerna_db // string: comma-separated rRNA reference FASTA(s)
    sortmerna_index // string: prebuilt index directory, or null to build one

    main:

    def ch_fastas = channel.value([
        [id: 'rrna_refs'],
        commaList(sortmerna_db).collect { entry -> file(entry, checkIfExists: true) },
    ])

    //
    // MODULE: Build the index once for the whole run.
    //
    // SortMeRNA indexes its references on every invocation unless handed an
    // `--idx-dir`, which on a cohort would repeat the same minutes of work per
    // sample. Note the index format changed in SortMeRNA 6.0 (CMPH -> BBHash),
    // so an index built by 4.x cannot be reused by 6.x/7.x - point
    // --sortmerna_index at one built by the same version, or let this build it.
    //
    def ch_index = channel.empty()

    if (sortmerna_index) {
        ch_index = channel.value([[id: 'rrna_refs'], file(sortmerna_index, checkIfExists: true)])
    }
    else {
        SORTMERNA_INDEX(
            channel.value([[id: 'rrna_refs'], []]),
            ch_fastas,
            channel.value([[:], []]),
        )
        ch_index = SORTMERNA_INDEX.out.index
    }

    //
    // MODULE: One filtering task per sample against the shared index.
    //
    SORTMERNA_READS(
        ch_reads,
        ch_fastas,
        ch_index.first(),
    )

    // Sort so mate 1 is always first, the same contract host depletion emits.
    def ch_filtered = SORTMERNA_READS.out.reads.map { meta, fastq ->
        def files = (fastq instanceof List ? fastq : [fastq]).sort { fq -> fq.name }
        [meta, files]
    }

    emit:
    reads = ch_filtered // channel: [ val(meta), [ path(fastq) ] ]
    log = SORTMERNA_READS.out.log // channel: [ val(meta), path(log) ]
    index = ch_index // channel: [ val(meta), path(idx) ]
    multiqc_files = SORTMERNA_READS.out.log.map { _meta, log -> log } // channel: path(log)
}

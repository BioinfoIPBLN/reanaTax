//
// What share of each taxon's reads is rRNA, beside the minimizer filter's
// verdict. Annotation only: SortMeRNA flags the rRNA reads among those Kraken2
// classified, and nothing is removed or re-filtered. A taxon seen mostly
// through its rRNA looks narrow to the minimizer filter - its reads pile onto
// one operon - so a low-coverage drop with a high rRNA share is a different
// story from one with none. See bin/rrna_share.py.
//

include { commaList                          } from '../utils_nfcore_reanatax_pipeline'
include { SORTMERNA as RRNA_ANNOTATION_INDEX } from '../../../modules/local/sortmerna/main'
include { RRNA_ANNOTATE                      } from '../../../modules/local/rrna/annotate/main'
include { RRNA_SHARE                         } from '../../../modules/local/rrna/share/main'
include { RRNA_SHARE_MERGE                   } from '../../../modules/local/rrna/merge/main'

workflow RRNA_ANNOTATION {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ] - the reads Kraken2 classified
    ch_assignments // channel: [ val(meta), path(txt) ] - Kraken2's per-read output
    ch_reports // channel: [ val(meta), path(txt) ] - Kraken2's per-sample reports
    ch_evidence // channel: [ val(meta), path(tsv) ] - the minimizer filter's evidence, or empty
    sortmerna_db // string: comma-separated rRNA reference FASTA(s)
    sortmerna_index // string: prebuilt index directory, or null to build one

    main:

    def fastas = commaList(sortmerna_db).collect { entry -> file(entry, checkIfExists: true) }

    def ch_index
    if (sortmerna_index) {
        ch_index = channel.value(file(sortmerna_index, checkIfExists: true))
    }
    else {
        RRNA_ANNOTATION_INDEX(
            channel.value([[id: 'rrna_refs'], []]),
            channel.value([[id: 'rrna_refs'], fastas]),
            channel.value([[:], []]),
        )
        ch_index = RRNA_ANNOTATION_INDEX.out.index.map { _meta, index -> index }.first()
    }

    RRNA_ANNOTATE(
        ch_reads,
        channel.value(fastas),
        ch_index,
    )

    // Keyed on meta.id: the Kraken2 route adds fields to the meta map, so the
    // whole maps would not match.
    RRNA_SHARE(
        RRNA_ANNOTATE.out.ids.map { meta, ids -> [meta.id, meta, ids] }
            .join(ch_assignments.map { meta, txt -> [meta.id, txt] })
            .join(ch_reports.map { meta, report -> [meta.id, report] })
            .map { _id, meta, ids, txt, report -> [meta, ids, txt, report] }
    )

    RRNA_SHARE_MERGE(
        RRNA_SHARE.out.tsv
            .map { _meta, tsv -> tsv }
            .collect(sort: true)
            .map { tables -> [[id: 'reanatax'], tables] },
        ch_evidence.map { _meta, tsv -> tsv }.collect().ifEmpty([]),
    )

    emit:
    tsv = RRNA_SHARE_MERGE.out.tsv // channel: [ val(meta), path(tsv) ]
    ids = RRNA_ANNOTATE.out.ids // channel: [ val(meta), path(txt.gz) ]
}

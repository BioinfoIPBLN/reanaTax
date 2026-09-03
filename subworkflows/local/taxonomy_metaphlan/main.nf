//
// MetaPhlAn profiling of the non-host reads, and one merged table for the run.
//

include { METAPHLAN_PROFILE } from '../../../modules/local/metaphlan/profile/main'
include { METAPHLAN_MERGE   } from '../../../modules/local/metaphlan/merge/main'

workflow TAXONOMY_METAPHLAN {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    metaphlan_db // string: directory holding the bowtie2 index and its .pkl
    metaphlan_index // string: index basename inside that directory

    main:

    def ch_db = channel.value(file(metaphlan_db, checkIfExists: true))

    METAPHLAN_PROFILE(ch_reads, ch_db, metaphlan_index)

    // Sorted so the column order of the merged table is reproducible.
    METAPHLAN_MERGE(
        METAPHLAN_PROFILE.out.profile
            .map { _meta, profile -> profile }
            .collect(sort: true)
            .map { profiles -> [[id: 'metaphlan_combined'], profiles] }
    )

    emit:
    profile = METAPHLAN_PROFILE.out.profile // channel: [ val(meta), path(txt) ]
    merged = METAPHLAN_MERGE.out.merged // channel: [ val(meta), path(tsv) ]
}

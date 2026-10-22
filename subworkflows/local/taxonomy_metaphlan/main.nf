//
// MetaPhlAn profiling of the non-host reads, and one merged table for the run.
//

include { METAPHLAN_PROFILE } from '../../../modules/local/metaphlan/profile/main'
include { METAPHLAN_MERGE   } from '../../../modules/local/metaphlan/merge/main'
include { METAPHLAN_MM2INDEX } from '../../../modules/local/metaphlan/mm2index/main'

workflow TAXONOMY_METAPHLAN {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    metaphlan_db // string: directory holding the bowtie2 index and its .pkl
    metaphlan_index // string: index basename inside that directory
    long_reads // boolean: the libraries are long reads, mapped with minimap2
    minimap2_ps // string: MetaPhlAn's --minimap2_ps, for the long-read index and mapping
    minimap2_index // string: a prebuilt index of the markers for those presets, or null to build one

    main:

    def ch_db = channel.value(file(metaphlan_db, checkIfExists: true))

    // MetaPhlAn's minimap2 index, built once for every long-read sample (see
    // METAPHLAN_MM2INDEX). Short reads go through Bowtie2 and take none.
    def ch_mm2_index = channel.value([])
    if (long_reads && minimap2_index) {
        ch_mm2_index = channel.value(file(minimap2_index, checkIfExists: true))
    }
    else if (long_reads) {
        METAPHLAN_MM2INDEX(ch_db, metaphlan_index, minimap2_ps)
        ch_mm2_index = METAPHLAN_MM2INDEX.out.index
    }

    METAPHLAN_PROFILE(ch_reads, ch_db, metaphlan_index, ch_mm2_index, minimap2_ps)

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

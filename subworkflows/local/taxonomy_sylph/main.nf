//
// sylph profiling of the non-host reads: per-sample genome calls, an optional
// sylph-tax taxonomy on top, and cohort tables for the run.
//

include { SYLPH_PROFILE } from '../../../modules/local/sylph/profile/main'
include { SYLPH_TAXPROF } from '../../../modules/local/sylph/taxprof/main'
include { SYLPH_MERGE   } from '../../../modules/local/sylph/merge/main'

workflow TAXONOMY_SYLPH {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    sylph_db // string: one or more .syldb / .syl2db databases, comma-separated
    sylph_taxonomy // string: sylph-tax metadata for those databases, comma-separated; null skips the taxonomy

    main:

    // Value channels. Every sample is profiled against the same databases, and
    // a queue channel holding one item would pair with the first sample only.
    def ch_db = channel.value(
        sylph_db.tokenize(',').collect { db -> file(db.trim(), checkIfExists: true) }
    )

    SYLPH_PROFILE(ch_reads, ch_db)

    def ch_taxprof = channel.empty()
    if (sylph_taxonomy) {
        def ch_taxonomy = channel.value(
            sylph_taxonomy.tokenize(',').collect { tax -> file(tax.trim(), checkIfExists: true) }
        )
        SYLPH_TAXPROF(SYLPH_PROFILE.out.profile, ch_taxonomy)
        ch_taxprof = SYLPH_TAXPROF.out.profile
    }

    // Sorted so the cohort tables come out in a reproducible order. The
    // sylph-tax list is wrapped before `combine` so it stays one element of
    // the tuple instead of having its files spliced in.
    SYLPH_MERGE(
        SYLPH_PROFILE.out.profile
            .map { _meta, profile -> profile }
            .collect(sort: true)
            .map { profiles -> [[id: 'sylph_combined'], profiles] }
            .combine(
                ch_taxprof
                    .map { _meta, profile -> profile }
                    .collect(sort: true)
                    .ifEmpty([])
                    .map { profiles -> [profiles] }
            )
    )

    emit:
    profile = SYLPH_PROFILE.out.profile // channel: [ val(meta), path(tsv) ]
    taxprof = ch_taxprof // channel: [ val(meta), path(sylphmpa) ]
    genomes = SYLPH_MERGE.out.genomes // channel: [ val(meta), path(tsv) ]
    merged = SYLPH_MERGE.out.merged // channel: [ val(meta), [ path(tsv) ] ]
}

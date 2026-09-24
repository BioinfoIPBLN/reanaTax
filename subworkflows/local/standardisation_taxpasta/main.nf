//
// Standardise every Taxpasta-readable profiler's per-sample outputs into one
// table per profiler: taxonomy_id plus a count column per sample.
//
// Four of this pipeline's profilers are covered - Kraken2, Bracken, KrakenUniq
// and MetaPhlAn - because those are the ones Taxpasta 0.7 reads. sylph, Metax
// and PathSeq are not, and write their own cohort tables instead.
//
// These are the RAW per-sample outputs, before this pipeline's evidence filters
// (host clade, host k-mers, minimizers, negative controls). Taxpasta reads each
// profiler's native per-sample format, which is what exists before filtering;
// the filtered views stay in each profiler's own combined tables.
//

include { TAXPASTA_MERGE as TAXPASTA_KRAKEN2    } from '../../../modules/local/taxpasta/merge/main'
include { TAXPASTA_MERGE as TAXPASTA_BRACKEN    } from '../../../modules/local/taxpasta/merge/main'
include { TAXPASTA_MERGE as TAXPASTA_KRAKENUNIQ } from '../../../modules/local/taxpasta/merge/main'
include { TAXPASTA_MERGE as TAXPASTA_METAPHLAN  } from '../../../modules/local/taxpasta/merge/main'

workflow STANDARDISATION_TAXPASTA {

    take:
    ch_kraken2 // channel: [ val(meta), path(kreport) ]
    ch_bracken // channel: [ val(meta), path(bracken tsv) ]
    ch_krakenuniq // channel: [ val(meta), path(report) ]
    ch_metaphlan // channel: [ val(meta), path(profile) ]
    taxonomy_dir // string: taxdump directory for --add-name and friends, or null

    main:

    def ch_taxonomy = taxonomy_dir ? channel.value(file(taxonomy_dir, checkIfExists: true)) : []

    TAXPASTA_KRAKEN2(gather(ch_kraken2, 'kraken2'), 'kraken2', ch_taxonomy)
    TAXPASTA_BRACKEN(gather(ch_bracken, 'bracken'), 'bracken', ch_taxonomy)
    TAXPASTA_KRAKENUNIQ(gather(ch_krakenuniq, 'krakenuniq'), 'krakenuniq', ch_taxonomy)
    TAXPASTA_METAPHLAN(gather(ch_metaphlan, 'metaphlan'), 'metaphlan', ch_taxonomy)

    emit:
    tables = TAXPASTA_KRAKEN2.out.table
        .mix(TAXPASTA_BRACKEN.out.table, TAXPASTA_KRAKENUNIQ.out.table, TAXPASTA_METAPHLAN.out.table) // channel: [ val(meta), path(table) ]
}

/**
* One Taxpasta input per profiler: every sample's output, sorted by sample id so
* the column order is reproducible, with the ids and files kept in step. A
* profiler that did not run yields nothing, so its Taxpasta task never starts.
*/
def gather(ch_profiles, tool) {
    return ch_profiles
        .map { meta, profile -> [meta.id, profile] }
        .toSortedList { a, b -> a[0] <=> b[0] }
        .filter { pairs -> !pairs.isEmpty() }
        .map { pairs -> [[id: tool], pairs.collect { pair -> pair[0] }, pairs.collect { pair -> pair[1] }] }
}

//
// Metax profiling of the non-host reads: per-sample profiles, filtered on
// genome-wide coverage, and cohort tables for the run.
//

include { METAX_PROFILE } from '../../../modules/local/metax/profile/main'
include { METAX_MERGE   } from '../../../modules/local/metax/merge/main'

workflow TAXONOMY_METAX {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    metax_db // string: the database's metax_db.json, or the directory holding it
    metax_dmp_dir // string: taxonomy dump directory (nodes.dmp, names.dmp, merged.dmp)
    save_classify // boolean: publish the per-read classify.txt

    main:

    // maCMD finds the index files through the JSON, so the directory is staged
    // whole and the JSON is named inside it rather than staged on its own.
    def db = file(metax_db, checkIfExists: true)
    def db_json = db.isDirectory() ? db.resolve('metax_db.json') : db
    if (!db_json.exists()) {
        error("--metax_db: no metax_db.json at ${db_json}. Point it at the JSON of a Metax database, or at the directory that holds it.")
    }

    // Value channels, so every sample pairs with the one database and taxonomy.
    METAX_PROFILE(
        ch_reads,
        channel.value(db_json.parent),
        db_json.name,
        channel.value(file(metax_dmp_dir, checkIfExists: true)),
        save_classify,
    )

    // Sorted so the column order of the cohort tables is reproducible.
    METAX_MERGE(
        METAX_PROFILE.out.profile
            .map { _meta, profile -> profile }
            .collect(sort: true)
            .map { profiles -> [[id: 'metax_combined'], profiles] }
    )

    emit:
    profile = METAX_PROFILE.out.profile // channel: [ val(meta), path(txt) ]
    classify = METAX_PROFILE.out.classify // channel: [ val(meta), path(txt.gz) ]
    combined = METAX_MERGE.out.combined // channel: [ val(meta), path(tsv) ]
    merged = METAX_MERGE.out.merged // channel: [ val(meta), [ path(tsv) ] ]
}

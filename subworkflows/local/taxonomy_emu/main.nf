//
// Emu profiling of the non-host reads of a long-read run: per-sample species
// profiles and cohort tables at the ranks asked for.
//

include { UNTAR as UNTAR_EMU } from '../../../modules/nf-core/untar/main'
include { EMU_ABUNDANCE      } from '../../../modules/local/emu/abundance/main'
include { EMU_COMBINE        } from '../../../modules/local/emu/combine/main'
include { isIndexArchive     } from '../utils_nfcore_reanatax_pipeline'

workflow TAXONOMY_EMU {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    emu_db // string: directory holding species_taxid.fasta and taxonomy.tsv, or a .tar/.tar.gz of one
    ranks // list: Emu ranks to write cohort tables at
    counts // boolean: EMU_ABUNDANCE keeps estimated counts, so count tables are written too

    main:

    // A value channel, so that every sample is profiled against it; Emu's
    // prebuilt databases are distributed as tar archives.
    def ch_db = channel.empty()
    if (isIndexArchive(emu_db)) {
        UNTAR_EMU(channel.value([[id: 'emu_db'], file(emu_db, checkIfExists: true)]))
        ch_db = UNTAR_EMU.out.untar.map { _meta, dir -> dir }
    }
    else {
        ch_db = channel.value(file(emu_db, checkIfExists: true))
    }

    EMU_ABUNDANCE(ch_reads, ch_db)

    // Sorted so the column order of the cohort tables is reproducible.
    EMU_COMBINE(
        EMU_ABUNDANCE.out.profile
            .map { _meta, profile -> profile }
            .collect(sort: true)
            .map { profiles -> [[id: 'emu_combined'], profiles] },
        ranks,
        counts,
    )

    emit:
    profile = EMU_ABUNDANCE.out.profile // channel: [ val(meta), path(tsv) ]
    combined = EMU_COMBINE.out.combined // channel: [ val(meta), [ path(tsv) ] ]
}

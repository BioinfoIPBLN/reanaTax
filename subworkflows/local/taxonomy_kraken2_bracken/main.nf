//
// Taxonomic classification of the non-host fraction with Kraken2, abundance
// re-estimation with Bracken and interactive Krona charts.
//

include { UNTAR                          } from '../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2                } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN                } from '../../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS  } from '../../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { KRAKENTOOLS_KREPORT2KRONA      } from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENTOOLS_COMBINEKREPORTS    } from '../../../modules/nf-core/krakentools/combinekreports/main'
include { KRONA_KTIMPORTTEXT             } from '../../../modules/nf-core/krona/ktimporttext/main'

workflow TAXONOMY_KRAKEN2_BRACKEN {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    kraken2_db // string: path to a Kraken2 database directory or tarball
    bracken_db // string: path to the Bracken kmer distributions, or null to reuse kraken2_db
    save_output_fastqs // boolean: keep the classified/unclassified FASTQs
    save_reads_assignment // boolean: keep the per-read assignment table
    skip_bracken // boolean
    skip_krona // boolean

    main:

    def ch_multiqc_files = channel.empty()

    //
    // Databases are large; accept either an unpacked directory (the usual case
    // on a cluster) or a tarball (handy for CI and for object storage).
    //
    def ch_kraken2_db = channel.empty()
    if (kraken2_db.endsWith('.tar.gz') || kraken2_db.endsWith('.tgz')) {
        UNTAR(channel.value([[id: 'kraken2_db'], file(kraken2_db, checkIfExists: true)]))
        ch_kraken2_db = UNTAR.out.untar.map { _meta, db -> db }
    }
    else {
        ch_kraken2_db = channel.value(file(kraken2_db, checkIfExists: true))
    }

    def ch_bracken_db = bracken_db
        ? channel.value(file(bracken_db, checkIfExists: true))
        : ch_kraken2_db

    //
    // MODULE: Kraken2
    //
    KRAKEN2_KRAKEN2(
        ch_reads,
        ch_kraken2_db,
        save_output_fastqs,
        save_reads_assignment,
    )
    ch_multiqc_files = ch_multiqc_files.mix(KRAKEN2_KRAKEN2.out.report.map { _meta, report -> report })

    //
    // A sample in which nothing was classified produces a report holding only
    // the `U` (unclassified) row. Both Bracken and combine_kreports.py crash on
    // those, so they are dropped here - loudly - rather than being allowed to
    // take the whole run down. The report itself is still published and still
    // reaches MultiQC.
    //
    def ch_report_classified = KRAKEN2_KRAKEN2.out.report.filter { meta, report ->
        def classified = hasClassifiedReads(report)
        if (!classified) {
            log.warn("Kraken2 classified no reads for '${meta.id}'; excluding it from Bracken and from the combined tables.")
        }
        classified
    }

    //
    // MODULE: Combine every sample's Kraken2 report into one table.
    // Sorting by sample id both makes the column order reproducible and lets the
    // sample names be handed to the tool in the same order as the files, so the
    // combined table is labelled with sample ids rather than file names.
    //
    KRAKENTOOLS_COMBINEKREPORTS(
        ch_report_classified
            .toSortedList { entry_a, entry_b -> entry_a[0].id <=> entry_b[0].id }
            .filter { entries -> entries }
            .map { entries ->
                [
                    [id: 'kraken2_combined', names: entries.collect { entry -> entry[0].id }.join(' ')],
                    entries.collect { entry -> entry[1] },
                ]
            }
    )

    def ch_bracken = channel.empty()
    def ch_bracken_combined = channel.empty()

    if (!skip_bracken) {
        //
        // MODULE: Bracken re-estimates abundances from the Kraken2 report
        //
        BRACKEN_BRACKEN(ch_report_classified, ch_bracken_db)
        ch_bracken = BRACKEN_BRACKEN.out.reports

        BRACKEN_COMBINEBRACKENOUTPUTS(
            BRACKEN_BRACKEN.out.reports
                .toSortedList { entry_a, entry_b -> entry_a[0].id <=> entry_b[0].id }
                .filter { entries -> entries }
                .map { entries ->
                    [
                        [id: 'bracken_combined', names: entries.collect { entry -> entry[0].id }.join(',')],
                        entries.collect { entry -> entry[1] },
                    ]
                }
        )
        ch_bracken_combined = BRACKEN_COMBINEBRACKENOUTPUTS.out.txt
    }

    def ch_krona = channel.empty()

    if (!skip_krona) {
        //
        // Krona is rendered from the Bracken-corrected report when Bracken ran,
        // because that is the abundance estimate users are meant to interpret.
        //
        // Krona copes with an all-unclassified report, so it keeps every sample.
        def ch_for_krona = skip_bracken ? KRAKEN2_KRAKEN2.out.report : BRACKEN_BRACKEN.out.txt

        KRAKENTOOLS_KREPORT2KRONA(ch_for_krona)
        KRONA_KTIMPORTTEXT(KRAKENTOOLS_KREPORT2KRONA.out.txt)
        ch_krona = KRONA_KTIMPORTTEXT.out.html
    }

    emit:
    report = KRAKEN2_KRAKEN2.out.report // channel: [ val(meta), path(report) ]
    report_combined = KRAKENTOOLS_COMBINEKREPORTS.out.txt // channel: [ val(meta), path(txt) ]
    bracken = ch_bracken // channel: [ val(meta), path(tsv) ]
    bracken_combined = ch_bracken_combined // channel: [ val(meta), path(txt) ]
    krona = ch_krona // channel: [ val(meta), path(html) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

//
// True if a Kraken2 report contains at least one taxon, i.e. anything beyond the
// `U` row. Matching on the rank-code column rather than a fixed index keeps this
// working with the extra columns `--report-minimizer-data` adds. A Kraken2
// report is one line per taxon, so even against a large database this is a few
// MB, and `any` stops scanning at the first classified row.
//
def hasClassifiedReads(report) {
    return report.readLines().any { line ->
        line.split('\t').any { column -> column.trim() ==~ /^[RDKPCOFGS][0-9]*$/ }
    }
}

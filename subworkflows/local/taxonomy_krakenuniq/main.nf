//
// Taxonomic classification with KrakenUniq.
//
// KrakenUniq's reason to exist is the `kmers`/`dup`/`cov` columns: an estimate,
// by HyperLogLog, of how many DISTINCT k-mers of a taxon the reads actually
// cover. A taxon called from reads piled on one conserved stretch has few
// distinct k-mers for its read count and low coverage of its reference, which
// abundance alone cannot show. Kraken2's `--report-minimizer-data` answers the
// same question with distinct MINIMIZERS - a subsample of roughly one in five
// at the default k=35/l=31 - so this route buys resolution at low abundance,
// where the minimizer count gets too thin to discriminate.
//
// It needs its own database: KrakenUniq cannot read a Kraken2 index. Hence
// --krakenuniq_db rather than a switch on the existing one, and hence this
// route being entirely opt-in.
//
// Bracken is deliberately NOT run here. Bracken's re-estimation needs the k-mer
// distributions built alongside a Kraken2 database, and pairing those with a
// KrakenUniq report would be estimating one database's abundances from
// another's statistics.
//

include { UNTAR as UNTAR_KRAKENUNIQ } from '../../../modules/nf-core/untar/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                 } from '../../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { KRAKENTOOLS_COMBINEKREPORTS as COMBINE_KRAKENUNIQ } from '../../../modules/nf-core/krakentools/combinekreports/main'
include { KRAKENTOOLS_KREPORT2KRONA as KREPORT2KRONA_UNIQ   } from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT as KRONA_KRAKENUNIQ            } from '../../../modules/nf-core/krona/ktimporttext/main'
include { ABUNDANCE_MINIMIZER as MINIMIZER_KRAKENUNIQ       } from '../../../modules/local/abundance/minimizer/main'
include { ABUNDANCE_FILTER as FILTER_KRAKENUNIQ             } from '../../../modules/local/abundance/filter/main'
include { hasClassifiedReads                                } from '../taxonomy_kraken2_bracken'

workflow TAXONOMY_KRAKENUNIQ {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    krakenuniq_db // string: path to a KrakenUniq database directory, or a .tar.gz of one
    save_output_fastqs // boolean
    save_reads_assignment // boolean
    skip_krona // boolean
    minimizer_filter // boolean: drop taxa whose evidence has no breadth
    minimizer_thresholds // map: min_reads, min_distinct, max_duplication, min_coverage
    min_rel_abundance // float
    min_samples // integer
    min_reads // number: ...and this many reads in total, 0 to disable

    main:

    // A tarball is unpacked first, as on the Kraken2 route - which is what lets
    // the test profiles use nf-core's packed SARS-CoV-2 KrakenUniq database.
    // .first() for the same reason given there: every batch needs the database,
    // not just the first.
    def ch_db = channel.empty()
    if (krakenuniq_db.endsWith('.tar.gz') || krakenuniq_db.endsWith('.tgz')) {
        UNTAR_KRAKENUNIQ(channel.value([[id: 'krakenuniq_db'], file(krakenuniq_db, checkIfExists: true)]))
        ch_db = UNTAR_KRAKENUNIQ.out.untar.map { _meta, db -> db }.first()
    }
    else {
        ch_db = channel.value(file(krakenuniq_db, checkIfExists: true))
    }

    //
    // MODULE: One task per layout, not per sample.
    //
    // The module classifies a whole batch after a single `--preload`, which is
    // the only sane way to use a database this size: loaded once, then every
    // sample in the batch classified against the resident copy. Single- and
    // paired-end libraries go in separate batches because the module builds a
    // different command for each and asserts on the count of files per sample.
    //
    def ch_batches = ch_reads
        .map { meta, reads -> [meta.single_end, meta.id, reads instanceof List ? reads : [reads]] }
        .groupTuple(by: 0)
        .map { single_end, ids, reads ->
            // Sort by sample id so the batch is byte-identical between resumes.
            def ordered = [ids, reads].transpose().sort { entry -> entry[0] }
            [
                [id: single_end ? 'krakenuniq_single' : 'krakenuniq_paired', single_end: single_end],
                ordered.collect { entry -> entry[1] }.flatten(),
                ordered.collect { entry -> entry[0] },
            ]
        }

    KRAKENUNIQ_PRELOADEDKRAKENUNIQ(
        ch_batches,
        'fastq',
        ch_db,
        save_output_fastqs,
        true,
        save_reads_assignment,
    )

    //
    // The batch comes back as one list of reports, so the per-sample identity
    // has to be recovered from the file names the module wrote.
    //
    def ch_report = KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report.flatMap { _meta, reports ->
        (reports instanceof List ? reports : [reports]).collect { report ->
            [[id: report.name.replaceAll(/\.krakenuniq\.report\.txt$/, '')], report]
        }
    }

    def ch_report_classified = ch_report.filter { meta, report ->
        def classified = hasClassifiedReads(report)
        if (!classified) {
            log.warn("KrakenUniq classified no reads for '${meta.id}'; excluding it from the combined table.")
        }
        classified
    }

    //
    // MODULE: Combine into one table. KrakenTools reads KrakenUniq's layout -
    // taxid before rank, ranks spelled out - through its own `map_kuniq` branch.
    //
    COMBINE_KRAKENUNIQ(
        ch_report_classified
            .toSortedList { entry_a, entry_b -> entry_a[0].id <=> entry_b[0].id }
            .filter { entries -> entries }
            .map { entries ->
                [
                    [id: 'krakenuniq_combined', names: entries.collect { entry -> entry[0].id }.join(' ')],
                    entries.collect { entry -> entry[1] },
                ]
            }
    )

    //
    // MODULE: Which taxa have no breadth behind their reads. KrakenUniq needs
    // no database lookup for this - it reports distinct k-mers, duplication and
    // reference coverage in the report itself.
    //
    def ch_drop_list = channel.value([])
    def ch_evidence = channel.empty()

    if (minimizer_filter) {
        MINIMIZER_KRAKENUNIQ(
            ch_report_classified
                .map { _meta, report -> report }
                .collect(sort: true)
                .map { reports -> [[id: 'krakenuniq_minimizer'], reports] },
            [],
            minimizer_thresholds.min_reads,
            minimizer_thresholds.min_distinct,
            minimizer_thresholds.max_duplication,
            minimizer_thresholds.min_coverage,
            minimizer_thresholds.distinct_scale,
        )
        ch_drop_list = MINIMIZER_KRAKENUNIQ.out.drop_list.map { _meta, list -> list }.first()
        ch_evidence = MINIMIZER_KRAKENUNIQ.out.evidence
    }

    FILTER_KRAKENUNIQ(COMBINE_KRAKENUNIQ.out.txt, ch_drop_list, min_rel_abundance, min_samples, min_reads)

    def ch_krona = channel.empty()

    if (!skip_krona) {
        KREPORT2KRONA_UNIQ(ch_report)
        KRONA_KRAKENUNIQ(KREPORT2KRONA_UNIQ.out.txt)
        ch_krona = KRONA_KRAKENUNIQ.out.html
    }

    emit:
    report = ch_report // channel: [ val(meta), path(report) ]
    report_combined = COMBINE_KRAKENUNIQ.out.txt // channel: [ val(meta), path(txt) ]
    report_combined_filtered = FILTER_KRAKENUNIQ.out.filtered // channel: [ val(meta), path(tsv) ]
    krona = ch_krona // channel: [ val(meta), path(html) ]
    minimizer_evidence = ch_evidence // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_report.map { _meta, report -> report } // channel: path(report)
}

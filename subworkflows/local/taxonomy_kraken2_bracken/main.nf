//
// Taxonomic classification of the non-host fraction with Kraken2, abundance
// re-estimation with Bracken and interactive Krona charts.
//

include { UNTAR                          } from '../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2                } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_DAEMON                } from '../../../modules/local/kraken2/daemon/main'
include { BRACKEN_BRACKEN                } from '../../../modules/nf-core/bracken/bracken/main'
include { BRACKEN_COMBINEBRACKENOUTPUTS  } from '../../../modules/nf-core/bracken/combinebrackenoutputs/main'
include { KRAKENTOOLS_KREPORT2KRONA      } from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { KRAKENTOOLS_COMBINEKREPORTS    } from '../../../modules/nf-core/krakentools/combinekreports/main'
include { KRONA_KTIMPORTTEXT             } from '../../../modules/nf-core/krona/ktimporttext/main'
include { ABUNDANCE_MINIMIZER               } from '../../../modules/local/abundance/minimizer/main'
include { HOSTKMER_SCAN                     } from '../../../modules/local/hostkmer/scan/main'
include { HOSTKMER_MERGE                    } from '../../../modules/local/hostkmer/merge/main'
include { DECONTAM_FILTER                   } from '../../../modules/local/decontam/main'
include { SHUFFLE_READS                     } from '../../../modules/local/shuffle/reads/main'
include { SHUFFLE_COMPARE                   } from '../../../modules/local/shuffle/compare/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_SHUFFLED } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_DAEMON as KRAKEN2_DAEMON_SHUFFLED } from '../../../modules/local/kraken2/daemon/main'
include { ABUNDANCE_CONTROL as CONTROL_KRAKEN2 } from '../../../modules/local/abundance/control/main'
include { ABUNDANCE_CONTROL as CONTROL_BRACKEN } from '../../../modules/local/abundance/control/main'
include { ABUNDANCE_FILTER as FILTER_KRAKEN2 } from '../../../modules/local/abundance/filter/main'
include { ABUNDANCE_FILTER as FILTER_BRACKEN } from '../../../modules/local/abundance/filter/main'
include { KRAKENBIOM                        } from '../../../modules/local/krakenbiom/main'
include { MEGAHIT                           } from '../../../modules/local/megahit/main'
include { SPADES_META                       } from '../../../modules/local/spades/meta/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CONTIGS } from '../../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_DAEMON as KRAKEN2_DAEMON_CONTIGS } from '../../../modules/local/kraken2/daemon/main'
include { CONTIG_EVIDENCE                   } from '../../../modules/local/contig/evidence/main'
include { HOSTCLADE_EXPAND                  } from '../../../modules/local/hostclade/main'

workflow TAXONOMY_KRAKEN2_BRACKEN {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    kraken2_db // string: path to a Kraken2 database directory or tarball
    bracken_db // string: path to the Bracken kmer distributions, or null to reuse kraken2_db
    save_output_fastqs // boolean: keep the classified/unclassified FASTQs
    save_reads_assignment // boolean: keep the per-read assignment table
    skip_bracken // boolean
    skip_krona // boolean
    minimizer_filter // boolean: drop taxa whose evidence has no breadth
    minimizer_thresholds // map: min_reads, min_distinct, max_duplication, min_coverage
    host_kmer_filter // boolean: drop taxa whose reads are mostly host k-mers
    host_kmer_settings // map: taxid, max_fraction, min_reads
    decontam_settings // map: metadata, neg_column, neg_value, conc_column, method, threshold, batch_column, batch_combine - or null
    shuffle_settings // map: method, seed, max_reads, max_ratio, min_reads - or null
    control_settings // map: controls, ratio, statistic, min_reads, floor_reads, prevalence, prevalence_min_reads - or null
    assembly_settings // map: assembler, pool, groups, min_reads, min_length, min_contigs, filter - or null
    host_clade // map: taxid, rank - or null
    use_daemon // boolean: classify through `k2 classify --use-daemon`
    min_rel_abundance // float: relative abundance a taxon must exceed...
    min_samples // integer: ...in at least this many samples
    min_reads // number: ...and this many reads in total, 0 to disable
    export_biom // boolean: also write the combined counts as a BIOM table
    biom_metadata // path: sample metadata to embed in the BIOM, or null

    main:

    def ch_multiqc_files = channel.empty()

    //
    // Databases are large; accept either an unpacked directory (the usual case
    // on a cluster) or a tarball (handy for CI and for object storage).
    //
    def ch_kraken2_db = channel.empty()
    if (kraken2_db.endsWith('.tar.gz') || kraken2_db.endsWith('.tgz')) {
        UNTAR(channel.value([[id: 'kraken2_db'], file(kraken2_db, checkIfExists: true)]))
        // .first() makes it a VALUE channel. Without it this is a queue channel
        // holding one item, and a process pairs a queue channel element-by-element
        // against the reads - so the cohort would classify its first sample and
        // silently stop. The directory branch below is already a value channel;
        // there are now three consumers (reads, shuffled reads, contigs) and all
        // of them need the database once per task, not once per run.
        ch_kraken2_db = UNTAR.out.untar.map { _meta, db -> db }.first()
    }
    else {
        ch_kraken2_db = channel.value(file(kraken2_db, checkIfExists: true))
    }

    def ch_bracken_db = bracken_db
        ? channel.value(file(bracken_db, checkIfExists: true))
        : ch_kraken2_db

    //
    // MODULE: Kraken2.
    //
    // The daemon variant hands the database to a resident background process
    // instead of loading it per task, which is the difference between paying
    // the index load once and paying it per sample. It needs the database as an
    // unstaged absolute path - see modules/local/kraken2/daemon - so it cannot
    // take the untarred channel, and a tarball is rejected up front.
    //
    def ch_kraken2_report = channel.empty()
    def ch_classifiedreads = channel.empty()

    if (use_daemon) {
        KRAKEN2_DAEMON(
            ch_reads,
            file(kraken2_db, checkIfExists: true).toAbsolutePath().toString(),
            save_output_fastqs,
            save_reads_assignment,
        )
        ch_kraken2_report = KRAKEN2_DAEMON.out.report
        ch_classifiedreads = KRAKEN2_DAEMON.out.classified_reads_assignment
    }
    else {
        KRAKEN2_KRAKEN2(
            ch_reads,
            ch_kraken2_db,
            save_output_fastqs,
            save_reads_assignment,
        )
        ch_kraken2_report = KRAKEN2_KRAKEN2.out.report
        ch_classifiedreads = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    }
    ch_multiqc_files = ch_multiqc_files.mix(ch_kraken2_report.map { _meta, report -> report })

    //
    // A sample in which nothing was classified produces a report holding only
    // the `U` (unclassified) row. Both Bracken and combine_kreports.py crash on
    // those, so they are dropped here - loudly - rather than being allowed to
    // take the whole run down. The report itself is still published and still
    // reaches MultiQC.
    //
    def ch_report_classified = ch_kraken2_report.filter { meta, report ->
        def classified = hasClassifiedReads(report)
        if (!classified) {
            log.warn("Kraken2 classified no reads for '${meta.id}'; excluding it from Bracken and from the combined tables.")
        }
        classified
    }

    //
    // MODULE: The shuffled-read negative control.
    //
    // The same libraries, classified twice: once as sequenced and once with
    // every read shuffled so that its length, GC content and dinucleotide
    // frequencies survive and none of its k-mers do. Whatever the database
    // reports from the second pass is what it produces from base composition
    // alone, and a taxon that scores comparably on both is not supported by
    // homology at all.
    //
    // This is the only check here with an external null. Every other filter
    // asks whether a taxon's evidence looks strong on its own terms; this one
    // measures what "strong" is worth against reads that contain no sequence.
    //
    // The shuffled reports deliberately keep the same file names as the real
    // ones - SHUFFLE_COMPARE pairs them by basename - and are staged into
    // separate directories to keep them apart. They are NOT mixed into MultiQC,
    // which keys on file name and would treat them as duplicate samples.
    //
    def ch_shuffle_drop = channel.empty()
    def ch_shuffle_evidence = channel.empty()

    if (shuffle_settings) {
        SHUFFLE_READS(
            ch_reads,
            shuffle_settings.method,
            shuffle_settings.seed,
            shuffle_settings.max_reads,
        )

        def ch_shuffled_report = channel.empty()
        if (use_daemon) {
            KRAKEN2_DAEMON_SHUFFLED(
                SHUFFLE_READS.out.reads,
                file(kraken2_db, checkIfExists: true).toAbsolutePath().toString(),
                false,
                false,
            )
            ch_shuffled_report = KRAKEN2_DAEMON_SHUFFLED.out.report
        }
        else {
            KRAKEN2_SHUFFLED(SHUFFLE_READS.out.reads, ch_kraken2_db, false, false)
            ch_shuffled_report = KRAKEN2_SHUFFLED.out.report
        }

        // Joined on sample id rather than zipped: a sample whose shuffled copy
        // classified nothing at all would otherwise shift every later pairing
        // by one and silently compare the wrong two libraries.
        SHUFFLE_COMPARE(
            ch_report_classified
                .map { meta, report -> [meta.id, report] }
                .join(ch_shuffled_report.map { meta, report -> [meta.id, report] })
                .toSortedList { entry_a, entry_b -> entry_a[0] <=> entry_b[0] }
                .filter { entries -> entries }
                .map { entries ->
                    [
                        [id: 'kraken2_shuffle'],
                        entries.collect { entry -> entry[1] },
                        entries.collect { entry -> entry[2] },
                    ]
                },
            shuffle_settings.max_ratio,
            shuffle_settings.min_reads,
        )
        ch_shuffle_drop = SHUFFLE_COMPARE.out.drop_list.map { _meta, list -> list }
        ch_shuffle_evidence = SHUFFLE_COMPARE.out.evidence
        ch_multiqc_files = ch_multiqc_files.mix(SHUFFLE_COMPARE.out.mqc.map { _meta, mqc -> mqc })
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

    //
    // MODULE: Drop taxa too sparse for any downstream statistic to speak about.
    // Applied to the COMBINED tables only: the per-sample reports stay intact,
    // because prevalence is a cross-sample property and filtering a sample in
    // isolation would just be a detection threshold.
    //
    //
    // MODULE: Which taxa have no breadth behind their reads.
    //
    // Reads on one conserved locus look abundant, so an abundance threshold
    // cannot see them; the distinct-minimizer columns can. The taxids it
    // condemns are handed to the abundance filter rather than removed here, so
    // one step owns every removal from the combined tables and one file records
    // them all.
    //
    def ch_minimizer_drop = channel.empty()
    def ch_hostkmer_drop = channel.empty()
    def ch_hostkmer_evidence = channel.empty()

    if (minimizer_filter) {
        ABUNDANCE_MINIMIZER(
            ch_report_classified
                .map { _meta, report -> report }
                .collect(sort: true)
                .map { reports -> [[id: 'kraken2_minimizer'], reports] },
            // The coverage denominator lives in the database directory; a
            // tarballed database has not been unpacked to a path we can name.
            kraken2_db.endsWith('.tar.gz') || kraken2_db.endsWith('.tgz')
                ? []
                : file("${kraken2_db}/inspect.txt").exists() ? file("${kraken2_db}/inspect.txt") : [],
            minimizer_thresholds.min_reads,
            minimizer_thresholds.min_distinct,
            minimizer_thresholds.max_duplication,
            minimizer_thresholds.min_coverage,
            minimizer_thresholds.distinct_scale,
        )
        ch_minimizer_drop = ABUNDANCE_MINIMIZER.out.drop_list.map { _meta, list -> list }
        ch_multiqc_files = ch_multiqc_files.mix(ABUNDANCE_MINIMIZER.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // MODULE: Which taxa are host leakage wearing a species name.
    //
    // Alignment-based depletion is not exhaustive, and what it misses is a
    // genuine read with genuine k-mers - invisible to every filter that judges
    // a taxon by its counts. Kraken2's read-level output says which taxon each
    // run of k-mers in a read was assigned to, so a read carrying host k-mers
    // can be recognised whatever the read as a whole was called.
    //
    // Scanned per sample because that file is one line per read; merged over
    // the cohort because "are this taxon's reads mostly host" is a question
    // about the taxon, not about one library.
    //
    if (host_kmer_filter) {
        HOSTKMER_SCAN(ch_classifiedreads, host_kmer_settings.taxid)
        HOSTKMER_MERGE(
            HOSTKMER_SCAN.out.table
                .map { _meta, table -> table }
                .collect(sort: true)
                .map { tables -> [[id: 'kraken2_host_kmer'], tables] },
            host_kmer_settings.taxid,
            host_kmer_settings.max_fraction,
            host_kmer_settings.min_reads,
        )
        ch_hostkmer_drop = HOSTKMER_MERGE.out.drop_list.map { _meta, list -> list }
        ch_hostkmer_evidence = HOSTKMER_MERGE.out.evidence
        ch_multiqc_files = ch_multiqc_files.mix(HOSTKMER_MERGE.out.mqc.map { _meta, mqc -> mqc })
    }

    // Union of the two lists: each condemns a taxon for its own reason, and
    // surviving one is no argument against the other. `collect` gives the
    //
    // MODULE: Bracken, and the combined table.
    //
    // Combined BEFORE the filters run rather than after, because decontam is
    // scored on that table and its verdict has to reach both filters. Bracken's
    // flat species-level counts are the right input for it: a combined kreport
    // is a hierarchy, so a clade's count already contains its children's and a
    // prevalence test on it would count the same reads at every rank.
    //
    def ch_bracken = channel.empty()
    def ch_bracken_report = channel.empty()
    def ch_bracken_combined = channel.empty()

    if (!skip_bracken) {
        BRACKEN_BRACKEN(ch_report_classified, ch_bracken_db)
        ch_bracken = BRACKEN_BRACKEN.out.reports
        ch_bracken_report = BRACKEN_BRACKEN.out.txt

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

    //
    // MODULE: Which taxa came out of the kit rather than the sample.
    //
    // The one filter here given an external measurement of the reagents, and so
    // the only one that can tell a contaminant from a rare organism at all. It
    // needs blanks (or DNA concentrations) that most public datasets do not
    // have, which is why it is opt-in and why it fails loudly rather than
    // quietly reporting nothing when they are missing.
    //
    def ch_decontam_drop = channel.empty()
    def ch_decontam_evidence = channel.empty()

    if (decontam_settings && !skip_bracken) {
        DECONTAM_FILTER(
            ch_bracken_combined,
            file(decontam_settings.metadata, checkIfExists: true),
            decontam_settings.neg_column ?: '',
            decontam_settings.neg_value,
            decontam_settings.conc_column ?: '',
            decontam_settings.method,
            decontam_settings.threshold,
            decontam_settings.batch_column ?: '',
            decontam_settings.batch_combine,
        )
        ch_decontam_drop = DECONTAM_FILTER.out.drop_list.map { _meta, list -> list }
        ch_decontam_evidence = DECONTAM_FILTER.out.evidence
        ch_multiqc_files = ch_multiqc_files.mix(DECONTAM_FILTER.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // MODULE: How much of this taxon arrives without a sample?
    //
    // The only step here that rewrites the table rather than naming taxids,
    // because its verdict is per LIBRARY: a taxon can be signal in one and
    // carryover in the next, and index hopping puts genuine reads of a genuine
    // organism into the wrong library, where nothing about the reads is wrong
    // for any evidence filter to find. It therefore runs BEFORE the abundance
    // filter, so the fractions that survive are recomputed against what is
    // left; the taxa it condemns everywhere still leave by the usual list.
    //
    // Exact on a Bracken table, which is flat and holds each read once. On the
    // combine_kreports hierarchy a zeroed clade leaves its ancestors' clade
    // counts stale - the same limitation drop_taxa() already documents, and
    // the same answer: docs/output.md says to take Bracken downstream.
    //
    //
    // MODULE: What the assembly says, which is the one thing read counts cannot.
    //
    // Libraries are POOLED before assembling. Per-library assembly is close to
    // useless on the cohorts this filter is for: on the CSI-Microbes plate the
    // false positives sit at 11-78 reads per well, and nothing assembles from
    // that. Pooling by group puts every library that should share an organism
    // into one assembly, which is where a real taxon finally has the depth to
    // produce a contig and a hopped read still does not.
    //
    // MEGAHIT rather than an RNA assembler even though the input is usually
    // RNA-seq - see modules/local/megahit for why the microbial fraction is the
    // wrong shape for an isoform model.
    //
    //
    // MODULE: the host's relatives, which are host reads under another name.
    //
    // Scored on the per-sample reports rather than the combined table because
    // the LINEAGE is what is needed, and the combined table flattens it away.
    //
    def ch_hostclade_drop = channel.empty()
    def ch_hostclade_evidence = channel.empty()

    if (host_clade) {
        HOSTCLADE_EXPAND(
            ch_report_classified
                .map { _meta, report -> report }
                .collect(sort: true)
                .map { reports -> [[id: 'kraken2_host_clade'], reports] },
            host_clade.taxid,
            host_clade.rank,
        )
        ch_hostclade_drop = HOSTCLADE_EXPAND.out.drop_list.map { _meta, list -> list }
        ch_hostclade_evidence = HOSTCLADE_EXPAND.out.evidence
        ch_multiqc_files = ch_multiqc_files.mix(HOSTCLADE_EXPAND.out.mqc.map { _meta, mqc -> mqc })
    }

    def ch_contig_drop = channel.empty()
    def ch_contig_evidence = channel.empty()
    def ch_contigs = channel.empty()

    if (assembly_settings) {
        def ch_pools = ch_reads
            .map { meta, reads ->
                def files = reads instanceof List ? reads : [reads]
                def key = assembly_settings.pool == 'sample'
                    ? meta.id
                    : assembly_settings.pool == 'group'
                        ? (assembly_settings.groups[meta.id] ?: 'ungrouped')
                        : 'all'
                [[id: "pool_${key}", single_end: meta.single_end], files]
            }
            .groupTuple()
            .map { meta, filesets ->
                [
                    meta,
                    filesets.collect { files -> files[0] },
                    meta.single_end ? [] : filesets.collect { files -> files[1] },
                ]
            }

        if (assembly_settings.assembler == 'metaspades') {
            SPADES_META(ch_pools)
            ch_contigs = SPADES_META.out.contigs
        }
        else {
            MEGAHIT(ch_pools)
            ch_contigs = MEGAHIT.out.contigs
        }

        // single_end on the way OUT, whatever the libraries were: a contig file
        // is one unpaired FASTA, and both Kraken2 modules read meta.single_end
        // to decide whether to pass --paired. Leaving a paired-end pool's meta
        // alone would hand Kraken2 one file and tell it there are two.
        ch_contigs = ch_contigs.map { meta, contigs -> [meta + [single_end: true], contigs] }

        // Kraken2 reads FASTA as happily as FASTQ, and the contigs are gzipped,
        // which is what the module's hardcoded --gzip-compressed needs. The
        // per-contig assignment is the point: the report alone gives counts, and
        // contig LENGTHS are what separate a 2 kb genome fragment from a 210 bp
        // fragment that happens to classify.
        def ch_contig_report = channel.empty()
        def ch_contig_assignments = channel.empty()

        if (use_daemon) {
            KRAKEN2_DAEMON_CONTIGS(
                ch_contigs,
                file(kraken2_db, checkIfExists: true).toAbsolutePath().toString(),
                false,
                true,
            )
            ch_contig_report = KRAKEN2_DAEMON_CONTIGS.out.report
            ch_contig_assignments = KRAKEN2_DAEMON_CONTIGS.out.classified_reads_assignment
        }
        else {
            KRAKEN2_CONTIGS(ch_contigs, ch_kraken2_db, false, true)
            ch_contig_report = KRAKEN2_CONTIGS.out.report
            ch_contig_assignments = KRAKEN2_CONTIGS.out.classified_reads_assignment
        }

        // Scored on the UNFILTERED combined kreport, like every other evidence
        // filter: its verdict is one of the inputs to the filtering, not a
        // comment on the result of it. The kreport and not Bracken, because a
        // contig assigned at genus has no species row to be scored against.
        CONTIG_EVIDENCE(
            KRAKENTOOLS_COMBINEKREPORTS.out.txt,
            ch_contig_assignments.map { _meta, table -> table }.collect(sort: true),
            ch_contig_report.map { _meta, report -> report }.collect(sort: true),
            assembly_settings.min_reads,
            assembly_settings.min_length,
            assembly_settings.min_contigs,
        )
        ch_contig_evidence = CONTIG_EVIDENCE.out.evidence
        ch_multiqc_files = ch_multiqc_files.mix(CONTIG_EVIDENCE.out.mqc.map { _meta, mqc -> mqc })

        // Evidence by default, a filter only on request. A contig is strong
        // positive evidence; its absence is weak negative evidence, because a
        // genuinely rare organism does not assemble either.
        if (assembly_settings.filter) {
            ch_contig_drop = CONTIG_EVIDENCE.out.drop_list.map { _meta, list -> list }
        }
    }

    def ch_control_drop = channel.empty()
    def ch_control_cells = channel.empty()
    def ch_control_evidence = channel.empty()
    def ch_kraken2_for_filter = KRAKENTOOLS_COMBINEKREPORTS.out.txt
    def ch_bracken_for_filter = ch_bracken_combined

    if (control_settings) {
        CONTROL_KRAKEN2(
            KRAKENTOOLS_COMBINEKREPORTS.out.txt,
            control_settings.controls ?: '',
            control_settings.ratio,
            control_settings.statistic,
            control_settings.min_reads,
            control_settings.floor_reads,
            control_settings.prevalence,
            control_settings.prevalence_min_reads,
        )
        ch_kraken2_for_filter = CONTROL_KRAKEN2.out.filtered
        ch_multiqc_files = ch_multiqc_files.mix(CONTROL_KRAKEN2.out.mqc.map { _meta, mqc -> mqc })

        // The per-library verdicts come from the KREPORT, never from Bracken,
        // and the difference is not cosmetic. Bracken's table is species-only;
        // the cell-by-taxon matrix carries whatever rank Kraken2 assigned each
        // read, and a great deal of it is genus. Measured on the CSI-Microbes
        // plate, every false positive the filter exists to remove sat at genus
        // Fusobacterium (848), which is absent from the Bracken table - so
        // Bracken-derived verdicts reached none of them and the specificity
        // was unchanged. The kreport holds both ranks and reaches both.
        ch_control_cells = CONTROL_KRAKEN2.out.drop_cells.map { _meta, list -> list }

        // Only Bracken's verdict joins the drop list. Both tables are scored,
        // but a taxid condemned twice would be no more removed than once, and
        // the flat table is the one whose "failed in every library" is a
        // statement about reads rather than about a clade sum.
        if (!skip_bracken) {
            CONTROL_BRACKEN(
                ch_bracken_combined,
                control_settings.controls ?: '',
                control_settings.ratio,
                control_settings.statistic,
                control_settings.min_reads,
                control_settings.floor_reads,
                control_settings.prevalence,
                control_settings.prevalence_min_reads,
            )
            ch_bracken_for_filter = CONTROL_BRACKEN.out.filtered
            // The whole-taxon verdict still comes from Bracken: that table is
            // flat and holds each read once, so "failed in every library" is a
            // statement about reads rather than about a clade sum.
            ch_control_drop = CONTROL_BRACKEN.out.drop_list.map { _meta, list -> list }
            ch_control_evidence = CONTROL_BRACKEN.out.evidence
        }
        else {
            ch_control_drop = CONTROL_KRAKEN2.out.drop_list.map { _meta, list -> list }
            ch_control_evidence = CONTROL_KRAKEN2.out.evidence
        }
    }

    // Union of every list: each condemns a taxon for its own reason, and
    // surviving one is no argument against the others. `collect` gives the
    // filter a single list, and emits an empty one when no filter ran.
    def ch_drop_list = minimizer_filter || host_kmer_filter || shuffle_settings || control_settings || host_clade || (assembly_settings && assembly_settings.filter) || (decontam_settings && !skip_bracken)
        ? ch_minimizer_drop.mix(ch_hostkmer_drop).mix(ch_decontam_drop).mix(ch_shuffle_drop).mix(ch_control_drop).mix(ch_contig_drop).mix(ch_hostclade_drop).collect(sort: true)
        : channel.value([])

    FILTER_KRAKEN2(ch_kraken2_for_filter, ch_drop_list, min_rel_abundance, min_samples, min_reads)

    def ch_bracken_combined_filtered = channel.empty()

    if (!skip_bracken) {
        FILTER_BRACKEN(ch_bracken_for_filter, ch_drop_list, min_rel_abundance, min_samples, min_reads)
        ch_bracken_combined_filtered = FILTER_BRACKEN.out.filtered
    }

    def ch_krona = channel.empty()

    if (!skip_krona) {
        //
        // Krona is rendered from the Bracken-corrected report when Bracken ran,
        // because that is the abundance estimate users are meant to interpret.
        //
        // Krona copes with an all-unclassified report, so it keeps every sample.
        def ch_for_krona = skip_bracken ? ch_kraken2_report : ch_bracken_report

        KRAKENTOOLS_KREPORT2KRONA(ch_for_krona)
        KRONA_KTIMPORTTEXT(KRAKENTOOLS_KREPORT2KRONA.out.txt)
        ch_krona = KRONA_KTIMPORTTEXT.out.html
    }

    //
    // MODULE: The same counts as a BIOM table, for tools outside this pipeline.
    //
    // Built from the per-sample reports Krona uses, for the same reason: when
    // Bracken ran, its redistributed report is the one that answers "how much",
    // and the raw Kraken2 report is only the right source when it did not.
    //
    // Restricted afterwards to the taxa in the FILTERED combined kreport, not
    // the Bracken one. Bracken's table is species-only, and kraken-biom emits a
    // row per rank; intersecting against a species-only list would silently
    // delete every genus, family and order observation - the same trap the
    // negative-control filter hit from the other direction.
    //
    def ch_biom = channel.empty()

    if (export_biom) {
        KRAKENBIOM(
            (skip_bracken ? ch_report_classified : ch_bracken_report)
                .toSortedList { entry_a, entry_b -> entry_a[0].id <=> entry_b[0].id }
                .filter { entries -> entries }
                .map { entries ->
                    [
                        [
                            id: skip_bracken ? 'kraken2_biom' : 'bracken_biom',
                            names: entries.collect { entry -> entry[0].id }.join(' '),
                        ],
                        entries.collect { entry -> entry[1] },
                    ]
                },
            FILTER_KRAKEN2.out.filtered.map { _meta, table -> table },
            biom_metadata ? file(biom_metadata, checkIfExists: true) : [],
        )
        ch_biom = KRAKENBIOM.out.biom
    }

    emit:
    report = ch_kraken2_report // channel: [ val(meta), path(report) ]
    report_combined = KRAKENTOOLS_COMBINEKREPORTS.out.txt // channel: [ val(meta), path(txt) ]
    bracken = ch_bracken // channel: [ val(meta), path(tsv) ]
    bracken_report = ch_bracken_report // channel: [ val(meta), path(txt) ] - Bracken's `-w` kreport
    bracken_combined = ch_bracken_combined // channel: [ val(meta), path(txt) ]
    report_combined_filtered = FILTER_KRAKEN2.out.filtered // channel: [ val(meta), path(tsv) ]
    bracken_combined_filtered = ch_bracken_combined_filtered // channel: [ val(meta), path(tsv) ]
    krona = ch_krona // channel: [ val(meta), path(html) ]
    biom = ch_biom // channel: [ val(meta), path(biom) ]
    contigs = ch_contigs // channel: [ val(meta), path(fa.gz) ]
    contig_evidence = ch_contig_evidence // channel: [ val(meta), path(tsv) ]
    host_clade_evidence = ch_hostclade_evidence // channel: [ val(meta), path(tsv) ]
    drop_list = ch_drop_list // channel: [ path(txt) ] - every filter's verdict, for consumers outside this subworkflow
    minimizer_evidence = minimizer_filter ? ABUNDANCE_MINIMIZER.out.evidence : channel.empty() // channel: [ val(meta), path(tsv) ]
    host_kmer_evidence = ch_hostkmer_evidence // channel: [ val(meta), path(tsv) ]
    decontam_evidence = ch_decontam_evidence // channel: [ val(meta), path(tsv) ]
    shuffle_evidence = ch_shuffle_evidence // channel: [ val(meta), path(tsv) ]
    control_evidence = ch_control_evidence // channel: [ val(meta), path(tsv) ]
    // Per-library verdicts, for the objects the combined tables do not cover -
    // the cell-by-taxon matrix is built from the per-read assignments, so the
    // cells zeroed here would otherwise survive there.
    control_cells = ch_control_cells.collect(sort: true).ifEmpty([]) // channel: [ path(tsv) ]
    classifiedreads = ch_classifiedreads // channel: [ val(meta), path(txt) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

//
// True if a Kraken-style report contains at least one taxon, i.e. anything
// beyond the unclassified row. Keyed on the TAXID rather than the rank code so
// it holds for both layouts this pipeline can produce: Kraken2 ends
// `... rank taxid name` and spells ranks as codes, KrakenUniq ends
// `... taxid rank name` and spells them out. Taxid 0 is unclassified in both,
// so anything else with reads means something was classified.
//
// A report is one line per taxon - a few MB even against a large database -
// and `any` stops at the first classified row.
//
def hasClassifiedReads(report) {
    return report.readLines().any { line ->
        if (line.startsWith('#') || line.startsWith('%')) {
            return false
        }
        def fields = line.split('\t')
        if (fields.size() < 5) {
            return false
        }
        def clade_reads = fields[1].trim().isInteger() ? fields[1].trim() as Integer : null
        if (!clade_reads) {
            return false
        }
        def tail = fields[-2].trim()
        def taxid = tail.isInteger() ? tail as Integer : (fields[-3].trim().isInteger() ? fields[-3].trim() as Integer : null)
        return taxid != null && taxid != 0
    }
}

//
// Resolve ENA/SRA accessions to their runs, download every run in parallel and
// (optionally) merge the runs that belong to the same experiment or sample.
//

include { FASTQDL_METADATA } from '../../../modules/local/fastqdl/metadata/main'
include { FASTQDL          } from '../../../modules/nf-core/fastqdl/main'
include { CAT_FASTQ        } from '../../../modules/nf-core/cat/fastq/main'

workflow FASTQ_DOWNLOAD_FASTQDL {

    take:
    ch_accessions // channel: [ val(meta), val(accession) ]
    group_runs_by //  string: 'run', 'experiment' or 'sample'

    main:

    //
    // Turn each (possibly umbrella) accession into its list of runs. Doing this
    // up front is what lets a 200-run BioProject be downloaded by 200 parallel
    // tasks rather than by one serial `fastq-dl` invocation.
    //
    FASTQDL_METADATA(ch_accessions)

    def ch_runs = FASTQDL_METADATA.out.runsheet
        .map { _meta, runsheet -> runsheet }
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                id: row.run_accession,
                run_accession: row.run_accession,
                experiment_accession: row.experiment_accession,
                sample_accession: row.sample_accession,
                study_accession: row.study_accession,
                title: row.sample_title,
                organism: row.scientific_name,
                // Provisional; re-derived below from what actually downloaded.
                single_end: (row.library_layout ?: '').toUpperCase() != 'PAIRED',
            ]
            [meta, row.run_accession]
        }
        // Overlapping queries (e.g. a BioProject plus one of its experiments)
        // must not download the same run twice.
        .unique { meta, _accession -> meta.run_accession }

    //
    // MODULE: Download one run per task
    //
    FASTQDL(ch_runs)

    //
    // Trust the files on disk over the archive's `library_layout` field: runs
    // are regularly mislabelled, and ENA additionally serves an unpaired
    // "orphan" FASTQ alongside _1/_2 for some submissions.
    //
    def ch_runs_reads = FASTQDL.out.fastq.map { meta, fastq ->
        def files = (fastq instanceof List ? fastq : [fastq]).sort { fq -> fq.name }
        def read1 = files.find { fq -> fq.name.endsWith('_1.fastq.gz') }
        def read2 = files.find { fq -> fq.name.endsWith('_2.fastq.gz') }
        read1 && read2
            ? [meta + [single_end: false], [read1, read2]]
            : [meta + [single_end: true], [files.first()]]
    }

    //
    // Group runs that belong together. `run` keeps every run as its own sample.
    //
    def ch_grouped = ch_runs_reads
        .map { meta, reads ->
            def key = group_runs_by == 'run'
                ? meta.run_accession
                : group_runs_by == 'sample' ? meta.sample_accession : meta.experiment_accession
            [key, meta, reads]
        }
        .groupTuple(by: 0)
        .map { key, metas, reads ->
            // Sort by run accession so that merged FASTQs are byte-identical
            // between resumes and between machines.
            def ordered = [metas, reads].transpose().sort { entry -> entry[0].run_accession }
            def ordered_metas = ordered.collect { entry -> entry[0] }
            def ordered_reads = ordered.collect { entry -> entry[1] }
            if (ordered_metas.collect { meta -> meta.single_end }.unique().size() > 1) {
                error("Cannot merge runs of '${key}': they are a mix of single- and paired-end (runs: ${ordered_metas.collect { meta -> meta.run_accession }.join(', ')}). Re-run with --group_runs_by run.")
            }
            def meta = ordered_metas.first() + [id: key, n_runs: ordered_metas.size()]
            [meta, ordered_reads]
        }
        .branch { _meta, reads ->
            single: reads.size() == 1
            merge: true
        }

    //
    // MODULE: Concatenate the FASTQs of multi-run experiments/samples
    //
    CAT_FASTQ(ch_grouped.merge.map { meta, reads -> [meta, reads.flatten()] })

    def ch_reads = ch_grouped.single
        .map { meta, reads -> [meta, reads.first()] }
        .mix(CAT_FASTQ.out.reads.map { meta, reads -> [meta, reads instanceof List ? reads : [reads]] })

    emit:
    reads = ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    runinfo = FASTQDL_METADATA.out.runinfo // channel: [ val(meta), path(tsv) ]
    runsheet = FASTQDL_METADATA.out.runsheet // channel: [ val(meta), path(csv) ]
}

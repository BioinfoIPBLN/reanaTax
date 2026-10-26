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
                // 'dna', 'rna' or null: decides the host aligner under
                // --host_aligner auto (Bowtie2 for DNA, HISAT2 otherwise).
                molecule: libraryMolecule(row.library_strategy, row.library_source),
                // ILLUMINA, OXFORD_NANOPORE, PACBIO_SMRT...: checked against --long_reads.
                platform: row.instrument_platform,
                // Provisional; re-derived below from what actually downloaded.
                single_end: (row.library_layout ?: '').toUpperCase() != 'PAIRED',
            ]
            [meta, row.run_accession]
        }
        // Overlapping queries (e.g. a BioProject plus one of its experiments)
        // must not download the same run twice.
        .unique { meta, _accession -> meta.run_accession }

    //
    // Platform audit before download. --long_reads is a whole-run mode (all
    // short or all long). Each run's platform is already on the run sheet here
    // (metadata only - no bytes), so a run that cannot analyse as requested is
    // rejected before its FASTQ is fetched. longReadLibrary is the later,
    // redundant guard.
    def ch_runs_validated = ch_runs.map { meta, accession ->
        def platform = (meta.platform ?: '').toString().toUpperCase()
        def archive_kind = [OXFORD_NANOPORE: 'ont', PACBIO_SMRT: 'pacbio'][platform]
        if (!params.long_reads) {
            if (archive_kind) {
                error("Run ${meta.run_accession ?: meta.id} was sequenced on ${platform}. Long reads take their own route: add --long_reads ${archive_kind}, and keep short-read runs to a separate run of the pipeline.")
            }
        } else if (platform && archive_kind != params.long_reads) {
            error("--long_reads ${params.long_reads}: run ${meta.run_accession ?: meta.id} was sequenced on ${platform}. Keep each platform to its own run of the pipeline.")
        }
        [meta, accession]
    }

    //
    // MODULE: Download one run per task
    //
    FASTQDL(ch_runs_validated)

    //
    // Trust the files on disk over the archive's `library_layout` field: runs
    // are regularly mislabelled, and ENA additionally serves an unpaired
    // "orphan" FASTQ alongside _1/_2 for some submissions.
    //
    def ch_runs_reads = FASTQDL.out.fastq.map { meta, fastq ->
        def files = (fastq instanceof List ? fastq : [fastq]).toSorted { fq -> fq.name }
        def read1 = files.find { fq -> fq.name.endsWith('_1.fastq.gz') }
        def read2 = files.find { fq -> fq.name.endsWith('_2.fastq.gz') }
        read1 && read2
            ? [meta + [single_end: false], [read1, read2]]
            : [meta + [single_end: true], [files.first()]]
    }

    //
    // Group runs that belong together. `run` keeps every run as its own sample.
    //
    // Each group is released as soon as its own runs have downloaded. A plain
    // groupTuple only emits when the channel closes, so trimming, alignment and
    // classification of every sample waited for the slowest download in the
    // cohort. The run sheet already says how many runs each group has, long
    // before the downloads finish, so groupKey can carry that count. An
    // incomplete group is still emitted at the end (remainder), as before.
    //
    def ch_expected_runs = ch_runs
        .map { meta, _accession -> [runGroup(meta, group_runs_by), meta.run_accession] }
        .groupTuple(by: 0)
        .map { key, runs -> [key, runs.size()] }

    def ch_grouped = ch_runs_reads
        .map { meta, reads -> [runGroup(meta, group_runs_by), meta, reads] }
        .combine(ch_expected_runs, by: 0)
        .map { key, meta, reads, expected -> [groupKey(key, expected), meta, reads] }
        .groupTuple(by: 0, remainder: true)
        .map { group, metas, reads ->
            def key = group.getGroupTarget()
            // Sort by run accession so that merged FASTQs are byte-identical
            // between resumes and between machines.
            def ordered = [metas, reads].transpose().sort { entry -> entry[0].run_accession }
            def ordered_metas = ordered.collect { entry -> entry[0] }
            def ordered_reads = ordered.collect { entry -> entry[1] }
            if (ordered_metas.collect { meta -> meta.single_end }.unique().size() > 1) {
                error("Cannot merge runs of '${key}': they are a mix of single- and paired-end (runs: ${ordered_metas.collect { meta -> meta.run_accession }.join(', ')}). Re-run with --group_runs_by run.")
            }
            if (ordered_metas.collect { meta -> meta.molecule }.unique().size() > 1) {
                error("Cannot merge runs of '${key}': they are a mix of DNA and RNA libraries (runs: ${ordered_metas.collect { meta -> "${meta.run_accession} ${meta.molecule ?: 'unknown'}" }.join(', ')}), which are depleted with different aligners. Re-run with --group_runs_by run.")
            }
            if (ordered_metas.collect { meta -> meta.platform }.unique().size() > 1) {
                error("Cannot merge runs of '${key}': they come from different platforms (runs: ${ordered_metas.collect { meta -> "${meta.run_accession} ${meta.platform ?: 'unknown'}" }.join(', ')}). Re-run with --group_runs_by run.")
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

// Whether a run's library is DNA or RNA, from the archive's own labels: 'dna',
// 'rna', or null when neither says. The library strategy decides first, because
// the source is the less reliable of the two - RNA-seq runs submitted with
// source GENOMIC are common (20 of PRISM's CLID set). The source decides only
// when the strategy is OTHER or missing, as for libraries submitted as
// GENOMIC/OTHER.
def libraryMolecule(strategy, source) {
    def rna_strategies = ['RNA-SEQ', 'SSRNA-SEQ', 'SNRNA-SEQ', 'SCRNA-SEQ', 'MIRNA-SEQ', 'NCRNA-SEQ',
                          'FL-CDNA', 'EST', 'RIBO-SEQ', 'RIP-SEQ', 'MERIP-SEQ']
    def dna_strategies = ['WGS', 'WGA', 'WXS', 'WCS', 'AMPLICON', 'TARGETED-CAPTURE', 'CLONE', 'POOLCLONE',
                          'CLONEEND', 'FINISHING', 'CHIP-SEQ', 'ATAC-SEQ', 'MNASE-SEQ', 'DNASE-HYPERSENSITIVITY',
                          'BISULFITE-SEQ', 'MRE-SEQ', 'MEDIP-SEQ', 'MBD-SEQ', 'HI-C', 'RAD-SEQ', 'GBS', 'TN-SEQ',
                          'FAIRE-SEQ', 'CHIA-PET', 'NOME-SEQ', 'CHM-SEQ', 'SYNTHETIC-LONG-READ']
    def strategy_key = (strategy ?: '').toString().trim().toUpperCase()
    if (strategy_key in rna_strategies) {
        return 'rna'
    }
    if (strategy_key in dna_strategies) {
        return 'dna'
    }
    def source_key = (source ?: '').toString().trim().toUpperCase()
    if (source_key in ['TRANSCRIPTOMIC', 'TRANSCRIPTOMIC SINGLE CELL', 'METATRANSCRIPTOMIC', 'VIRAL RNA']) {
        return 'rna'
    }
    if (source_key in ['GENOMIC', 'GENOMIC SINGLE CELL', 'METAGENOMIC']) {
        return 'dna'
    }
    return null
}

// The sample a run belongs to under --group_runs_by.
def runGroup(meta, group_runs_by) {
    return group_runs_by == 'run'
        ? meta.run_accession
        : group_runs_by == 'sample' ? meta.sample_accession : meta.experiment_accession
}

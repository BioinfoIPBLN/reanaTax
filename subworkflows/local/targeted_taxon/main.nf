//
// One organism from the classification, followed all the way to differential
// expression of ITS OWN genes.
//
// The rest of the pipeline answers "who is there, and how much". This answers
// "and what were they doing", for one taxon chosen after the fact. Kraken2
// already recorded which reads it assigned where, so the reads for a taxon can
// be pulled back out, aligned to that organism's genome, counted against its
// annotation, and tested between the same two groups as everything else.
//
// Almost nothing here is new. Extraction is the only new step; the reference
// preparation, alignment, sorting and gene counting are the same subworkflow
// the host pass uses, and the merge-and-test tail is the same one the host gene
// counts go through. They are aliased rather than copied so there is one
// implementation of each to keep correct.
//
// A caveat worth stating plainly: reads reach this branch because a k-mer
// classifier assigned them, and a classifier's false positives arrive looking
// exactly like true positives. The alignment rate in the published HISAT2
// summary is the check - reads that genuinely came from the target organism
// align to its genome, and reads Kraken2 mis-assigned do not.
//

include { KRAKENTOOLS_EXTRACTREADS                       } from '../../../modules/local/krakentools/extractreads/main'
include { PREPARE_HOST_REFERENCE as PREPARE_TARGET_REFERENCE } from '../prepare_host_reference'
include { HOST_DEPLETION_HISAT2 as TARGET_ALIGNMENT      } from '../host_depletion_hisat2'
include { HOST_EXPRESSION as TARGET_EXPRESSION           } from '../host_expression'
include { QUANTIFY_KALLISTO as TARGET_KALLISTO           } from '../quantify_kallisto'
include { HOSTCOUNTS_MERGE as TARGET_COUNTS_MERGE        } from '../../../modules/local/hostcounts/merge/main'

workflow TARGETED_TAXON {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ] - the non-host fraction
    ch_classifiedreads // channel: [ val(meta), path(txt) ] - Kraken2 per-read output
    ch_report // channel: [ val(meta), path(report) ] - Kraken2 report, for --include-children
    taxid // string: taxon to extract, comma-separated for several
    reference // map: [kind, value] naming the target genome, from classifyHostReference()
    gtf // string: annotation for that genome, or null to align only
    transcripts // string: its transcriptome FASTA - quantify with kallisto instead of aligning
    tx2gene // string: two-column transcript->gene TSV for the above, or null
    de_methods // list: subset of ['deseq2', 'edger'], empty to stop after counting
    metadata // string: --da_metadata, or null
    ncbi_group // string: passed through to the reference download

    main:

    def ch_multiqc_files = channel.empty()

    //
    // MODULE: the reads Kraken2 put on this taxon.
    //
    // Joined rather than combined: the three inputs are per sample and must line
    // up sample by sample. `join` keys on the meta and drops anything unmatched,
    // which is the behaviour wanted - a sample whose classification is missing
    // must not be silently paired with another sample's reads.
    //
    KRAKENTOOLS_EXTRACTREADS(
        ch_reads.join(ch_classifiedreads).join(ch_report),
        taxid,
    )
    ch_multiqc_files = ch_multiqc_files.mix(KRAKENTOOLS_EXTRACTREADS.out.log.map { _meta, log_file -> log_file })

    //
    // Two routes to per-gene numbers for this taxon. With --target_transcripts
    // the extracted reads are pseudoaligned straight to the transcriptome and
    // no genome, index or GTF is needed at all; otherwise they are aligned and
    // counted exactly as host reads are.
    //
    def ch_counts = channel.empty()
    def ch_bam = channel.empty()
    def ch_summary = channel.empty()

    if (transcripts) {
        TARGET_KALLISTO(KRAKENTOOLS_EXTRACTREADS.out.reads, transcripts, tx2gene)
        ch_counts = TARGET_KALLISTO.out.counts
        ch_multiqc_files = ch_multiqc_files.mix(TARGET_KALLISTO.out.multiqc_files)
    }
    else {

    //
    // SUBWORKFLOW: the target's genome, prepared exactly as a host genome is.
    //
    PREPARE_TARGET_REFERENCE(
        reference.kind == 'fasta' ? reference.value : null,
        reference.kind == 'index' ? reference.value : null,
        reference.kind == 'accession' ? reference.value : null,
        reference.kind == 'taxid' ? reference.value : null,
        ncbi_group,
        gtf,
        'hisat2',
        null,
    )

    //
    // SUBWORKFLOW: align and count.
    //
    // The host-depletion subworkflow is reused whole. It splits aligned from
    // unaligned and counts the aligned side against the GTF, which is what is
    // wanted here - only the interesting half is the opposite one. Its
    // `reads` output (what did NOT align to the target) is deliberately dropped.
    //
    // Qualimap is skipped: its output directory is named from the sample and
    // would collide with the host pass's report.
    //
    TARGET_ALIGNMENT(
        KRAKENTOOLS_EXTRACTREADS.out.reads,
        PREPARE_TARGET_REFERENCE.out.index,
        PREPARE_TARGET_REFERENCE.out.fasta,
        true,
        true,
        null,
        gtf,
        true,
        0,
    )
    ch_multiqc_files = ch_multiqc_files.mix(TARGET_ALIGNMENT.out.multiqc_files)
    ch_bam = TARGET_ALIGNMENT.out.bam
    ch_summary = TARGET_ALIGNMENT.out.summary

    if (gtf) {
        TARGET_COUNTS_MERGE(
            TARGET_ALIGNMENT.out.host_counts
                .map { _meta, table -> table }
                .collect(sort: true)
                .map { tables -> [[id: 'reanatax'], tables] }
        )
        ch_counts = TARGET_COUNTS_MERGE.out.counts
    }

    }

    //
    // SUBWORKFLOW: test the counts, whichever route produced them.
    //
    // Same models as the host side, for the same reason: these are gene counts,
    // not a composition. Correlation is off - it asks about host versus microbe,
    // and both sides here are the microbe.
    //
    def ch_results = channel.empty()

    if (de_methods && metadata) {
        TARGET_EXPRESSION(
            ch_counts,
            channel.empty(),
            metadata,
            de_methods,
            false,
        )
        ch_results = TARGET_EXPRESSION.out.results
        ch_multiqc_files = ch_multiqc_files.mix(TARGET_EXPRESSION.out.multiqc_files)
    }

    emit:
    reads = KRAKENTOOLS_EXTRACTREADS.out.reads // channel: [ val(meta), [ path(fastq.gz) ] ]
    bam = ch_bam // channel: [ val(meta), path(bam) ] - empty on the kallisto route
    summary = ch_summary // channel: [ val(meta), path(log) ]
    counts = ch_counts // channel: [ val(meta), path(tsv) ]
    results = ch_results // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

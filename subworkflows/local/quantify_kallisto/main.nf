//
// Pseudoalignment quantification against a transcriptome.
//
// An alternative to align-then-count wherever the pipeline QUANTIFIES: host
// gene expression, and the targeted taxon. It is not offered as an alternative
// for host DEPLETION, and that is a deliberate limit rather than an omission -
// a transcriptome index can only remove reads that came from an annotated
// transcript, so intronic, intergenic, repeat and satellite reads would all
// survive depletion and go on to be classified as spurious microbes. Those are
// exactly the reads the second host assembly exists to catch.
//
// The output is the same counts matrix HOSTCOUNTS_MERGE produces, so the
// differential-expression tail is shared and there is one of it to keep right.
//
include { KALLISTO_INDEX } from '../../../modules/local/kallisto/index/main'
include { KALLISTO_QUANT } from '../../../modules/local/kallisto/quant/main'
include { KALLISTO_MERGE } from '../../../modules/local/kallisto/merge/main'

workflow QUANTIFY_KALLISTO {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    transcripts // string: transcriptome FASTA
    tx2gene // string: two-column transcript->gene TSV, or null to stay at transcript level

    main:

    KALLISTO_INDEX(file(transcripts, checkIfExists: true))

    KALLISTO_QUANT(ch_reads, KALLISTO_INDEX.out.index)

    //
    // Collected with a fixed sort so the column order of the matrix - and with
    // it the task hash - is the same on a resumed run as on the original.
    //
    KALLISTO_MERGE(
        KALLISTO_QUANT.out.abundance
            .map { _meta, abundance -> abundance }
            .collect(sort: true)
            .map { files -> [[id: 'reanatax'], files] },
        tx2gene ? file(tx2gene, checkIfExists: true) : [],
    )

    emit:
    counts = KALLISTO_MERGE.out.counts // channel: [ val(meta), path(tsv) ]
    tpm = KALLISTO_MERGE.out.tpm // channel: [ val(meta), path(tsv) ]
    lengths = KALLISTO_MERGE.out.lengths // channel: [ val(meta), path(tsv) ]
    abundance = KALLISTO_QUANT.out.abundance // channel: [ val(meta), path(tsv) ]
    multiqc_files = KALLISTO_QUANT.out.run_info.map { _meta, info -> info } // channel: path(json)
}

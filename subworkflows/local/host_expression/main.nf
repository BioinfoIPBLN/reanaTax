//
// The host half of the split, treated as data rather than as a QC by-product.
//
// --quantify_host already counted host reads against the annotation, but until
// now those counts were published per sample and nothing read them. Three
// things happen here:
//
//   1. the caller merges the per-sample counts into one gene x sample matrix -
//      from featureCounts, or from kallisto, whichever quantified;
//   2. that matrix is tested for differential expression between the SAME two
//      groups the microbial side is tested on, so one --da_* design describes
//      both halves of the library;
//   3. optionally, host expression is correlated against microbial abundance -
//      Monteleone et al.'s third layer, and the only analysis here that needs
//      both halves at once.
//
// The models are negative binomial, not compositional. ALDEx2 and ANCOM-BC2
// exist because a microbial profile has no absolute scale; gene counts do, up to
// a library-size factor that DESeq2 and edgeR estimate directly. Running the
// microbial methods on gene counts would work mechanically and answer a
// different question badly.
//

include { GENE_DE_DESEQ2          } from '../../../modules/local/genede/deseq2/main'
include { GENE_DE_EDGER           } from '../../../modules/local/genede/edger/main'
include { HOSTMICROBE_CORRELATION } from '../../../modules/local/hostmicrobe/correlation/main'

workflow HOST_EXPRESSION {

    take:
    ch_counts // channel: [ val(meta), path(tsv) ] - a merged gene x sample matrix
    ch_microbial // channel: [ val(meta), path(combined_abundance) ] - or empty
    metadata // string: path to the sample metadata TSV
    methods // list: subset of ['deseq2', 'edger']
    correlate // boolean: also correlate host expression against microbial abundance

    main:

    def ch_multiqc_files = channel.empty()

    def ch_metadata = channel.value(file(metadata, checkIfExists: true))
    def ch_results = channel.empty()

    if (methods.contains('deseq2')) {
        GENE_DE_DESEQ2(ch_counts, ch_metadata)
        ch_results = ch_results.mix(GENE_DE_DESEQ2.out.results)
        ch_multiqc_files = ch_multiqc_files.mix(GENE_DE_DESEQ2.out.mqc.map { _meta, mqc -> mqc })
    }

    if (methods.contains('edger')) {
        GENE_DE_EDGER(ch_counts, ch_metadata)
        ch_results = ch_results.mix(GENE_DE_EDGER.out.results)
        ch_multiqc_files = ch_multiqc_files.mix(GENE_DE_EDGER.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // MODULE: host expression against microbial abundance.
    //
    // Joined on the meta rather than combined blindly: both sides are
    // cohort-level tables keyed [id: 'reanatax'], and a join makes that
    // assumption explicit and fails loudly if it ever stops holding.
    //
    def ch_correlation = channel.empty()

    if (correlate) {
        HOSTMICROBE_CORRELATION(ch_counts.join(ch_microbial))
        ch_correlation = HOSTMICROBE_CORRELATION.out.results
        ch_multiqc_files = ch_multiqc_files.mix(HOSTMICROBE_CORRELATION.out.mqc.map { _meta, mqc -> mqc })
    }

    emit:
    results = ch_results // channel: [ val(meta), path(tsv) ]
    correlation = ch_correlation // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files // channel: path(mqc)
}

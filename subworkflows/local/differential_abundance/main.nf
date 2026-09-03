//
// Differential abundance between two groups of samples.
//
// Group membership comes from --da_metadata, a TSV whose first column is the
// sample id and whose remaining columns are factors. It is deliberately a
// separate file rather than a samplesheet column: the accession and folder
// input routes have no samplesheet at all, and this way one metadata file
// serves all three.
//

include { DIFFERENTIAL_ABUNDANCE_ALDEX2  } from '../../../modules/local/differentialabundance/aldex2/main'
include { DIFFERENTIAL_ABUNDANCE_ANCOMBC } from '../../../modules/local/differentialabundance/ancombc/main'

workflow DIFFERENTIAL_ABUNDANCE {

    take:
    ch_counts // channel: [ val(meta), path(bracken_combined) ]
    metadata // string: path to the sample metadata TSV
    methods // list: subset of ['aldex2', 'ancombc']

    main:

    def ch_metadata = channel.value(file(metadata, checkIfExists: true))
    def ch_results = channel.empty()
    def ch_multiqc_files = channel.empty()

    if (methods.contains('aldex2')) {
        DIFFERENTIAL_ABUNDANCE_ALDEX2(ch_counts, ch_metadata)
        ch_results = ch_results.mix(DIFFERENTIAL_ABUNDANCE_ALDEX2.out.results)
        ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_ABUNDANCE_ALDEX2.out.mqc.map { _meta, mqc -> mqc })
    }

    if (methods.contains('ancombc')) {
        DIFFERENTIAL_ABUNDANCE_ANCOMBC(ch_counts, ch_metadata)
        ch_results = ch_results.mix(DIFFERENTIAL_ABUNDANCE_ANCOMBC.out.results)
        ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_ABUNDANCE_ANCOMBC.out.mqc.map { _meta, mqc -> mqc })
    }

    emit:
    results = ch_results // channel: [ val(meta), path(tsv) ]
    multiqc_files = ch_multiqc_files // channel: path(mqc)
}

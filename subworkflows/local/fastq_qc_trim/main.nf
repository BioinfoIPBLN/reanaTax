//
// Read QC before and after adapter/quality trimming.
//

include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'

workflow FASTQ_QC_TRIM {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    adapter_fasta // path(fasta) or [] - a plain value, not a channel
    skip_fastqc // boolean
    skip_trimming // boolean
    save_trimmed_fail // boolean

    main:

    def ch_multiqc_files = channel.empty()

    if (!skip_fastqc) {
        FASTQC_RAW(ch_reads)
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.map { _meta, zip -> zip })
    }

    def ch_trimmed = ch_reads
    def ch_fastp_json = channel.empty()

    if (!skip_trimming) {
        FASTP(
            ch_reads.map { meta, reads -> [meta, reads, adapter_fasta ?: []] },
            false,
            save_trimmed_fail,
            false,
        )
        ch_fastp_json = FASTP.out.json
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.map { _meta, json -> json })

        // fastp emits a single-element list for single-end input; normalise so
        // downstream processes always see a list.
        ch_trimmed = FASTP.out.reads.map { meta, reads -> [meta, reads instanceof List ? reads : [reads]] }

        if (!skip_fastqc) {
            FASTQC_TRIM(ch_trimmed)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIM.out.zip.map { _meta, zip -> zip })
        }
    }

    emit:
    reads = ch_trimmed // channel: [ val(meta), [ path(fastq) ] ]
    fastp_json = ch_fastp_json // channel: [ val(meta), path(json) ]
    multiqc_files = ch_multiqc_files // channel: path(file)
}

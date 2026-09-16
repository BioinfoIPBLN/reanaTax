//
// Barcode-aware host alignment for single-cell data, and the non-host reads
// that come out of it with their cell barcode attached.
//
// This REPLACES HOST_DEPLETION_HISAT2 on the single-cell route rather than
// supplementing it. In bulk, host depletion discards the host reads; in single
// cell they are the other half of the experiment - the cell-by-gene matrix -
// and the cell barcode is the only thing that joins them to the microbial
// reads. One barcode-aware aligner has to do both jobs in one pass, which is
// what STARsolo is for and why SAHMI uses it too.
//
// What comes out is an ordinary single-end FASTQ of non-host reads, so the rest
// of the pipeline - Kraken2, Bracken, Krona, the evidence filters, the
// exploreMetaTax bundle - runs on it unchanged. The barcode is carried in the
// read name and recovered afterwards by SCTAXA_COUNTS, so the cell-by-taxon
// matrix is built from the same classification the pseudobulk profile is.
//

include { STARSOLO          } from '../../../modules/local/starsolo/align/main'
include { STARSOLO_UNMAPPED } from '../../../modules/local/starsolo/unmapped/main'

workflow SINGLECELL_STARSOLO {

    take:
    ch_reads // channel: [ val(meta), [ barcode_fastq, cdna_fastq ] ]
    ch_index // channel: [ val(meta), path(star_index) ]
    ch_gtf // channel: [ val(meta), path(gtf) ]
    whitelist // string: barcode whitelist path, or 'None'

    main:

    def ch_whitelist = whitelist && whitelist.toString().toLowerCase() != 'none'
        ? channel.value(file(whitelist, checkIfExists: true))
        : channel.value([])

    // .first() on both, and it is not cosmetic. PREPARE_HOST_REFERENCE emits
    // QUEUE channels holding one item each, and a process pairs a queue channel
    // element-by-element against the reads - so without this STARsolo aligns
    // the FIRST sample and silently drops every other one. Measured on the
    // CSI-Microbes droplet cohort: 3 libraries in, "STARSOLO | 1 of 1", two
    // thirds of the cohort gone with no error. The bulk route already does this
    // at every HISAT2_ALIGN* call site; this path had been exercised only with
    // one sample at a time.
    STARSOLO(ch_reads, ch_index.first(), ch_gtf.first(), ch_whitelist)
    STARSOLO_UNMAPPED(STARSOLO.out.bam)

    //
    // Single-end from here on. STARsolo writes one BAM record per cDNA read;
    // the barcode read carries no biological sequence, so there is nothing to
    // classify in it and nothing is lost by dropping the pairing.
    //
    def ch_nonhost = STARSOLO_UNMAPPED.out.reads.map { meta, reads ->
        [meta + [single_end: true], reads]
    }

    emit:
    reads = ch_nonhost // channel: [ val(meta), path(fastq) ]
    solo = STARSOLO.out.solo // channel: [ val(meta), path(Solo.out) ]
    bam = STARSOLO.out.bam // channel: [ val(meta), path(bam) ]
    barcode_stats = STARSOLO_UNMAPPED.out.stats // channel: [ val(meta), path(tsv) ]
    multiqc_files = STARSOLO.out.log_final.map { _meta, log_file -> log_file } // channel: path(log)
}

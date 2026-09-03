//
// GATK PathSeq: alignment-based classification of the non-host fraction.
//
// A third route alongside Kraken2 and KrakenUniq, and the only one here that is
// not k-mer based. PathSeq aligns the surviving reads to a microbe reference
// with BWA and scores taxa from those alignments, so it fails differently from
// the k-mer classifiers: it is not fooled by a composition-driven chance match,
// and it can place a read at species level that Kraken2 could only place at a
// genus. Where the two agree, two unrelated failure modes have been avoided.
//
// Additive, like every other classifier here. `--skip_kraken2` uses it alone.
//
// The resource bundle is the price of entry: a BWA index image, a sequence
// dictionary and a taxonomy database for the microbes, plus optionally a BWA
// image and a k-mer file for the host. Broad distributes prebuilt ones; this
// pipeline does not build them, because they are large, slow to build and
// specific to a reference choice the user has to make anyway.
//

include { GATK4_FASTQTOSAM } from '../../../modules/local/gatk4/fastqtosam/main'
include { GATK4_PATHSEQ    } from '../../../modules/local/gatk4/pathseq/main'
include { PATHSEQ_COMBINE  } from '../../../modules/local/pathseq/combine/main'

workflow TAXONOMY_PATHSEQ {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ]
    microbe_bwa_image // string: path to the microbe BWA index image
    microbe_dict // string: path to the microbe sequence dictionary
    taxonomy_db // string: path to the PathSeq taxonomy database
    host_bwa_image // string: host BWA index image, or null
    host_kmers // string: host k-mer file, or null
    save_bam // boolean: keep the per-read PathSeq BAM
    rank // string: PathSeq `type` to keep in the combined tables, or null

    main:

    //
    // PathSeq reads a SAM/BAM and nothing else, so the FASTQs are converted
    // first. No alignment happens in that step and no read is lost by it.
    //
    GATK4_FASTQTOSAM(ch_reads)

    GATK4_PATHSEQ(
        GATK4_FASTQTOSAM.out.bam,
        channel.value(file(microbe_bwa_image, checkIfExists: true)),
        channel.value(file(microbe_dict, checkIfExists: true)),
        channel.value(file(taxonomy_db, checkIfExists: true)),
        host_bwa_image ? channel.value(file(host_bwa_image, checkIfExists: true)) : channel.value([]),
        host_kmers ? channel.value(file(host_kmers, checkIfExists: true)) : channel.value([]),
        save_bam,
    )

    //
    // Sorted by sample id so the column order is reproducible across resumes -
    // `collect` makes no promise about the order files arrive in.
    //
    PATHSEQ_COMBINE(
        GATK4_PATHSEQ.out.scores
            .map { _meta, scores -> scores }
            .collect(sort: true)
            .map { tables -> [[id: 'pathseq'], tables] },
        rank,
    )

    emit:
    scores = GATK4_PATHSEQ.out.scores // channel: [ val(meta), path(txt) ]
    bam = GATK4_PATHSEQ.out.bam // channel: [ val(meta), path(bam) ]
    reads = PATHSEQ_COMBINE.out.reads // channel: [ val(meta), path(tsv) ]
    unambiguous = PATHSEQ_COMBINE.out.unambiguous // channel: [ val(meta), path(tsv) ]
    score = PATHSEQ_COMBINE.out.score // channel: [ val(meta), path(tsv) ]
    multiqc_files = PATHSEQ_COMBINE.out.mqc.map { _meta, mqc -> mqc } // channel: path(file)
}

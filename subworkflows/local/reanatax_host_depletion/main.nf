//
// Depletion of the host fraction, bulk and single-cell routes.
//
// The reads that survive every reference pass are the non-host fraction that
// everything downstream classifies. The single-cell route replaces bulk
// depletion rather than adding to it: STARsolo has to produce the cell-by-gene
// matrix and the unmapped reads in one pass, because the cell barcode is the
// only thing that joins them.
//

include { PREPARE_HOST_REFERENCE } from '../prepare_host_reference'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_FIRST } from '../prepare_host_reference'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_EXTRA1 } from '../prepare_host_reference'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_EXTRA2 } from '../prepare_host_reference'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FIRST } from '../host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_EXTRA1 } from '../host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_EXTRA2 } from '../host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FINAL } from '../host_depletion_hisat2'
include { SINGLECELL_STARSOLO } from '../singlecell_starsolo'
include { hostReferences } from '../utils_nfcore_reanatax_pipeline'
include { bowtie2References } from '../utils_nfcore_reanatax_pipeline'
include { hostMaskSettings } from '../utils_nfcore_reanatax_pipeline'
include { hostAligner } from '../host_depletion_hisat2'
include { minimap2Preset } from '../host_depletion_hisat2'

workflow HOST_DEPLETION {

    take:
    ch_trimmed_reads // channel: [ meta, [ reads ] ] after trimming (and optional rRNA removal)

    main:
    // MultiQC files mixed in by each depletion pass over the run.
    def ch_multiqc_files = channel.empty()

    def ch_nonhost_reads = ch_trimmed_reads
    def ch_qualimap = channel.empty()
    def ch_host_counts = channel.empty()
    def ch_hisat2_summaries = channel.empty()

    // References are applied in turn: the reads that survive one pass are the
    // input to the next, so a read has to fail against EVERY reference to be
    // called non-host. That is the point of pairing GRCh38 with T2T-CHM13 - the
    // second assembly holds the centromeric, satellite and structurally variant
    // sequence the first is missing, and those are exactly the regions whose
    // reads otherwise surface as spurious microbes. The same chain takes a
    // second organism the library really contains (a blood meal, a graft) or a
    // synthetic reference such as UniVec; validateInputParameters() caps the
    // count at maxHostPasses().
    // --host_bundle concatenates several genomes into one combined reference
    // that is depleted against ONCE, so it is handed in as a single 'fasta' host
    // reference - the same object the multi-pass --host logic already feeds the
    // primary index-build from. --host comma entries stay separate passes and
    // are not mixed with the bundle.
    def host_refs = params.host_bundle
        ? [[kind: 'fasta', value: params.host_bundle]]
        : hostReferences()
    def bowtie2_refs = bowtie2References()
    // Each library's aligner: HISAT2 for RNA, Bowtie2 for DNA under
    // --host_aligner auto. The reference subworkflows build only the indexes
    // some library actually needs.
    def ch_host_aligners = ch_trimmed_reads.map { meta, _reads -> hostAligner(meta, params.host_aligner) }
    // The minimap2 preset of each long-read library, so that one index is
    // built per preset in use.
    def ch_minimap2_presets = ch_trimmed_reads
        .filter { meta, _reads -> hostAligner(meta, params.host_aligner) == 'minimap2' }
        .map { meta, _reads -> minimap2Preset(meta, params.minimap2_preset) }
    // --host_mask, applied to every host pass alike.
    def host_mask = hostMaskSettings()
    def ch_solo = channel.empty()
    def ch_barcode_stats = channel.empty()

    if (params.single_cell) {
        //
        // The single-cell route replaces host depletion rather than adding to
        // it: STARsolo has to produce the cell-by-gene matrix and the unmapped
        // reads in one pass, because the cell barcode is the only thing that
        // joins them. Only the first --host reference is used - a second
        // depletion pass would have to re-align reads that no longer carry a
        // barcode read alongside them.
        //
        def primary = host_refs ? host_refs.first() : null
        if (!primary) {
            error("--single_cell needs a host reference: give --host (an accession, taxid, FASTA or index) or --star_index.")
        }
        if (host_refs.size() > 1) {
            log.warn("--single_cell uses only the first --host reference; the multi-pass depletion the bulk route does is not available here.")
        }

        PREPARE_HOST_REFERENCE(
            primary.kind == 'fasta' ? primary.value : null,
            null,
            primary.kind == 'accession' ? primary.value : null,
            primary.kind == 'taxid' ? primary.value : null,
            params.ncbi_group,
            params.gtf,
            'star',
            params.star_index ?: (primary.kind == 'index' ? primary.value : null),
            null,
            channel.empty(),
            channel.empty(),
            host_mask,
        )

        SINGLECELL_STARSOLO(
            ch_trimmed_reads,
            PREPARE_HOST_REFERENCE.out.index,
            PREPARE_HOST_REFERENCE.out.gtf,
            params.sc_whitelist,
        )
        ch_nonhost_reads = SINGLECELL_STARSOLO.out.reads
        ch_solo = SINGLECELL_STARSOLO.out.solo
        ch_barcode_stats = SINGLECELL_STARSOLO.out.barcode_stats
        ch_multiqc_files = ch_multiqc_files.mix(SINGLECELL_STARSOLO.out.multiqc_files)
    }
    else if (!params.skip_host_removal && host_refs) {
        def multi_pass = host_refs.size() > 1
        // The first reference is the primary one: --gtf describes it, so it is
        // the one indexed splice-aware, QC'd with Qualimap and counted. Every
        // later reference is a bare depletion pass - a GTF describes one
        // assembly, and handing the primary's annotation to a different one
        // would build a nonsense splice index, so only the primary gets it.
        // A Qualimap report of one reference's leftovers aligned to another
        // answers no question either, so the extra passes skip it.
        def primary = host_refs.first()

        if (multi_pass) {
            PREPARE_HOST_REFERENCE_FIRST(
                primary.kind == 'fasta' ? primary.value : null,
                primary.kind == 'index' ? primary.value : null,
                primary.kind == 'accession' ? primary.value : null,
                primary.kind == 'taxid' ? primary.value : null,
                params.ncbi_group,
                params.gtf,
                'hisat2',
                null,
                bowtie2_refs ? bowtie2_refs[0] : null,
                ch_host_aligners,
                ch_minimap2_presets,
                host_mask,
            )
            HOST_DEPLETION_FIRST(
                ch_trimmed_reads,
                PREPARE_HOST_REFERENCE_FIRST.out.index,
                PREPARE_HOST_REFERENCE_FIRST.out.fasta,
                params.save_host_bam,
                params.skip_qualimap,
                params.qualimap_gff,
                params.quantify_host ? params.gtf : null,
                params.require_both_mates_unmapped,
                params.hisat2_chunk_size,
                params.alignment_output_format,
                params.cleanup_intermediates,
                PREPARE_HOST_REFERENCE_FIRST.out.bowtie2_index,
                params.host_aligner,
                PREPARE_HOST_REFERENCE_FIRST.out.minimap2_index,
                params.minimap2_preset,
            )
            ch_nonhost_reads = HOST_DEPLETION_FIRST.out.reads
            ch_qualimap = HOST_DEPLETION_FIRST.out.qualimap
            ch_host_counts = HOST_DEPLETION_FIRST.out.host_counts
            ch_hisat2_summaries = HOST_DEPLETION_FIRST.out.summary
            ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_FIRST.out.multiqc_files)
        }

        // Everything between the primary and the last reference. Nextflow needs
        // a distinct alias per invocation, so these are unrolled rather than
        // looped; maxHostPasses() is the number of aliases declared above.
        if (host_refs.size() > 2) {
            def extra1 = host_refs[1]
            PREPARE_HOST_REFERENCE_EXTRA1(
                extra1.kind == 'fasta' ? extra1.value : null,
                extra1.kind == 'index' ? extra1.value : null,
                extra1.kind == 'accession' ? extra1.value : null,
                extra1.kind == 'taxid' ? extra1.value : null,
                params.ncbi_group,
                null,
                'hisat2',
                null,
                bowtie2_refs ? bowtie2_refs[1] : null,
                ch_host_aligners,
                ch_minimap2_presets,
                host_mask,
            )
            HOST_DEPLETION_EXTRA1(
                ch_nonhost_reads,
                PREPARE_HOST_REFERENCE_EXTRA1.out.index,
                PREPARE_HOST_REFERENCE_EXTRA1.out.fasta,
                params.save_host_bam,
                true,
                params.qualimap_gff,
                null,
                params.require_both_mates_unmapped,
                params.hisat2_chunk_size,
                params.alignment_output_format,
                params.cleanup_intermediates,
                PREPARE_HOST_REFERENCE_EXTRA1.out.bowtie2_index,
                params.host_aligner,
                PREPARE_HOST_REFERENCE_EXTRA1.out.minimap2_index,
                params.minimap2_preset,
            )
            ch_nonhost_reads = HOST_DEPLETION_EXTRA1.out.reads
            ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_EXTRA1.out.summary)
            ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_EXTRA1.out.multiqc_files)
        }

        if (host_refs.size() > 3) {
            def extra2 = host_refs[2]
            PREPARE_HOST_REFERENCE_EXTRA2(
                extra2.kind == 'fasta' ? extra2.value : null,
                extra2.kind == 'index' ? extra2.value : null,
                extra2.kind == 'accession' ? extra2.value : null,
                extra2.kind == 'taxid' ? extra2.value : null,
                params.ncbi_group,
                null,
                'hisat2',
                null,
                bowtie2_refs ? bowtie2_refs[2] : null,
                ch_host_aligners,
                ch_minimap2_presets,
                host_mask,
            )
            HOST_DEPLETION_EXTRA2(
                ch_nonhost_reads,
                PREPARE_HOST_REFERENCE_EXTRA2.out.index,
                PREPARE_HOST_REFERENCE_EXTRA2.out.fasta,
                params.save_host_bam,
                true,
                params.qualimap_gff,
                null,
                params.require_both_mates_unmapped,
                params.hisat2_chunk_size,
                params.alignment_output_format,
                params.cleanup_intermediates,
                PREPARE_HOST_REFERENCE_EXTRA2.out.bowtie2_index,
                params.host_aligner,
                PREPARE_HOST_REFERENCE_EXTRA2.out.minimap2_index,
                params.minimap2_preset,
            )
            ch_nonhost_reads = HOST_DEPLETION_EXTRA2.out.reads
            ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_EXTRA2.out.summary)
            ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_EXTRA2.out.multiqc_files)
        }

        // The final pass is always the one whose leftovers get classified, so it
        // is always this alias - which is what lets conf/modules.config publish
        // the non-host FASTQs from exactly one place regardless of how many
        // references were given. With a single reference nothing above ran and
        // ch_nonhost_reads is still the trimmed reads, so this pass is also the
        // primary one and takes the GTF and Qualimap.
        def last = host_refs.last()
        PREPARE_HOST_REFERENCE(
            last.kind == 'fasta' ? last.value : null,
            last.kind == 'index' ? last.value : null,
            last.kind == 'accession' ? last.value : null,
            last.kind == 'taxid' ? last.value : null,
            params.ncbi_group,
            multi_pass ? null : params.gtf,
            'hisat2',
            null,
            bowtie2_refs ? bowtie2_refs[host_refs.size() - 1] : null,
            ch_host_aligners,
            ch_minimap2_presets,
            host_mask,
        )
        HOST_DEPLETION_FINAL(
            ch_nonhost_reads,
            PREPARE_HOST_REFERENCE.out.index,
            PREPARE_HOST_REFERENCE.out.fasta,
            params.save_host_bam,
            multi_pass ? true : params.skip_qualimap,
            params.qualimap_gff,
            multi_pass ? null : (params.quantify_host ? params.gtf : null),
            params.require_both_mates_unmapped,
            params.hisat2_chunk_size,
            params.alignment_output_format,
            params.cleanup_intermediates,
            PREPARE_HOST_REFERENCE.out.bowtie2_index,
            params.host_aligner,
            PREPARE_HOST_REFERENCE.out.minimap2_index,
            params.minimap2_preset,
        )
        ch_nonhost_reads = HOST_DEPLETION_FINAL.out.reads
        ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_FINAL.out.summary)
        ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_FINAL.out.multiqc_files)
        if (!multi_pass) {
            ch_qualimap = HOST_DEPLETION_FINAL.out.qualimap
            ch_host_counts = HOST_DEPLETION_FINAL.out.host_counts
        }
    }
    emit:
    nonhost_reads    = ch_nonhost_reads     // channel: reads that failed every host reference
    qualimap         = ch_qualimap          // channel: Qualimap reports (bulk route)
    host_counts      = ch_host_counts       // channel: host-aligned read counts
    hisat2_summaries = ch_hisat2_summaries  // channel: per-pass alignment summaries
    solo             = ch_solo              // channel: STARsolo cell-by-gene matrix
    barcode_stats    = ch_barcode_stats     // channel: STARsolo barcode statistics
    multiqc_files    = ch_multiqc_files     // channel: files for MultiQC
}

//
// The classifiers that run beside the Kraken2 route rather than inside it:
// PRISM on the trimmed reads, GATK PathSeq and Metax on the non-host reads.
// None of them feeds the Kraken2 route's own filters; what comes back is their
// own tables, PRISM's optional drop pass over the Bracken table, and MultiQC
// sections. Kept out of the main workflow so that its body stays under the
// Groovy method-size limit.
//

include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_PRISM } from '../prepare_host_reference'
include { MINIMAP2_INDEX                     } from '../../../modules/local/minimap2/index/main'
include { PRISM_RUN                          } from '../../../modules/local/prism/run/main'
include { PRISM_COMBINE                      } from '../../../modules/local/prism/combine/main'
include { ABUNDANCE_FILTER as FILTER_PRISM   } from '../../../modules/local/abundance/filter/main'
include { TAXONOMY_PATHSEQ                   } from '../taxonomy_pathseq'
include { TAXONOMY_METAX                     } from '../taxonomy_metax'
include { externalToolProblem                } from '../utils_nfcore_reanatax_pipeline'
include { hostReferences                     } from '../utils_nfcore_reanatax_pipeline'

workflow ADDITIONAL_CLASSIFIERS {

    take:
    ch_trimmed_reads // channel: [ val(meta), [ path(fastq) ] ] - after trimming, before host removal
    ch_nonhost_reads // channel: [ val(meta), [ path(fastq) ] ] - after host removal
    ch_taxonomy_drop // channel: path(txt) - the Kraken2 route's taxonomy verdicts, or empty
    ch_bracken_filtered // channel: [ val(meta), path(tsv) ] - the filtered combined Bracken table, or empty

    main:

    def ch_multiqc_files = channel.empty()

    //
    // MODULE: PRISM, on the TRIMMED reads.
    //
    // The one tool here that must not be handed the non-host fraction. PRISM
    // does its own host removal with STAR and minimap2, and several of its
    // forty features are statistics of that removal - how much of the library
    // was host, how much was unclassified, how the k-mers of a candidate taxon
    // distribute across ranks. Give it reads that have already had the host
    // taken out and those columns empty, which changes what its model is
    // scoring without changing anything that would look wrong in the output.
    //
    // It is also not a classifier: it takes candidate taxa and tries to
    // CONFIRM them with independent evidence - full-length BLAST against nt,
    // the GenBank annotation the alignments land on, host mapping and k-mer
    // composition - and returns a probability per taxon.
    //
    // Not installed by the pipeline: when it cannot be reached, startup has
    // already warned and the run goes on without it.
    //
    if (params.run_prism && !externalToolProblem('prism')) {
        //
        // PRISM does its own host removal with STAR and minimap2, on the
        // TRIMMED reads. --prism_star_genome_dir / --prism_minimap2_index exist
        // for callers that already have the indices; when they are absent the
        // indices are built here from the run's host genome - the --host_bundle
        // if one was given, else the --host FASTA/accession/taxid - so a run
        // needs no hand-built PRISM indices. The built STAR/minimap2 indexes
        // live in the run's work directory and stay valid for its duration.
        //
        def ch_prism_star_genome
        def ch_prism_minimap2_index
        if (params.prism_star_genome_dir) {
            ch_prism_star_genome = channel.value(params.prism_star_genome_dir)
            ch_prism_minimap2_index = channel.value(params.prism_minimap2_index)
        }
        else {
            def prism_host = params.host_bundle
                ? [kind: 'fasta', value: params.host_bundle]
                : hostReferences().find { ref -> ref.kind in ['fasta', 'accession', 'taxid'] }
            PREPARE_HOST_REFERENCE_PRISM(
                prism_host.kind == 'fasta' ? prism_host.value : null,
                null,
                prism_host.kind == 'accession' ? prism_host.value : null,
                prism_host.kind == 'taxid' ? prism_host.value : null,
                params.ncbi_group,
                null,
                'star',
                null,
                null,
                channel.empty(),
                channel.empty(),
                null,
            )
            ch_prism_star_genome = PREPARE_HOST_REFERENCE_PRISM.out.index
                .map { _meta, index -> index.toString() }
                .first()
            // PREPARE with aligner 'star' builds no minimap2 index, so build one
            // for PRISM's short-read host pass (-x sr) from the same host FASTA.
            MINIMAP2_INDEX(
                PREPARE_HOST_REFERENCE_PRISM.out.fasta
                    .map { meta, genome -> [meta + [preset: 'sr'], genome] }
            )
            ch_prism_minimap2_index = MINIMAP2_INDEX.out.index
                .map { _meta, index -> index.toString() }
                .first()
        }
        PRISM_RUN(
            ch_trimmed_reads,
            channel.value(file(params.prism_path, checkIfExists: true)),
            params.prism_kraken_db ?: params.kraken2_db,
            params.prism_blast_db,
            ch_prism_star_genome,
            ch_prism_minimap2_index,
            params.prism_min_qcovs,
            params.prism_min_read_per,
            params.prism_min_uniq_frac,
            params.prism_max_sample,
            params.prism_barcode_only,
            params.prism_bin_dir ?: '',
            params.prism_decoy_taxids ?: (params.host_kmer_taxid ?: params.host_carryover_taxid).toString(),
        )
        PRISM_COMBINE(
            PRISM_RUN.out.counts
                .map { _meta, table -> table }
                .collect(sort: true)
                .map { tables -> [[id: 'reanatax'], tables] },
            params.prism_score_threshold,
            params.prism_min_reads,
        )
        ch_multiqc_files = ch_multiqc_files.mix(PRISM_COMBINE.out.mqc.map { _meta, mqc -> mqc })

        //
        // A second filter pass, for the same reason the gene-diversity one
        // needs one: PRISM runs beside the Kraken2 route rather than inside
        // it, so its verdict cannot reach the ABUNDANCE_FILTER in the taxonomy
        // subworkflow. The thresholds are zeroed so this pass removes PRISM's
        // list and nothing else, and the table is named for what was applied.
        //
        if (params.prism_filter && !params.skip_kraken2 && !params.skip_bracken) {
            FILTER_PRISM(
                ch_bracken_filtered,
                PRISM_COMBINE.out.drop_list.map { _meta, list -> list }.collect(sort: true),
                [],
                0,
                1,
                0,
            )
        }
    }

    //
    // SUBWORKFLOW: GATK PathSeq, on the same non-host reads.
    //
    // The only classifier here that aligns rather than matching k-mers, so it
    // fails differently from the other two - and agreement between an
    // alignment-based and a k-mer-based route is a much stronger statement
    // about a taxon than either one alone. Additive; `--skip_kraken2` uses it
    // on its own.
    //
    if (params.pathseq_microbe_bwa_image) {
        TAXONOMY_PATHSEQ(
            ch_nonhost_reads,
            params.pathseq_microbe_bwa_image,
            params.pathseq_microbe_dict,
            params.pathseq_taxonomy_db,
            params.pathseq_host_bwa_image,
            params.pathseq_host_kmers,
            params.pathseq_save_bam,
            params.pathseq_rank,
            ch_taxonomy_drop,
        )
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_PATHSEQ.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: Metax, on the same non-host reads.
    //
    // Additive like MetaPhlAn and sylph, with `--run_metax --skip_kraken2` as
    // the replacement. Metax does align every read and writes a per-read
    // classification, but in its own format and not on the path the cell,
    // host-k-mer and minimizer logic read, so a Metax-only run carries the same
    // `--skip_kraken2` restrictions. The pre-built databases include the human
    // reference where applicable (the 2022 one names GRCh38), and with it in
    // the index, host reads under `--skip_host_removal` land on Homo sapiens
    // rather than on whatever microbe they happen to resemble. Like PRISM it is
    // not installed by the pipeline, and skipped with a startup warning when it
    // cannot be reached.
    //
    if (params.run_metax && !externalToolProblem('metax')) {
        TAXONOMY_METAX(
            ch_nonhost_reads,
            params.metax_db,
            params.metax_dmp_dir,
            params.metax_save_readclassifications,
            ch_taxonomy_drop,
            params.metax_bin_dir ?: '',
        )
    }

    emit:
    multiqc_files = ch_multiqc_files // channel: path(file)
}

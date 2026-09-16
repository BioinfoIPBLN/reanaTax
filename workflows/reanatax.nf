/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_DOWNLOAD_FASTQDL    } from '../subworkflows/local/fastq_download_fastqdl'
include { CLEANUP_INTERMEDIATES     } from '../modules/local/cleanup/main'
include { PREPARE_HOST_REFERENCE    } from '../subworkflows/local/prepare_host_reference'
include { FASTQ_QC_TRIM             } from '../subworkflows/local/fastq_qc_trim'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_FIRST } from '../subworkflows/local/prepare_host_reference'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_EXTRA1 } from '../subworkflows/local/prepare_host_reference'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_EXTRA2 } from '../subworkflows/local/prepare_host_reference'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FIRST } from '../subworkflows/local/host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_EXTRA1 } from '../subworkflows/local/host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_EXTRA2 } from '../subworkflows/local/host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FINAL } from '../subworkflows/local/host_depletion_hisat2'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_UNIVEC } from '../subworkflows/local/prepare_host_reference'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_UNIVEC } from '../subworkflows/local/host_depletion_hisat2'
include { TAXONOMY_KRAKEN2_BRACKEN  } from '../subworkflows/local/taxonomy_kraken2_bracken'
include { TAXONOMY_KRAKENUNIQ      } from '../subworkflows/local/taxonomy_krakenuniq'
include { TAXONOMY_METAPHLAN       } from '../subworkflows/local/taxonomy_metaphlan'
include { TAXONOMY_PATHSEQ         } from '../subworkflows/local/taxonomy_pathseq'
include { DIFFERENTIAL_ABUNDANCE    } from '../subworkflows/local/differential_abundance'
include { HOST_EXPRESSION           } from '../subworkflows/local/host_expression'
include { TARGETED_TAXON            } from '../subworkflows/local/targeted_taxon'
include { QUANTIFY_KALLISTO         } from '../subworkflows/local/quantify_kallisto'
include { FUNCTIONAL_HUMANN         } from '../subworkflows/local/functional_humann'
include { AI_ANNOTATE_REPORTS       } from '../subworkflows/local/ai_annotate_reports'
include { RRNA_REMOVAL_SORTMERNA    } from '../subworkflows/local/rrna_removal_sortmerna'
include { SINGLECELL_STARSOLO       } from '../subworkflows/local/singlecell_starsolo'
include { LLM_INSIGHT               } from '../modules/local/llm/insight/main'
include { READ_ACCOUNTING           } from '../modules/local/read/accounting/main'
include { EXPLOREMETATAX_BUNDLE     } from '../modules/local/exploremetatax/bundle/main'
include { POLYA_CARRYOVER           } from '../modules/local/polya/carryover/main'
include { POLYA_MERGE               } from '../modules/local/polya/merge/main'
include { DIVERSITY                 } from '../modules/local/diversity/main'
include { PRISM_RUN                 } from '../modules/local/prism/run/main'
include { PRISM_COMBINE             } from '../modules/local/prism/combine/main'
include { SPARCC_PREPARE            } from '../modules/local/sparcc/prepare/main'
include { SPARCC_CORRELATION        } from '../modules/local/sparcc/correlation/main'
include { SPARCC_NETWORK            } from '../modules/local/sparcc/network/main'
include { ABUNDANCE_GENEDIV as GENEDIV_GENE } from '../modules/local/abundance/genediv/main'
include { ABUNDANCE_GENEDIV as GENEDIV_PRODUCT } from '../modules/local/abundance/genediv/main'
include { ABUNDANCE_FILTER as FILTER_GENEDIV } from '../modules/local/abundance/filter/main'
include { ABUNDANCE_FILTER as FILTER_PRISM } from '../modules/local/abundance/filter/main'
include { HOSTCOUNTS_MERGE          } from '../modules/local/hostcounts/merge/main'
include { SCTAXA_COUNTS             } from '../modules/local/sctaxa/counts/main'
include { SCTAXA_DENOISE            } from '../modules/local/sctaxa/denoise/main'
include { SCTAXA_ENRICHMENT         } from '../modules/local/sctaxa/enrichment/main'
include { SCTAXA_FILTER             } from '../modules/local/sctaxa/filter/main'
include { SCTAXA_HOSTDE             } from '../modules/local/sctaxa/hostde/main'
include { SCTAXA_AMBIENT            } from '../modules/local/sctaxa/ambient/main'
include { SCTAXA_PLATE              } from '../modules/local/sctaxa/plate/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { aiInsightOptions          } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { hostReferences            } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { minimizerThresholds      } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { hostKmerSettings         } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { decontamSettings         } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { shuffleSettings         } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { controlSettings          } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { assemblySettings         } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { hostCladeSettings        } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { scChemistry              } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { meanReadLength            } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { brackenDistributions      } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { resolveBrackenReadLength  } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { daMethods                } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { hostDeMethods            } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { targetReference          } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow REANATAX {

    take:
    ch_samplesheet // channel: [ val(meta), [ path(fastq) ] ] from --input / --input_dir
    ch_accessions // channel: [ val(meta), val(accession) ] from --input_accessions
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // Which AI annotations to run. Empty unless --llm_endpoint is set, so the
    // pipeline never contacts an external service by default.
    def ai_options = aiInsightOptions(params.ai_insights, params.llm_endpoint)

    //
    // SUBWORKFLOW: Fetch raw reads from ENA/SRA when accessions were given
    //
    FASTQ_DOWNLOAD_FASTQDL(
        ch_accessions,
        params.group_runs_by,
    )

    def ch_raw_reads = ch_samplesheet.mix(FASTQ_DOWNLOAD_FASTQDL.out.reads)

    //
    // SUBWORKFLOW: Raw read QC, adapter/quality trimming, post-trim QC
    //
    FASTQ_QC_TRIM(
        ch_raw_reads,
        params.adapter_fasta ? file(params.adapter_fasta, checkIfExists: true) : [],
        params.skip_fastqc,
        params.skip_trimming,
        params.save_trimmed_fail,
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM.out.multiqc_files)

    //
    // --cleanup_intermediates: drop the downloaded FASTQ once trimming and raw
    // FastQC have both read it.
    //
    // Only the DOWNLOADED reads, never the samplesheet's. A downloaded FASTQ is
    // a copy of something still in the archive and the pipeline can fetch it
    // again; a samplesheet's FASTQ is the user's own file, and the module
    // refuses to touch anything outside the work directory in any case. The
    // join on `raw_finished` is what makes the deletion wait, and that channel
    // is empty when trimming is skipped - the raw reads are then the working
    // read set, and nothing here is finished with them.
    //
    if (params.cleanup_intermediates) {
        CLEANUP_INTERMEDIATES(
            FASTQ_DOWNLOAD_FASTQDL.out.reads.join(FASTQ_QC_TRIM.out.raw_finished)
        )
    }

    //
    // SUBWORKFLOW: Strip rRNA before host depletion.
    //
    // Everything downstream reads `ch_trimmed_reads` rather than the trimming
    // subworkflow's output directly, so inserting this step cannot leave a
    // later stage silently reading the unfiltered reads.
    //
    def ch_trimmed_reads = FASTQ_QC_TRIM.out.reads
    def ch_sortmerna_logs = channel.empty()

    if (params.remove_rrna) {
        RRNA_REMOVAL_SORTMERNA(
            FASTQ_QC_TRIM.out.reads,
            params.sortmerna_db,
            params.sortmerna_index,
        )
        ch_trimmed_reads = RRNA_REMOVAL_SORTMERNA.out.reads
        ch_sortmerna_logs = RRNA_REMOVAL_SORTMERNA.out.log
        ch_multiqc_files = ch_multiqc_files.mix(RRNA_REMOVAL_SORTMERNA.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: Split host from non-host reads
    //
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
    def host_refs = hostReferences()
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
        )
        ch_nonhost_reads = HOST_DEPLETION_FINAL.out.reads
        ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_FINAL.out.summary)
        ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_FINAL.out.multiqc_files)
        if (!multi_pass) {
            ch_qualimap = HOST_DEPLETION_FINAL.out.qualimap
            ch_host_counts = HOST_DEPLETION_FINAL.out.host_counts
        }
    }

    //
    // SUBWORKFLOW: Strip synthetic sequence (UniVec).
    //
    // Its own pass rather than another --host entry, for one reason: a vector
    // database needs DIFFERENT alignment stringency from a host genome.
    // --very-sensitive expands to `--score-min L,0,-1`, which lets a 100 bp
    // read carry roughly twenty-five mismatches and still align. Against a
    // genome that is what catches diverged and repetitive host sequence;
    // against a 3 kb plasmid backbone it deletes real reads. So this pass has
    // its own arguments, end-to-end and at HISAT2's own default score
    // threshold, and does not inherit --hisat2_args.
    //
    // Runs LAST, on whatever produced the non-host reads - so it also works
    // with --skip_host_removal (screen vectors with no host genome at all) and
    // after the single-cell route.
    //
    if (params.univec) {
        PREPARE_HOST_REFERENCE_UNIVEC(
            params.univec_source,
            null,
            null,
            null,
            params.ncbi_group,
            null,
            'hisat2',
            null,
        )
        HOST_DEPLETION_UNIVEC(
            ch_nonhost_reads,
            PREPARE_HOST_REFERENCE_UNIVEC.out.index,
            PREPARE_HOST_REFERENCE_UNIVEC.out.fasta,
            params.univec_save_bam,
            true,
            null,
            null,
            params.require_both_mates_unmapped,
            params.hisat2_chunk_size,
            params.alignment_output_format,
            params.cleanup_intermediates,
        )
        ch_nonhost_reads = HOST_DEPLETION_UNIVEC.out.reads
        ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_UNIVEC.out.summary)
        ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_UNIVEC.out.multiqc_files)
    }

    //
    // MODULE: Is the microbial signal genuine poly(A) capture or carry-over?
    //
    // Only meaningful for poly(A)-selected libraries, which is why it is opt-in:
    // on a metagenome or an rRNA-depleted library the question does not arise.
    //
    def ch_polya_mqc = channel.empty()

    if (params.run_polya_check) {
        POLYA_CARRYOVER(ch_nonhost_reads)
        POLYA_MERGE(
            POLYA_CARRYOVER.out.tsv
                .map { _meta, tsv -> tsv }
                .collect(sort: true)
                .map { tables -> [[id: 'reanatax'], tables] }
        )
        ch_polya_mqc = POLYA_MERGE.out.mqc
        ch_multiqc_files = ch_multiqc_files.mix(POLYA_MERGE.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // SUBWORKFLOW: Taxonomic classification of the non-host fraction
    //
    def ch_kraken2_report = channel.empty()
    def ch_bracken = channel.empty()
    def ch_bracken_combined = channel.empty()

    if (!params.skip_kraken2) {
        //
        // Bracken's -r has to match the k-mer distribution the database was
        // built with, not merely the reads. Getting it wrong does not fail, it
        // silently returns the wrong abundances - so the length is measured from
        // fastp rather than assumed, then snapped to a distribution the database
        // actually ships. It rides on the meta so it stays per-sample: a cohort
        // pulled from several studies will not share one read length.
        //
        def bracken_dists = brackenDistributions(params.bracken_db ?: params.kraken2_db)
        def ch_reads_for_tax = ch_nonhost_reads
            .join(FASTQ_QC_TRIM.out.fastp_json, remainder: true)
            // remainder: true exists so a sample with reads but NO fastp JSON
            // (--skip_trimming) still gets classified. It also lets the mirror
            // case through - a sample with a JSON and no reads - as a null path,
            // which surfaces as "Path value cannot be null" on a process whose
            // own error message then fails to render. Keeping the first case and
            // dropping the second is what this filter does; a sample that lost
            // its reads upstream is an upstream bug, and it should not be
            // rediscovered here as an unrecoverable crash.
            .filter { _meta, reads, _json -> reads }
            .map { meta, reads, json ->
                def observed = json ? meanReadLength(json) : null
                [meta + [bracken_r: resolveBrackenReadLength(observed, bracken_dists, params.bracken_read_length)], reads]
            }

        TAXONOMY_KRAKEN2_BRACKEN(
            ch_reads_for_tax,
            params.kraken2_db,
            params.bracken_db,
            params.kraken2_save_reads,
            params.kraken2_save_readclassifications,
            params.skip_bracken,
            params.skip_krona,
            params.minimizer_filter,
            minimizerThresholds(),
            params.host_kmer_filter,
            hostKmerSettings(),
            decontamSettings(),
            shuffleSettings(),
            controlSettings(),
            assemblySettings(),
            hostCladeSettings(),
            params.drop_host_taxon,
            params.host_carryover_taxid,
            params.kraken2_use_daemon,
            params.min_rel_abundance,
            params.min_samples,
            params.min_reads,
            params.export_biom,
            params.biom_metadata ?: params.da_metadata,
        )
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_KRAKEN2_BRACKEN.out.multiqc_files)
        ch_kraken2_report = TAXONOMY_KRAKEN2_BRACKEN.out.report
        ch_bracken = TAXONOMY_KRAKEN2_BRACKEN.out.bracken
        ch_bracken_combined = TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined
    }

    //
    // MODULE: The cell-by-taxon matrix.
    //
    // Built from the SAME classification the pseudobulk profile is: the reads
    // went through Kraken2 once, and the barcode was carried through in the
    // read name. So the two views cannot disagree, and every filter that acts
    // on the pseudobulk tables acts on the reads behind the matrix too.
    //
    def ch_cell_taxa = channel.empty()

    if (params.single_cell && !params.skip_kraken2) {
        // Joined on sample id rather than zipped: the three channels are
        // produced by different processes and their emission order is not
        // guaranteed to match.
        def ch_sc_input = TAXONOMY_KRAKEN2_BRACKEN.out.classifiedreads
            .map { meta, reads -> [meta.id, meta, reads] }
            .join(TAXONOMY_KRAKEN2_BRACKEN.out.report.map { meta, report -> [meta.id, report] })
            .join(ch_nonhost_reads.map { meta, reads -> [meta.id, reads] })
            .map { _id, meta, reads, report, fastq -> [meta, reads, report, fastq] }

        SCTAXA_COUNTS(
            ch_sc_input,
            params.sc_host_taxid ?: params.host_kmer_taxid ?: params.host_carryover_taxid,
            params.sc_min_frac,
            params.sc_max_homopolymer,
            params.sc_ranks,
            params.sc_umi_dedup,
            params.sc_min_umis,
        )
        ch_cell_taxa = SCTAXA_COUNTS.out.counts

        if (params.sc_kmer_denoise) {
            SCTAXA_DENOISE(
                ch_sc_input.map { meta, reads, _report, fastq -> [meta, reads, fastq] },
                params.sc_host_taxid ?: params.host_kmer_taxid ?: params.host_carryover_taxid,
                params.sc_kmer_len,
                params.sc_min_barcodes,
                params.sc_correlation_p,
                params.sc_adjust,
                params.sc_max_barcodes_per_taxon,
            )
            ch_multiqc_files = ch_multiqc_files.mix(SCTAXA_DENOISE.out.mqc.map { _meta, mqc -> mqc })
        }

        //
        // MODULE: The host half of SAHMI.
        //
        // The cell-by-taxon matrix is the instrument, not the result. What the
        // paper actually reports is what a cell's own transcriptome does when
        // it is carrying something, and both halves of that comparison come
        // out of the same STARsolo pass keyed on the same corrected barcode.
        //
        //
        // MODULE: Cell association against the ambient pool.
        //
        // An annotation, deliberately, not a filter. Empty droplets hold
        // reagent contaminants AND genuine extracellular organisms from the
        // tissue, and nothing at the droplet level separates them - so this
        // reports what it can measure (is the taxon where the cells are?) and
        // leaves the contamination question to --decontam, which is the only
        // step here given an external measurement of the kit.
        //
        if (params.sc_ambient) {
            SCTAXA_AMBIENT(
                SINGLECELL_STARSOLO.out.solo
                    .map { meta, solo -> [meta.id, meta, solo] }
                    .join(SCTAXA_COUNTS.out.counts.map { meta, table -> [meta.id, table] })
                    .map { _id, meta, solo, table -> [meta, solo, table] },
                params.sc_solo_features.tokenize(',')[0].trim(),
                params.sc_min_umis,
                params.sc_ambient_min_empty_umis,
                params.sc_ambient_max_empty_umis,
                params.sc_ambient_min_droplets,
                params.sc_ambient_enriched_ratio,
                params.sc_ambient_depleted_ratio,
                params.sc_enrichment_p,
                params.sc_ambient_drop,
            )
            ch_multiqc_files = ch_multiqc_files.mix(SCTAXA_AMBIENT.out.mqc.map { _meta, mqc -> mqc })
        }

        if (params.sc_host_de) {
            SCTAXA_HOSTDE(
                SINGLECELL_STARSOLO.out.solo
                    .map { meta, solo -> [meta.id, meta, solo] }
                    .join(SCTAXA_COUNTS.out.counts.map { meta, table -> [meta.id, table] })
                    .map { _id, meta, solo, table -> [meta, solo, table] },
                params.sc_cell_metadata ? file(params.sc_cell_metadata, checkIfExists: true) : [],
                params.sc_solo_features.tokenize(',')[0].trim(),
                params.sc_min_umis,
                params.sc_host_de_min_cells,
                params.sc_host_de_min_pct,
                params.sc_host_de_logfc,
                params.sc_host_de_top_taxa,
                params.sc_enrichment_p,
                params.sc_host_de_force_pooled,
            )
            ch_multiqc_files = ch_multiqc_files.mix(SCTAXA_HOSTDE.out.mqc.map { _meta, mqc -> mqc })
        }
    }

    //
    // MODULE: The same matrix, for plate-based single-cell data.
    //
    // In Smart-seq2 every well is its own library, so there is no barcode and
    // no STARsolo: each cell went through the ordinary bulk route and the
    // matrix is those per-cell reports stacked. CSI-Microbes keeps this as a
    // separate route for the same reason, and everything downstream of the
    // matrix is shared with the droplet one.
    //
    if (params.sc_plate_based && !params.skip_kraken2) {
        SCTAXA_PLATE(
            TAXONOMY_KRAKEN2_BRACKEN.out.report
                .map { _meta, report -> report }
                .collect(sort: true)
                .map { reports -> [[id: 'reanatax'], reports] },
            file(params.sc_cell_metadata, checkIfExists: true),
            params.sc_ranks,
            params.sc_plate_min_reads,
        )
        ch_cell_taxa = SCTAXA_PLATE.out.counts
    }

    //
    // MODULE: Which host cell types carry which microbes.
    //
    // Cell types are the one thing this pipeline cannot supply itself -
    // STARsolo produces a count matrix, not clusters, and a plate-based run
    // produces one library per well - so this waits for annotations rather
    // than guessing at them. Fed by whichever route built the matrix.
    //
    //
    // MODULE: the evidence filters, applied to the matrix as well as the tables.
    //
    // Without this the drop lists reach only the cohort tables and the
    // enrichment test still tests every background taxon - paying the
    // multiple-testing cost for taxa the pipeline has already judged false.
    //
    if (params.sc_apply_drop_list && (params.single_cell || params.sc_plate_based) && !params.skip_kraken2) {
        SCTAXA_FILTER(
            ch_cell_taxa,
            TAXONOMY_KRAKEN2_BRACKEN.out.drop_list,
            // --negative_controls reaches the matrix here or nowhere: the
            // cell-by-taxon table is built from the per-read assignments, not
            // from the combined tables the filter rewrote, so without this the
            // cohort profile and the single-cell profile would disagree about
            // every cell it zeroed. Empty on a run without controls.
            TAXONOMY_KRAKEN2_BRACKEN.out.control_cells,
        )
        ch_cell_taxa = SCTAXA_FILTER.out.counts
        ch_multiqc_files = ch_multiqc_files.mix(SCTAXA_FILTER.out.mqc.map { _meta, mqc -> mqc })
    }

    if (params.sc_cell_metadata && (params.single_cell || params.sc_plate_based) && !params.skip_kraken2) {
        SCTAXA_ENRICHMENT(
            ch_cell_taxa
                .map { _meta, table -> table }
                .collect(sort: true)
                .map { tables -> [[id: 'reanatax'], tables] },
            file(params.sc_cell_metadata, checkIfExists: true),
            params.sc_cooccurrence_min_cells,
            params.sc_enrichment_p,
            // The negative controls, kept out of the test entirely. They are
            // not a cell type, and after --negative_controls they are the one
            // group that was not filtered - a control cannot be judged against
            // itself - so leaving them in makes every reagent contaminant look
            // hugely enriched in the blanks. Measured on the CSI-Microbes
            // plate: 41 of 43 significant results were exactly that.
            controlSettings()?.controls ?: '',
        )
        ch_multiqc_files = ch_multiqc_files.mix(SCTAXA_ENRICHMENT.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // SUBWORKFLOW: KrakenUniq, on the same non-host reads.
    //
    // Its own database, so its own route rather than a switch on the Kraken2
    // one. Additive by default: with Kraken2 left on you get both profiles over
    // one set of reads, which is the only way to see whether a taxon Kraken2
    // called is backed by distinct k-mers or by reads on a single conserved
    // stretch. `--skip_kraken2` uses it on its own.
    //
    def ch_krakenuniq_report = channel.empty()

    if (params.krakenuniq_db) {
        TAXONOMY_KRAKENUNIQ(
            ch_nonhost_reads,
            params.krakenuniq_db,
            params.krakenuniq_save_reads,
            params.krakenuniq_save_readclassifications,
            params.skip_krona,
            params.minimizer_filter,
            minimizerThresholds(),
            params.min_rel_abundance,
            params.min_samples,
            params.min_reads,
        )
        ch_krakenuniq_report = TAXONOMY_KRAKENUNIQ.out.report
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_KRAKENUNIQ.out.multiqc_files)
    }

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
    if (params.run_prism) {
        PRISM_RUN(
            ch_trimmed_reads,
            channel.value(file(params.prism_path, checkIfExists: true)),
            params.prism_kraken_db ?: params.kraken2_db,
            params.prism_blast_db,
            params.prism_star_genome_dir,
            params.prism_minimap2_index,
            params.prism_min_qcovs,
            params.prism_min_read_per,
            params.prism_min_uniq_frac,
            params.prism_max_sample,
            params.prism_barcode_only,
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
                TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined_filtered,
                PRISM_COMBINE.out.drop_list.map { _meta, list -> list }.collect(sort: true),
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
        )
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_PATHSEQ.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: MetaPhlAn, on the same non-host reads.
    //
    // Additive rather than exclusive: `--run_metaphlan --skip_kraken2` uses it
    // instead of Kraken2, and leaving Kraken2 on runs both over one set of
    // reads, which is the only way to see where a marker-gene profile and a
    // whole-read classifier disagree on the same library.
    //
    def ch_metaphlan_merged = channel.empty()
    def ch_metaphlan_profile = channel.empty()
    def ch_humann_tables = channel.empty()
    def ch_genediv_drop = channel.empty()
    if (params.run_metaphlan) {
        TAXONOMY_METAPHLAN(
            ch_nonhost_reads,
            params.metaphlan_db,
            params.metaphlan_index,
        )
        ch_metaphlan_merged = TAXONOMY_METAPHLAN.out.merged
        ch_metaphlan_profile = TAXONOMY_METAPHLAN.out.profile

    }
    //
    // SUBWORKFLOW: What the community is doing, not just who is in it.
    //
    if (params.run_humann) {
        // Bracken's kreport, not Kraken2's, whenever Bracken ran. HUMAnN
        // selects pangenomes from the `s__` lines of the profile and ignores
        // everything above species, and re-estimating abundances down to
        // species is exactly what Bracken is for: a read Kraken2 could only
        // place at a genus is invisible to HUMAnN until Bracken pushes it
        // down. Same ordering the FGCZ HUMAnN SUSHI app uses.
        //
        // The two differ in coverage as well as in numbers: Bracken only runs
        // on reports with something classified, so a library where Kraken2
        // classified nothing drops out here rather than reaching HUMAnN with
        // an empty profile.
        def ch_humann_profile = params.skip_bracken
            ? TAXONOMY_KRAKEN2_BRACKEN.out.report
            : TAXONOMY_KRAKEN2_BRACKEN.out.bracken_report

        FUNCTIONAL_HUMANN(
            ch_nonhost_reads,
            ch_humann_profile,
            params.humann_nucleotide_db,
            params.humann_protein_db,
            params.humann_utility_db,
            params.humann_regroup,
            params.humann_renorm,
        )
        // Only the tables exploreMetaTax can name: the renormalised one is
        // left out because what it holds depends on --humann_regroup.
        ch_humann_tables = FUNCTIONAL_HUMANN.out.pathabundance
            .mix(FUNCTIONAL_HUMANN.out.genefamilies)
            .mix(FUNCTIONAL_HUMANN.out.regrouped)

        //
        // MODULE: Breadth of a taxon's evidence across its own genome.
        //
        // PRISM's argument, answered from a table HUMAnN has already written.
        // Every filter in the taxonomy subworkflow counts reads or the k-mers
        // behind them; none can see WHERE those reads landed, and a thousand
        // reads on one conserved gene passes all of them. HUMAnN's
        // species-stratified gene families say exactly that, at no extra cost.
        //
        // Twice, when a regrouping was asked for: gene families and the
        // products they encode are different questions, and a taxon can reach
        // dozens of UniRef90 families that all encode the same thing.
        //
        if (params.gene_diversity_filter && !params.skip_bracken) {
            def ch_genediv_bracken = TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined
                .map { _meta, table -> table }

            GENEDIV_GENE(
                FUNCTIONAL_HUMANN.out.genefamilies
                    .map { _meta, table -> table }
                    .collect(sort: true)
                    .map { tables -> [[id: 'reanatax'], tables] },
                ch_genediv_bracken,
                'gene',
                params.gene_diversity_min_families,
                params.gene_diversity_max_top_fraction,
                params.gene_diversity_min_abundance,
            )
            ch_genediv_drop = GENEDIV_GENE.out.drop_list.map { _meta, list -> list }
            ch_multiqc_files = ch_multiqc_files.mix(GENEDIV_GENE.out.mqc.map { _meta, mqc -> mqc })

            if (params.humann_regroup) {
                GENEDIV_PRODUCT(
                    FUNCTIONAL_HUMANN.out.regrouped
                        .map { _meta, table -> table }
                        .collect(sort: true)
                        .map { tables -> [[id: 'reanatax'], tables] },
                    ch_genediv_bracken,
                    'product',
                    params.gene_diversity_min_families,
                    params.gene_diversity_max_top_fraction,
                    params.gene_diversity_min_abundance,
                )
                ch_genediv_drop = ch_genediv_drop.mix(GENEDIV_PRODUCT.out.drop_list.map { _meta, list -> list })
                ch_multiqc_files = ch_multiqc_files.mix(GENEDIV_PRODUCT.out.mqc.map { _meta, mqc -> mqc })
            }

            //
            // A second filter pass rather than a contribution to the first.
            // HUMAnN runs downstream of the taxonomy subworkflow, so its
            // verdict cannot travel back to the ABUNDANCE_FILTER inside it; the
            // thresholds are zeroed here so this pass removes the drop list and
            // nothing else, and the table it produces is named for what was
            // applied to it.
            //
            FILTER_GENEDIV(
                TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined_filtered,
                ch_genediv_drop.collect(sort: true),
                0,
                1,
                0,
            )
        }
    }

    //
    // MODULE: Which taxa travel together.
    //
    // A co-occurrence network over the cohort, inferred by SparCC. Taxonomic
    // profiles are compositional - the counts sum to a library size with no
    // biological meaning - so an ordinary correlation matrix manufactures
    // negative associations everywhere and a spurious positive block among the
    // rare taxa. SparCC infers the correlations of the underlying basis from
    // log-ratio variances instead, which is the same reason ANCOM-BC2 and
    // ALDEx2 are this pipeline's differential-abundance methods.
    //
    // Off by default and refused outright on a small cohort: SPARCC_PREPARE
    // will not start a run whose bootstrap p-values could not reach the BH
    // threshold for the number of pairs being tested, because an empty network
    // from such a run means "this cohort cannot be asked" and not "there are no
    // associations".
    //
    if (params.run_sparcc && !params.skip_kraken2 && !params.skip_bracken) {
        SPARCC_PREPARE(
            TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined_filtered,
            params.sparcc_min_prevalence,
            params.sparcc_min_reads,
            params.sparcc_top_taxa,
            params.sparcc_permutations,
            params.sparcc_p_threshold,
            params.sparcc_force,
        )

        SPARCC_CORRELATION(
            SPARCC_PREPARE.out.otu,
            params.sparcc_iterations,
            params.sparcc_exclusion_iterations,
            params.sparcc_permutations,
        )

        SPARCC_NETWORK(
            SPARCC_CORRELATION.out.correlation
                .join(SPARCC_CORRELATION.out.pvalues)
                .join(SPARCC_PREPARE.out.taxa),
            params.sparcc_permutations,
            params.sparcc_min_correlation,
            params.sparcc_p_threshold,
        )
        ch_multiqc_files = ch_multiqc_files.mix(SPARCC_NETWORK.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // MODULE: Ask the LLM to narrate the combined taxonomic profile.
    //
    // This runs BEFORE MultiQC on purpose: its `*_mqc.html` becomes a
    // MultiQC custom-content section, so the summary lands at the top of the
    // report the user already opens instead of in a file they never see.
    // Running here also keeps it off the critical path of every other AI
    // call - see subworkflows/local/ai_annotate_reports on why LLM calls
    // are strictly serialised.
    //
    if (ai_options.contains('taxonomy') && !params.skip_kraken2) {
        // The filtered tables, deliberately: the unfiltered long tail of
        // single-read taxa is exactly the material an LLM would narrate as
        // if it meant something.
        def ch_ai_tables = TAXONOMY_KRAKEN2_BRACKEN.out.report_combined_filtered
            .mix(TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined_filtered)
            .map { _meta, table -> table }
            .collect(sort: true)
            .map { tables -> [[id: 'reanatax_taxonomy'], tables] }

        LLM_INSIGHT(
            ch_ai_tables,
            params.llm_endpoint,
            params.llm_model,
            params.llm_api_key,
        )
        ch_multiqc_files = ch_multiqc_files.mix(LLM_INSIGHT.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // SUBWORKFLOW: Which taxa differ between two groups.
    //
    // Runs on the UNFILTERED combined Bracken table on purpose: ALDEx2 and
    // ANCOM-BC2 apply their own prevalence and abundance filters, and feeding
    // them a table already thinned by --min_rel_abundance would filter twice
    // with two different rules.
    //
    if (params.da_metadata && !params.skip_bracken && !params.skip_kraken2) {
        def da_missing = ['da_grouping', 'da_sample_group', 'da_ref_group'].findAll { key -> !params[key] }
        if (da_missing) {
            error("--da_metadata was given, so ${da_missing.collect { key -> '--' + key }.join(', ')} must be too.")
        }
        DIFFERENTIAL_ABUNDANCE(
            TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined,
            params.da_metadata,
            daMethods(params.da_method),
        )
        ch_multiqc_files = ch_multiqc_files.mix(DIFFERENTIAL_ABUNDANCE.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: The host half of the library.
    //
    // featureCounts already ran under --quantify_host; without this its output
    // was published per sample and read by nothing. The contrast is the same
    // --da_* design the microbial side uses, so one metadata file describes both.
    //
    // Two ways to get a host count matrix. --host_transcripts pseudoaligns the
    // TRIMMED reads against a transcriptome, which is independent of the
    // depletion and of --quantify_host; otherwise the matrix is merged from the
    // featureCounts tables the host BAM already produced.
    def ch_host_matrix = channel.empty()
    def host_quantified = false

    if (params.host_transcripts) {
        QUANTIFY_KALLISTO(ch_trimmed_reads, params.host_transcripts, params.host_tx2gene)
        ch_host_matrix = QUANTIFY_KALLISTO.out.counts
        host_quantified = true
        ch_multiqc_files = ch_multiqc_files.mix(QUANTIFY_KALLISTO.out.multiqc_files)
    }
    else if (params.quantify_host) {
        HOSTCOUNTS_MERGE(
            ch_host_counts
                .map { _meta, table -> table }
                .collect(sort: true)
                .map { tables -> [[id: 'reanatax'], tables] }
        )
        ch_host_matrix = HOSTCOUNTS_MERGE.out.counts
        host_quantified = true
    }

    if (host_quantified && params.da_metadata && params.host_de_method) {
        HOST_EXPRESSION(
            ch_host_matrix,
            TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined,
            params.da_metadata,
            hostDeMethods(params.host_de_method),
            params.host_microbe_correlation,
        )
        ch_multiqc_files = ch_multiqc_files.mix(HOST_EXPRESSION.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: one taxon, followed to its own differential expression.
    //
    if (params.target_taxid && !params.skip_kraken2) {
        TARGETED_TAXON(
            ch_nonhost_reads,
            TAXONOMY_KRAKEN2_BRACKEN.out.classifiedreads,
            TAXONOMY_KRAKEN2_BRACKEN.out.report,
            params.target_taxid,
            targetReference(),
            params.target_gtf,
            params.target_transcripts,
            params.target_tx2gene,
            hostDeMethods(params.target_de_method),
            params.da_metadata,
            params.ncbi_group,
        )
        ch_multiqc_files = ch_multiqc_files.mix(TARGETED_TAXON.out.multiqc_files)
    }

    //
    // MODULE: Alpha diversity, and what explains the variation between samples.
    //
    // On the FILTERED table, unlike differential abundance above, and for the
    // opposite reason: ALDEx2 and ANCOM-BC2 filter internally, diversity
    // indices do not. Observed richness on an unfiltered Bracken table is very
    // largely a count of database artefacts, and Shannon inherits that.
    //
    // Metadata is optional here: without it there is still alpha diversity and
    // an ordination, just nothing to partition the variance by. It falls back
    // to --da_metadata so one file can serve both when they agree.
    //
    if (params.run_diversity && !params.skip_bracken && !params.skip_kraken2) {
        def diversity_metadata = params.diversity_metadata ?: params.da_metadata
        DIVERSITY(
            TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined_filtered,
            diversity_metadata ? file(diversity_metadata, checkIfExists: true) : [],
            params.diversity_permutations,
        )
        ch_multiqc_files = ch_multiqc_files.mix(DIVERSITY.out.mqc.map { _meta, mqc -> mqc })
    }

    //
    // MODULE: One archive of the taxonomic tables, named so the exploreMetaTax
    // Shiny app can load them.
    //
    // Not a tar of the results directory: the app filters an uploaded folder by
    // FILENAME against per-format regexes, and several of the names published
    // here either match nothing or match the wrong format. The helper renames
    // them and refuses to seal an archive that would load wrongly - see
    // bin/exploremetatax_bundle.py.
    //
    if (params.exploremetatax_bundle) {
        // Each category is gathered into one sorted list so the process runs
        // once for the cohort; `ifEmpty` covers the steps that did not run.
        EXPLOREMETATAX_BUNDLE(
            [id: 'reanatax'],
            ch_kraken2_report.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_krakenuniq_report.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_bracken.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_bracken_combined.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_metaphlan_profile.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_metaphlan_merged.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            ch_humann_tables.map { _meta, file_ -> file_ }.flatten().collect(sort: true).ifEmpty([]),
            params.da_metadata ? file(params.da_metadata, checkIfExists: true) : [],
        )
    }

    //
    // MODULE: One table saying what happened to every read, and how much host
    // survived depletion. Everything it needs already exists; nothing joins it.
    //
    if (!params.skip_read_accounting) {
        // Each of the three lists is wrapped in a single-element tuple before
        // combining: `combine` concatenates tuples, so combining the bare lists
        // would flatten all three into one long tuple instead of nesting them.
        def ch_accounting_fastp = FASTQ_QC_TRIM.out.fastp_json
            .map { _meta, json -> json }
            .collect(sort: true)
            .map { files -> [files] }
            .ifEmpty([[]])
        def ch_accounting_sortmerna = ch_sortmerna_logs
            .map { _meta, log_file -> log_file }
            .collect(sort: true)
            .map { files -> [files] }
            .ifEmpty([[]])
        def ch_accounting_hisat2 = ch_hisat2_summaries
            .map { _meta, log_file -> log_file }
            .collect(sort: true)
            .map { files -> [files] }
            .ifEmpty([[]])
        // Either classifier's report answers the same three questions, and
        // read_kraken2() reads both layouts. Kraken2 wins when both ran, so a
        // run with both does not count the same reads twice.
        def ch_accounting_reports = params.skip_kraken2 ? ch_krakenuniq_report : ch_kraken2_report
        def ch_accounting_kraken2 = ch_accounting_reports
            .map { _meta, report -> report }
            .collect(sort: true)
            .map { files -> [files] }
            .ifEmpty([[]])

        READ_ACCOUNTING(
            ch_accounting_fastp
                .combine(ch_accounting_sortmerna)
                .combine(ch_accounting_hisat2)
                .combine(ch_accounting_kraken2)
                .map { fastp, sortmerna, hisat2, kraken2 ->
                    [[id: 'reanatax'], fastp, sortmerna, hisat2, kraken2]
                },
            params.host_carryover_taxid,
        )
        ch_multiqc_files = ch_multiqc_files.mix(READ_ACCOUNTING.out.mqc.map { _meta, mqc -> mqc }.flatten())
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'reanatax_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    def ch_multiqc_report = channel.empty()
    def ch_multiqc_data = channel.empty()

    if (!params.skip_multiqc) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
        def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
        ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        def ch_multiqc_custom_methods_description = multiqc_methods_description
            ? file(multiqc_methods_description, checkIfExists: true)
            : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
        ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        MULTIQC(
            ch_multiqc_files.flatten().collect().map { files ->
                [
                    [id: 'reanatax'],
                    files,
                    multiqc_config
                        ? file(multiqc_config, checkIfExists: true)
                        : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                    multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                    [],
                    [],
                ]
            }
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_multiqc_data = MULTIQC.out.data
    }

    //
    // SUBWORKFLOW: AI annotation of the finished reports.
    //
    // Entered whenever an endpoint is configured, even if no annotation was
    // selected, because it is also where the endpoint gets scrubbed out of the
    // MultiQC report (the run's parameter summary contains it, whether or not
    // MultiQC's own AI feature was used).
    //
    if (params.llm_endpoint) {
        AI_ANNOTATE_REPORTS(
            ch_multiqc_report,
            ch_multiqc_data,
            ch_qualimap,
            params.llm_endpoint,
            params.llm_model,
            params.llm_api_key,
            ai_options.contains('multiqc'),
            ai_options.contains('qualimap'),
        )
        // The annotated, redacted report is the one to link from the completion
        // email - the raw one is not even published in this case.
        ch_multiqc_report = AI_ANNOTATE_REPORTS.out.multiqc_report
    }

    emit:
    multiqc_report = ch_multiqc_report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

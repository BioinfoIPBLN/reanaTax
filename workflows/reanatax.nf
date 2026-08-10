/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_DOWNLOAD_FASTQDL    } from '../subworkflows/local/fastq_download_fastqdl'
include { PREPARE_HOST_REFERENCE    } from '../subworkflows/local/prepare_host_reference'
include { FASTQ_QC_TRIM             } from '../subworkflows/local/fastq_qc_trim'
include { PREPARE_HOST_REFERENCE as PREPARE_HOST_REFERENCE_FIRST } from '../subworkflows/local/prepare_host_reference'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FIRST } from '../subworkflows/local/host_depletion_hisat2'
include { HOST_DEPLETION_HISAT2 as HOST_DEPLETION_FINAL } from '../subworkflows/local/host_depletion_hisat2'
include { TAXONOMY_KRAKEN2_BRACKEN  } from '../subworkflows/local/taxonomy_kraken2_bracken'
include { FUNCTIONAL_HUMANN         } from '../subworkflows/local/functional_humann'
include { AI_ANNOTATE_REPORTS       } from '../subworkflows/local/ai_annotate_reports'
include { LLM_INSIGHT               } from '../modules/local/llm/insight/main'
include { READ_ACCOUNTING           } from '../modules/local/read/accounting/main'
include { POLYA_CARRYOVER           } from '../modules/local/polya/carryover/main'
include { POLYA_MERGE               } from '../modules/local/polya/merge/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { aiInsightOptions          } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { hostReferences            } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { meanReadLength            } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { brackenDistributions      } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { resolveBrackenReadLength  } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'

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
    // SUBWORKFLOW: Split host from non-host reads
    //
    def ch_nonhost_reads = FASTQ_QC_TRIM.out.reads
    def ch_qualimap = channel.empty()
    def ch_host_counts = channel.empty()
    def ch_hisat2_summaries = channel.empty()

    // Up to two references, applied in turn. The reads that survive the first
    // pass are the input to the second, so a read has to fail against BOTH
    // assemblies to be called non-host - which is the point of pairing GRCh38
    // with T2T-CHM13: the second assembly holds the centromeric, satellite and
    // structurally variant sequence that the first is missing, and those are
    // exactly the regions whose reads otherwise surface as spurious microbes.
    def host_refs = hostReferences()

    if (!params.skip_host_removal && host_refs) {
        def two_pass = host_refs.size() > 1
        // The first reference is the primary one: --gtf describes it, so it is
        // the one indexed splice-aware, QC'd with Qualimap and counted.
        def primary = host_refs.first()
        def secondary = two_pass ? host_refs[1] : null

        if (two_pass) {
            PREPARE_HOST_REFERENCE_FIRST(
                primary.kind == 'fasta' ? primary.value : null,
                primary.kind == 'index' ? primary.value : null,
                primary.kind == 'accession' ? primary.value : null,
                primary.kind == 'taxid' ? primary.value : null,
                params.ncbi_group,
                params.gtf,
            )
            HOST_DEPLETION_FIRST(
                FASTQ_QC_TRIM.out.reads,
                PREPARE_HOST_REFERENCE_FIRST.out.index,
                PREPARE_HOST_REFERENCE_FIRST.out.fasta,
                params.save_host_bam,
                params.skip_qualimap,
                params.qualimap_gff,
                params.quantify_host ? params.gtf : null,
            )
            ch_nonhost_reads = HOST_DEPLETION_FIRST.out.reads
            ch_qualimap = HOST_DEPLETION_FIRST.out.qualimap
            ch_host_counts = HOST_DEPLETION_FIRST.out.host_counts
            ch_hisat2_summaries = HOST_DEPLETION_FIRST.out.summary
            ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_FIRST.out.multiqc_files)
        }

        // The final pass is always the one whose leftovers get classified, so it
        // is always this alias - which is what lets conf/modules.config publish
        // the non-host FASTQs from exactly one place regardless of how many
        // references were given.
        def last = secondary ?: primary
        PREPARE_HOST_REFERENCE(
            last.kind == 'fasta' ? last.value : null,
            last.kind == 'index' ? last.value : null,
            last.kind == 'accession' ? last.value : null,
            last.kind == 'taxid' ? last.value : null,
            params.ncbi_group,
            // A GTF describes one assembly. Handing the primary reference's
            // annotation to a second, different assembly would build a
            // nonsense splice index, so the second pass never gets it.
            two_pass ? null : params.gtf,
        )
        HOST_DEPLETION_FINAL(
            two_pass ? HOST_DEPLETION_FIRST.out.reads : FASTQ_QC_TRIM.out.reads,
            PREPARE_HOST_REFERENCE.out.index,
            PREPARE_HOST_REFERENCE.out.fasta,
            params.save_host_bam,
            // Qualimap on the primary reference only; a second report of the
            // leftovers aligned to another assembly answers no question.
            two_pass ? true : params.skip_qualimap,
            params.qualimap_gff,
            two_pass ? null : (params.quantify_host ? params.gtf : null),
        )
        ch_nonhost_reads = HOST_DEPLETION_FINAL.out.reads
        ch_hisat2_summaries = ch_hisat2_summaries.mix(HOST_DEPLETION_FINAL.out.summary)
        ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_FINAL.out.multiqc_files)
        if (!two_pass) {
            ch_qualimap = HOST_DEPLETION_FINAL.out.qualimap
            ch_host_counts = HOST_DEPLETION_FINAL.out.host_counts
        }
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
            params.min_rel_abundance,
            params.min_samples,
        )
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_KRAKEN2_BRACKEN.out.multiqc_files)
        ch_kraken2_report = TAXONOMY_KRAKEN2_BRACKEN.out.report

        //
        // SUBWORKFLOW: What the community is doing, not just who is in it.
        //
        if (params.run_humann) {
            FUNCTIONAL_HUMANN(
                ch_nonhost_reads,
                TAXONOMY_KRAKEN2_BRACKEN.out.report,
                params.humann_nucleotide_db,
                params.humann_protein_db,
                params.humann_utility_db,
                params.humann_regroup,
                params.humann_renorm,
            )
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
        if (ai_options.contains('taxonomy')) {
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
    }

    //
    // MODULE: One table saying what happened to every read, and how much host
    // survived depletion. Everything it needs already exists; nothing joins it.
    //
    if (!params.skip_read_accounting) {
        READ_ACCOUNTING(
            FASTQ_QC_TRIM.out.fastp_json.map { _meta, json -> json }.collect(sort: true).ifEmpty([])
                .combine(ch_hisat2_summaries.map { _meta, log_file -> log_file }.collect(sort: true).ifEmpty([]))
                .combine(ch_kraken2_report.map { _meta, report -> report }.collect(sort: true).ifEmpty([]))
                .map { fastp, hisat2, kraken2 -> [[id: 'reanatax'], fastp, hisat2, kraken2] },
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_DOWNLOAD_FASTQDL    } from '../subworkflows/local/fastq_download_fastqdl'
include { PREPARE_HOST_REFERENCE    } from '../subworkflows/local/prepare_host_reference'
include { FASTQ_QC_TRIM             } from '../subworkflows/local/fastq_qc_trim'
include { HOST_DEPLETION_HISAT2     } from '../subworkflows/local/host_depletion_hisat2'
include { TAXONOMY_KRAKEN2_BRACKEN  } from '../subworkflows/local/taxonomy_kraken2_bracken'
include { AI_ANNOTATE_REPORTS       } from '../subworkflows/local/ai_annotate_reports'
include { LLM_INSIGHT               } from '../modules/local/llm/insight/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'
include { aiInsightOptions          } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'

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

    if (!params.skip_host_removal) {
        PREPARE_HOST_REFERENCE(
            params.fasta,
            params.hisat2_index,
            params.host_accession,
            params.host_taxid,
            params.ncbi_group,
            params.gtf,
        )

        HOST_DEPLETION_HISAT2(
            FASTQ_QC_TRIM.out.reads,
            PREPARE_HOST_REFERENCE.out.index,
            PREPARE_HOST_REFERENCE.out.fasta,
            params.save_host_bam,
            params.skip_qualimap,
            params.qualimap_gff,
        )
        ch_nonhost_reads = HOST_DEPLETION_HISAT2.out.reads
        ch_qualimap = HOST_DEPLETION_HISAT2.out.qualimap
        ch_multiqc_files = ch_multiqc_files.mix(HOST_DEPLETION_HISAT2.out.multiqc_files)
    }

    //
    // SUBWORKFLOW: Taxonomic classification of the non-host fraction
    //
    if (!params.skip_kraken2) {
        TAXONOMY_KRAKEN2_BRACKEN(
            ch_nonhost_reads,
            params.kraken2_db,
            params.bracken_db,
            params.kraken2_save_reads,
            params.kraken2_save_readclassifications,
            params.skip_bracken,
            params.skip_krona,
        )
        ch_multiqc_files = ch_multiqc_files.mix(TAXONOMY_KRAKEN2_BRACKEN.out.multiqc_files)

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
            def ch_ai_tables = TAXONOMY_KRAKEN2_BRACKEN.out.report_combined
                .mix(TAXONOMY_KRAKEN2_BRACKEN.out.bracken_combined)
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

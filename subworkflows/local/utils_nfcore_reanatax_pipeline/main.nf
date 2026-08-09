//
// Subworkflow with functionality specific to the BioinfoIPBLN/reanatax pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input_accessions PRJNA123456 --kraken2_db /path/to/db --host_accession GCF_000001405.40 --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Build the read channel from whichever of the three input routes was used.
    // Exactly one of them is set (enforced by validateInputParameters).
    //
    def ch_samplesheet = channel.empty()

    if (input) {
        ch_samplesheet = channel
            .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
            .map { meta, fastq_1, fastq_2 ->
                fastq_2
                    ? [ meta.id, meta + [ single_end: false ], [ fastq_1, fastq_2 ] ]
                    : [ meta.id, meta + [ single_end: true ], [ fastq_1 ] ]
            }
            .groupTuple()
            .map { samplesheet -> validateInputSamplesheet(samplesheet) }
            .map { meta, fastqs -> [ meta, fastqs.flatten() ] }
    }
    else if (params.input_dir) {
        ch_samplesheet = fastqsFromDir(params.input_dir, params.fastq_pattern, params.single_end)
    }

    //
    // Accessions are resolved to runs (and downloaded) inside the main workflow.
    //
    def ch_accessions = params.input_accessions
        ? channel.fromList(parseAccessions(params.input_accessions)).map { accession -> [ [ id: accession ], accession ] }
        : channel.empty()

    emit:
    samplesheet = ch_samplesheet
    accessions  = ch_accessions
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check that the combination of parameters actually describes a runnable job.
// These are cross-parameter rules that the JSON schema cannot express.
//
def validateInputParameters() {
    def routes = ['--input': params.input, '--input_dir': params.input_dir, '--input_accessions': params.input_accessions]
    def given = routes.findAll { _flag, value -> value }.keySet()
    if (given.size() != 1) {
        error(
            given.isEmpty()
                ? "No input given. Provide exactly one of --input (samplesheet), --input_dir (folder of FASTQs) or --input_accessions (ENA/SRA accessions)."
                : "Provide exactly one input route, but ${given.join(' and ')} were given."
        )
    }

    if (!params.skip_host_removal && !params.fasta && !params.hisat2_index && !params.host_accession && !params.host_taxid) {
        error("Host depletion is enabled but no host genome was given. Provide one of --fasta, --hisat2_index, --host_accession or --host_taxid, or disable the step with --skip_host_removal.")
    }

    if (params.host_accession && params.host_taxid) {
        error("--host_accession and --host_taxid are mutually exclusive.")
    }

    if (!params.skip_kraken2 && !params.kraken2_db) {
        error("--kraken2_db is required. Point it at a Kraken2 database directory (or .tar.gz), or disable classification with --skip_kraken2.")
    }

    // `--report-minimizer-data` adds two columns to the Kraken2 report, which
    // both Bracken and MultiQC's kraken parser choke on.
    if (params.kraken2_report_minimizer_data && !params.skip_bracken) {
        error("--kraken2_report_minimizer_data changes the Kraken2 report layout and cannot be read by Bracken. Add --skip_bracken, or drop the minimizer columns.")
    }

    if (!['run', 'experiment', 'sample'].contains(params.group_runs_by)) {
        error("--group_runs_by must be one of 'run', 'experiment' or 'sample' (got '${params.group_runs_by}').")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Turn a folder of FASTQ files into the same [ meta, [ reads ] ] channel shape a
// samplesheet produces.
//
// Mates are detected from the file name: a trailing `_1`/`_2`, `_R1`/`_R2` or
// `_R1_001`/`_R2_001` (also with `.` as the separator) marks a mate, anything
// else is treated as single-end. `--single_end` skips detection entirely, which
// is the escape hatch for sample names that legitimately end in `_1`.
//
def fastqsFromDir(input_dir, pattern, single_end) {
    def ch_fastq = channel.fromPath("${input_dir}/${pattern}", checkIfExists: true)

    if (single_end) {
        return ch_fastq.map { fastq ->
            [ [ id: fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, ''), single_end: true ], [ fastq ] ]
        }
    }

    return ch_fastq
        .map { fastq ->
            def stem = fastq.name.replaceFirst(/\.(fastq|fq)(\.gz)?$/, '')
            def matcher = stem =~ /^(.+)[._](?:R)?([12])(?:[._]001)?$/
            matcher.matches()
                ? [ matcher.group(1), matcher.group(2) as Integer, fastq ]
                : [ stem, 1, fastq ]
        }
        .groupTuple(by: 0)
        .map { id, mates, fastqs ->
            def ordered = [ mates, fastqs ].transpose().sort { entry -> entry[0] }
            if (ordered.size() > 2) {
                error("Found ${ordered.size()} FASTQ files for sample '${id}' in ${input_dir}: ${ordered.collect { entry -> entry[1].name }.join(', ')}. Use a samplesheet (--input) for datasets with multiple runs per sample.")
            }
            if (ordered.size() == 2 && ordered.collect { entry -> entry[0] } != [1, 2]) {
                error("Sample '${id}' in ${input_dir} does not have exactly one R1 and one R2 file: ${ordered.collect { entry -> entry[1].name }.join(', ')}.")
            }
            [ [ id: id, single_end: ordered.size() == 1 ], ordered.collect { entry -> entry[1] } ]
        }
}

//
// Accept accessions either as a comma-separated string or as a file with one
// accession per line (`#` comments and blank lines are ignored).
//
def parseAccessions(input_accessions) {
    def raw = input_accessions.toString()
    def candidate = file(raw)
    def text = !raw.contains(',') && candidate.exists() && !candidate.isDirectory() ? candidate.text : raw

    def accessions = text
        .split(/[,\r\n]+/)
        .collect { entry -> entry.trim() }
        .findAll { entry -> entry && !entry.startsWith('#') }
        .unique()

    if (accessions.isEmpty()) {
        error("No accessions found in --input_accessions '${input_accessions}'.")
    }

    // fastq-dl resolves ENA/SRA/DDBJ accessions only. GEO series (GSE/GSM) have
    // to be given as the linked SRA study, which is a common first mistake.
    def supported = ~/^(PRJ(EB|NA|DB)|[EDS]RP|SAM[DEN]|[EDS]RS|[EDS]RX|[EDS]RR)[A-Z]?[0-9]+$/
    def unsupported = accessions.findAll { accession -> !(accession ==~ supported) }
    if (unsupported) {
        def hint = unsupported.any { accession -> accession ==~ /^GS[EM][0-9]+$/ }
            ? " GEO accessions are not supported by fastq-dl: look up the linked SRA study (SRP.../PRJNA...) on the GEO page and use that instead."
            : ""
        error("Unsupported accession(s) for fastq-dl: ${unsupported.join(', ')}.${hint}")
    }

    return accessions
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    def citation_text = [
            "Tools used in the workflow included:",
            params.input_accessions ? "fastq-dl (Petit and Read 2023)," : "",
            params.skip_fastqc ? "" : "FastQC (Andrews 2010),",
            params.skip_trimming ? "" : "fastp (Chen et al. 2018),",
            params.skip_host_removal ? "" : "HISAT2 (Kim et al. 2019),",
            params.skip_host_removal || !params.save_host_bam ? "" : "SAMtools (Danecek et al. 2021),",
            params.skip_kraken2 ? "" : "Kraken2 (Wood et al. 2019),",
            params.skip_kraken2 || params.skip_bracken ? "" : "Bracken (Lu et al. 2017),",
            params.skip_kraken2 || params.skip_krona ? "" : "KronaTools (Ondov et al. 2011),",
            "and MultiQC (Ewels et al. 2016)",
            "."
        ].findAll { entry -> entry }.join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    def reference_text = [
            params.input_accessions ? "<li>Petit III, R. A., & Read, T. D. (2023). fastq-dl: efficiently download FASTQ files from SRA or ENA repositories. doi: 10.5281/zenodo.8051230</li>" : "",
            params.skip_fastqc ? "" : "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/</li>",
            params.skip_trimming ? "" : "<li>Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884-i890. doi: 10.1093/bioinformatics/bty560</li>",
            params.skip_host_removal ? "" : "<li>Kim, D., Paggi, J. M., Park, C., Bennett, C., & Salzberg, S. L. (2019). Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. Nature Biotechnology, 37(8), 907-915. doi: 10.1038/s41587-019-0201-4</li>",
            params.skip_host_removal || !params.save_host_bam ? "" : "<li>Danecek, P., Bonfield, J. K., Liddle, J., et al. (2021). Twelve years of SAMtools and BCFtools. GigaScience, 10(2), giab008. doi: 10.1093/gigascience/giab008</li>",
            params.skip_kraken2 ? "" : "<li>Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. Genome Biology, 20(1), 257. doi: 10.1186/s13059-019-1891-0</li>",
            params.skip_kraken2 || params.skip_bracken ? "" : "<li>Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. PeerJ Computer Science, 3, e104. doi: 10.7717/peerj-cs.104</li>",
            params.skip_kraken2 || params.skip_krona ? "" : "<li>Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. BMC Bioinformatics, 12, 385. doi: 10.1186/1471-2105-12-385</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047-3048. doi: 10.1093/bioinformatics/btw354</li>"
        ].findAll { entry -> entry }.join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()

    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

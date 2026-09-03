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
        //
        // The template dumps every parameter to pipeline_info/params_*.json, so
        // that file would carry the LLM endpoint and API key into the results
        // directory in plain text. Scrub it here, once, at the end of the run.
        // (The MultiQC report is handled by MULTIQC_AI; the parameter summary
        // never gets them in the first place, see validation.summary.hideParams.)
        //
        redactPipelineInfo(outdir)

        //
        // Stop the Kraken2 daemon. It runs outside a PID namespace so that it
        // survives the task that started it - that is what lets every later
        // task reuse the resident index instead of reloading 344 GB per sample
        // - and the flip side is that nothing reaps it when the run ends. Left
        // alone it sits on the node holding the entire index; a completed run
        // was measured still holding 316 GB an hour and a half later.
        //
        stopKraken2Daemon()

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
        // A failed run strands the daemon exactly as effectively as a
        // successful one, so this runs on both paths.
        stopKraken2Daemon()
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

    if (params.kraken2_use_daemon) {
        if (params.kraken2_memory_mapping) {
            error("--kraken2_use_daemon and --kraken2_memory_mapping are contradictory: the daemon exists to hold the index in RAM, memory-mapping exists to avoid doing so. Pick one.")
        }
        if (params.kraken2_db && (params.kraken2_db.toString().endsWith('.tar.gz') || params.kraken2_db.toString().endsWith('.tgz'))) {
            error("--kraken2_use_daemon needs an unpacked database directory: the daemon keys its resident index on the --db path, and an untarred copy lives at a different path in every task. Point --kraken2_db at the extracted directory.")
        }
        if (workflow.containerEngine == 'docker') {
            error("--kraken2_use_daemon cannot work under Docker: every task gets its own PID namespace and its own /tmp, so the daemon and the FIFOs it is addressed through do not survive the task that started them. Use -profile singularity/apptainer, or drop the flag.")
        }
        log.warn("--kraken2_use_daemon leaves a background process holding the index after this run finishes. Stop it with 'k2 clean --stop-daemon' on the execution node when you are done.")
    }

    if (params.single_cell) {
        def chemistry = scChemistry()
        if (params.input_accessions) {
            log.warn("--single_cell with --input_accessions: the archives rarely label which run is the barcode read, so check that fastq_1 really is the barcode+UMI read and fastq_2 the cDNA read. A samplesheet makes this explicit.")
        }
        if (!chemistry.cb_len || !chemistry.umi_len || !chemistry.umi_start) {
            error("--sc_chemistry '${params.sc_chemistry}' has no preset, so --sc_cb_len, --sc_umi_start and --sc_umi_len must all be given.")
        }
        if (!params.sc_whitelist) {
            error("--single_cell needs --sc_whitelist: the barcode whitelists ship with Cell Ranger, not with STAR. For ${params.sc_chemistry} that is '${chemistry.expected_whitelist}'. Pass --sc_whitelist None to accept every observed barcode uncorrected, which lets one sequencing error in a barcode create a new cell.")
        }
        if (!params.gtf && !params.star_index) {
            error("--single_cell needs --gtf: STARsolo assigns reads to genes at alignment time, so an index built without an annotation cannot produce a cell-by-gene matrix. Supply --gtf, or a --star_index that was built with one.")
        }
        if (params.skip_kraken2) {
            error("--single_cell with --skip_kraken2 leaves nothing to classify the non-host reads with, and the cell-by-taxon matrix is built from Kraken2's read-level output.")
        }
        if (!params.kraken2_save_readclassifications) {
            error("--single_cell needs Kraken2's read-level output to recover the cell barcode from each read name, which is only written with --kraken2_save_readclassifications.")
        }
        if (params.remove_rrna) {
            error("--remove_rrna cannot be combined with --single_cell. 10x chemistry is poly(A) capture, and bacterial mRNA is polyadenylated rarely and with short tails, so poly(A) capture of bacteria is effectively rRNA capture: CSI-Microbes measured 82-95% of captured bacterial reads as rRNA (plexWell 82%, 10x 5' 92%, 10x 3'v3 95%). Depleting rRNA here would delete most of the microbial signal before Kraken2 ever sees it. rRNA depletion is correct for total-RNA/ribo-depleted libraries and catastrophic for poly(A) ones.")
        }
    }

    if (params.sc_ambient) {
        if (!params.single_cell) {
            error("--sc_ambient needs --single_cell. It compares cell-containing droplets against empty ones, and only droplet chemistry has empty droplets - in plate-based data (--sc_plate_based) every well is a library and there is no ambient pool to compare against.")
        }
        if (params.sc_cell_filter.toString().toLowerCase() == 'none') {
            error("--sc_ambient cannot be combined with --sc_cell_filter None. It needs BOTH halves of STARsolo's Solo.out tree - filtered/ to name the cells and raw/ to name the empty droplets - and with no cell filter STARsolo writes only one matrix, so there is nothing to compare cells against.")
        }
    }

    if (params.sc_host_de) {
        if (!params.single_cell) {
            error("--sc_host_de needs --single_cell: the comparison is between cells of one library, and it reads the cell-by-gene matrix STARsolo produces. On plate-based data (--sc_plate_based) every cell is its own library, so the equivalent question is asked by the ordinary --host_de_method contrast with the cell type as the grouping variable.")
        }
        if (!params.sc_cell_metadata && !params.sc_host_de_force_pooled) {
            error("--sc_host_de needs --sc_cell_metadata. Which cell types carry a taxon is itself a result - it is what the enrichment step measures - so infection is not distributed at random over cell types, and a pooled infected-vs-uninfected test recovers the difference BETWEEN cell types and reports it as a response TO infection. Supply the annotations, or pass --sc_host_de_force_pooled to accept a pooled comparison; every row of the output is then stamped ALL_POOLED.")
        }
    }

    if (params.sc_plate_based) {
        if (params.single_cell) {
            error("--sc_plate_based and --single_cell are two different experiments, not two settings of one. In droplet data one library holds thousands of barcoded cells and needs STARsolo; in plate-based data every well is its own library and needs none of that machinery. Pick the one that matches the samplesheet.")
        }
        if (!params.sc_cell_metadata) {
            error("--sc_plate_based needs --sc_cell_metadata: with one library per cell, the pipeline has no other way to know which cells came from which patient or what type they are. The file must carry the cell id in the first column plus a sample/patient/donor/plate column - without the latter every cell becomes its own stratum and the enrichment test collapses to the pooled one it exists to avoid.")
        }
        if (params.skip_kraken2) {
            error("--sc_plate_based builds its cell-by-taxon matrix from the per-cell Kraken2 reports, so it cannot be combined with --skip_kraken2.")
        }
    }

    if (params.univec && !params.univec_source) {
        error("--univec needs --univec_source: it names the vector database to deplete against, and defaults to NCBI's UniVec_Core. Set it to a local FASTA to run offline.")
    }

    if (params.run_prism) {
        if (!params.prism_path) {
            error("--run_prism needs --prism_path, pointing at a prepared clone of sjdlabgroup/PRISM. The clone alone is not enough: it must also hold `genbank/` (the processed GenBank RDS files, distributed separately by the authors) and `sorted_accession_map.txt` (built once with blastdbcmd against the BLAST database). See docs/usage.md.")
        }
        def missing = [
            prism_blast_db: params.prism_blast_db,
            prism_star_genome_dir: params.prism_star_genome_dir,
            prism_minimap2_index: params.prism_minimap2_index,
        ].findAll { _name, value -> !value }.keySet()
        if (missing) {
            error("--run_prism also needs ${missing.collect { name -> '--' + name }.join(', ')}. PRISM shells out to BLAST, STAR and minimap2 and validates every path before it starts, so a missing one fails the task rather than degrading the result.")
        }
        if (!(params.prism_kraken_db ?: params.kraken2_db)) {
            error("--run_prism needs a Kraken2 database: give --prism_kraken_db, or --kraken2_db, which it falls back to.")
        }
        if (params.single_cell && !params.prism_barcode_only) {
            error("--run_prism with --single_cell needs --prism_barcode_only. On 10x chemistry read 1 is barcode plus UMI with no biological sequence in it, so without that flag PRISM classifies and BLASTs it as if it were cDNA. Setting it drops the barcode read before PRISM sees anything, so Kraken2, STAR, minimap2 and BLAST all work on the cDNA read alone.")
        }
        if (workflow.containerEngine in ['docker', 'singularity', 'apptainer', 'podman'] && !params.prism_container) {
            error("--run_prism under -profile ${workflow.containerEngine} needs --prism_container. PRISM is an R script that shells out to Kraken2, minimap2, STAR, BLAST and SeqKit, and no public image carries that set - so there is no default to fall back to. Build one from modules/local/prism/run/environment.yml, or use -profile conda, which builds that environment directly.")
        }
    }

    if (params.gene_diversity_filter) {
        if (!params.run_humann) {
            error("--gene_diversity_filter needs --run_humann. It asks whether a taxon's reads spread over its genome or pile onto one locus, and the only table in this pipeline that says so is HUMAnN's species-stratified gene-family output. Nothing else here records which gene a read landed on.")
        }
        if (params.skip_bracken) {
            error("--gene_diversity_filter needs Bracken: HUMAnN names a species by its clade name, and the combined Bracken table is what turns that name back into the taxid the abundance filter removes by.")
        }
    }

    if (params.pathseq_microbe_bwa_image) {
        def missing = [
            pathseq_microbe_dict: params.pathseq_microbe_dict,
            pathseq_taxonomy_db: params.pathseq_taxonomy_db,
        ].findAll { _name, value -> !value }.keySet()
        if (missing) {
            error("--pathseq_microbe_bwa_image also needs ${missing.collect { name -> '--' + name }.join(' and ')}. PathSeq scores a taxon from the alignments of its reads, so it needs the sequence dictionary to know what was aligned to and the taxonomy database to know where that sits on the tree; the index image alone cannot produce a result. Broad distributes all three prebuilt as one resource bundle.")
        }
        if ((params.pathseq_host_bwa_image ? 1 : 0) + (params.pathseq_host_kmers ? 1 : 0) == 1) {
            error("--pathseq_host_bwa_image and --pathseq_host_kmers go together: the k-mer file is what makes the alignment-based host filter affordable, and PathSeq's filter stage expects the pair. Give both to run PathSeq's own host filter as a second pass after HISAT2, or neither to skip it - which is reasonable here, because the reads reaching PathSeq are the ones host depletion could not place.")
        }
    }
    else if (params.pathseq_microbe_dict || params.pathseq_taxonomy_db || params.pathseq_host_bwa_image || params.pathseq_host_kmers) {
        error("PathSeq resources were given but --pathseq_microbe_bwa_image was not, so the PathSeq route is off and they would be ignored.")
    }

    if (params.decontam) {
        def decontam_metadata = params.decontam_metadata ?: params.da_metadata
        if (!decontam_metadata) {
            error("--decontam needs sample metadata naming the negative controls: give --decontam_metadata (or --da_metadata, which it falls back to).")
        }
        if (params.skip_kraken2 || params.skip_bracken) {
            error("--decontam is scored on the combined Bracken table, so it cannot be combined with --skip_kraken2 or --skip_bracken.")
        }
        def needs_neg = params.decontam_method in ['prevalence', 'combined', 'either', 'minimum']
        def needs_conc = params.decontam_method in ['frequency', 'combined', 'either', 'minimum']
        if (needs_neg && !params.decontam_neg_column) {
            error("--decontam_method '${params.decontam_method}' compares against negative controls, so --decontam_neg_column must name the metadata column that marks them.")
        }
        if (needs_conc && !params.decontam_conc_column) {
            error("--decontam_method '${params.decontam_method}' regresses abundance on DNA concentration, so --decontam_conc_column must name the metadata column holding it.")
        }
    }

    if (params.run_diversity && (params.skip_kraken2 || params.skip_bracken)) {
        error("--run_diversity works from the combined Bracken table, so it cannot be combined with --skip_kraken2 or --skip_bracken.")
    }

    if (params.host_kmer_filter) {
        if (params.skip_kraken2) {
            error("--host_kmer_filter reads Kraken2's per-read output, so it cannot be combined with --skip_kraken2.")
        }
        if (!params.kraken2_save_readclassifications) {
            error("--host_kmer_filter needs Kraken2's read-level output, which is only written with --kraken2_save_readclassifications. Add that flag, or drop --host_kmer_filter.")
        }
        if (!(params.host_kmer_taxid ?: params.host_carryover_taxid)) {
            error("--host_kmer_filter needs a host taxid: set --host_kmer_taxid, or leave --host_carryover_taxid at a value that names your host.")
        }
    }

    if (params.remove_rrna && !params.single_cell) {
        log.warn("--remove_rrna is correct for total-RNA / ribo-depleted libraries and wrong for poly(A)-selected ones, where most captured bacterial reads are rRNA (82-95% in CSI-Microbes' measurements). On a poly(A) library this deletes most of the microbial signal. Check the library type.")
    }

    if (params.minimizer_correlation && !params.minimizer_filter) {
        error("--minimizer_correlation is part of the minimizer filter; add --minimizer_filter to enable it.")
    }

    if (params.minimizer_filter && !params.skip_kraken2 && !params.kraken2_report_minimizer_data) {
        error("--minimizer_filter needs the distinct-minimizer columns, which Kraken2 only writes with --kraken2_report_minimizer_data. Add that flag, or run the filter on a KrakenUniq route alone (--krakenuniq_db --skip_kraken2).")
    }

    if (params.remove_rrna && !params.sortmerna_db) {
        error("--remove_rrna needs rRNA references: give --sortmerna_db (e.g. the smr_*_db.fasta shipped in SortMeRNA's database.tar.gz), or drop --remove_rrna.")
    }

    def host_refs = hostReferences()

    if (!params.skip_host_removal && host_refs.isEmpty()) {
        error("Host depletion is enabled but no host genome was given. Provide one of --fasta, --hisat2_index, --host_accession or --host_taxid, or disable the step with --skip_host_removal.")
    }

    // The limit is how many depletion passes the workflow declares, not a
    // property of the data. Nextflow needs one alias per invocation of a
    // subworkflow, so the passes are named at compile time and their number is
    // fixed; four covers the cases that actually arise - two assemblies of one
    // species (GRCh38 + T2T-CHM13), a second organism the library genuinely
    // contains (a blood meal, a graft, a co-cultured cell line), and a
    // synthetic reference such as UniVec or PhiX.
    if (host_refs.size() > maxHostPasses()) {
        error("At most ${maxHostPasses()} host references are supported, but ${host_refs.size()} were given (${host_refs.collect { ref -> ref.value }.join(', ')}). Reads are depleted against each in turn. To go higher, add another HOST_DEPLETION_HISAT2 alias in workflows/reanatax.nf, its PREPARE_HOST_REFERENCE alias, its naming rules in conf/modules.config, and raise maxHostPasses().")
    }

    // Measured, not assumed: with `--very-sensitive -k 1` on the command line
    // HISAT2 still reported up to 50 alignments per read (identical NH-tag
    // distributions to a `-k 50` run, and identical alignment rates). The
    // preset wins over a later -k, so --hisat2_max_alignments is silently
    // inert while a preset is in hisat2_args. Expand the preset by hand to
    // make it bite: --hisat2_args '--bowtie2-dp 2 --score-min L,0,-1'.
    if (params.hisat2_max_alignments && params.hisat2_args =~ /--(very-)?(sensitive|fast)/) {
        log.warn("--hisat2_max_alignments ${params.hisat2_max_alignments} has no effect: the preset in --hisat2_args ('${params.hisat2_args}') overrides -k. See docs/usage.md.")
    }

    if (params.quantify_host && !params.gtf) {
        error("--quantify_host needs an annotation: add --gtf. Counts are taken from the host BAM of the first reference given.")
    }

    // Validate the method names up front rather than at the point of use, so a
    // typo fails in seconds instead of after the classification has run.
    hostDeMethods(params.host_de_method)

    if (params.host_de_method && params.da_metadata && !params.quantify_host) {
        log.warn("--host_de_method is set and --da_metadata was given, but --quantify_host is not: there are no host gene counts to test. Add --quantify_host (which also needs --gtf), or set --host_de_method null.")
    }

    if (params.host_microbe_correlation && !(params.quantify_host && params.da_metadata)) {
        error("--host_microbe_correlation needs host gene counts and the sample metadata: add --quantify_host (with --gtf) and --da_metadata.")
    }

    if (params.target_taxid) {
        hostDeMethods(params.target_de_method)
        if (!params.target_reference) {
            error("--target_taxid needs --target_reference: extracting a taxon's reads is only useful if they can then be aligned to that organism's genome. Give an accession (GCF_*/GCA_*), a taxid, a FASTA or a HISAT2 index.")
        }
        if (!params.kraken2_save_readclassifications) {
            error("--target_taxid needs --kraken2_save_readclassifications: the per-read output is the record of which read Kraken2 put on which taxon, and without it there is nothing to extract from.")
        }
        if (params.skip_kraken2) {
            error("--target_taxid needs Kraken2; it cannot be combined with --skip_kraken2.")
        }
        if (params.target_de_method && !params.target_gtf) {
            log.warn("--target_taxid is set without --target_gtf: the branch will align the extracted reads and stop at the BAM, because gene counts need an annotation.")
        }
        if (params.target_de_method && params.target_gtf && !params.da_metadata) {
            log.warn("--target_taxid will count genes but not test them: differential expression needs --da_metadata and a --da_grouping contrast.")
        }
        if (params.target_include_parents) {
            log.warn("--target_include_parents also takes reads Kraken2 could only place ABOVE the target taxon. Those are consistent with the target but are not evidence for it, and they will be quantified as if they were.")
        }
    }

    if (params.host_transcripts && params.quantify_host) {
        log.warn("--host_transcripts and --quantify_host both give host gene counts; kallisto's are used and featureCounts' are published but not tested. Drop --quantify_host unless you want both tables.")
    }

    if (params.host_tx2gene && !params.host_transcripts) {
        error("--host_tx2gene has nothing to map without --host_transcripts.")
    }

    if (params.target_tx2gene && !params.target_transcripts) {
        error("--target_tx2gene has nothing to map without --target_transcripts.")
    }

    if (params.target_transcripts && params.target_gtf) {
        log.warn("--target_transcripts is set, so the targeted taxon is quantified by pseudoalignment and --target_gtf is not used. Drop one of them.")
    }

    if (params.run_humann && !(params.humann_nucleotide_db && params.humann_protein_db)) {
        error("--run_humann requires both --humann_nucleotide_db (ChocoPhlAn) and --humann_protein_db (UniRef).")
    }

    if (params.run_humann && params.skip_kraken2) {
        error("--run_humann needs the Kraken2 report to build HUMAnN's taxonomic profile; it cannot be combined with --skip_kraken2.")
    }

    if (params.skip_kraken2 && !params.krakenuniq_db && !params.run_metaphlan) {
        log.warn("--skip_kraken2 with neither --krakenuniq_db nor --run_metaphlan leaves nothing to classify the non-host reads with.")
    }

    if (!params.skip_kraken2 && !params.kraken2_db) {
        error("--kraken2_db is required. Point it at a Kraken2 database directory (or .tar.gz), or disable classification with --skip_kraken2.")
    }

    // `--report-minimizer-data` inserts two columns in the MIDDLE of the Kraken2
    // report, and every consumer here was checked against that: Bracken 3.1,
    // combine_kreports.py, kreport2mpa.py and kreport2krona.py all index the
    // rank/taxid/name from the END of the row and the counts from the front, so
    // the insertion passes between them untouched; MultiQC's kraken module
    // unpacks the 8-field layout explicitly and plots the duplication it adds.
    // read_accounting.py and hasClassifiedReads() key on the taxid, which is
    // stable across all three layouts this pipeline can produce. An earlier
    // guard here refused the combination on the assumption Bracken could not
    // read it; that was wrong, and removing it is what makes --minimizer_filter
    // usable without --skip_bracken.

    if (!['run', 'experiment', 'sample'].contains(params.group_runs_by)) {
        error("--group_runs_by must be one of 'run', 'experiment' or 'sample' (got '${params.group_runs_by}').")
    }

    // Parses and therefore validates --ai_insights; also warns about the
    // combinations that silently do nothing.
    def ai_options = aiInsightOptions(params.ai_insights, params.llm_endpoint)

    if (params.llm_endpoint && !params.llm_model) {
        error("--llm_endpoint was given but --llm_model was not. The endpoint has no way to know which model to serve.")
    }

    if (!params.llm_endpoint && (params.llm_model || params.multiqc_ai_builtin)) {
        log.warn("An LLM option was set but --llm_endpoint was not, so every AI annotation is disabled.")
    }

    if (params.multiqc_ai_builtin && !params.llm_endpoint) {
        error("--multiqc_ai_builtin requires --llm_endpoint (MultiQC's own AI feature needs an endpoint to call).")
    }

    if (ai_options.contains('qualimap') && (params.skip_qualimap || params.skip_host_removal || !params.save_host_bam)) {
        log.warn("--ai_insights includes 'qualimap' but no Qualimap report will be produced (--skip_qualimap / --skip_host_removal / --save_host_bam false).")
    }

    if (ai_options.contains('taxonomy') && params.skip_kraken2) {
        log.warn("--ai_insights includes 'taxonomy' but classification is disabled with --skip_kraken2, so there is nothing to summarise.")
    }

    if (!(params.alignment_output_format in ['bam', 'cram'])) {
        error("--alignment_output_format must be 'bam' or 'cram', not '${params.alignment_output_format}'.")
    }
    if (params.alignment_output_format == 'cram' && params.save_host_bam && !params.skip_qualimap) {
        error("--alignment_output_format cram cannot run with Qualimap: qualimap bamqc takes -bam and cannot open a CRAM. Add --skip_qualimap to keep CRAM - samtools stats, flagstat and idxstats all read CRAM and still run - or leave the format as bam.")
    }
    if (params.alignment_output_format == 'cram' && !params.host) {
        error("--alignment_output_format cram needs the reference it is encoded against: CRAM stores differences from a reference rather than the sequence itself, so the file is unreadable without it. Pass --host, or use bam.")
    }
    // Subread reads SAM and BAM. featureCounts on a CRAM fails per sample, at
    // the very end of a run, having already paid for the alignment.
    if (params.alignment_output_format == 'cram' && params.quantify_host) {
        error("--alignment_output_format cram cannot run with --quantify_host: featureCounts (Subread) reads SAM and BAM only. Drop --quantify_host to keep CRAM, or leave the format as bam.")
    }
    // The chunk BAMs are deleted from inside HISAT2_MERGECHUNKS, which only
    // runs when the library was chunked. Saying so is better than leaving the
    // user to infer that half the flag did nothing.
    if (params.cleanup_intermediates && !params.hisat2_chunk_size && !params.input_accessions) {
        log.warn("--cleanup_intermediates has nothing to delete on this run: the chunk BAMs it removes only exist with --hisat2_chunk_size, and the raw FASTQs it removes only exist for a downloaded cohort. A samplesheet's own FASTQs are never touched.")
    }
    if (params.cleanup_intermediates && params.skip_trimming) {
        log.warn("--cleanup_intermediates will not remove the downloaded FASTQs with --skip_trimming: they are then the working read set, and everything downstream is still reading them. The chunk BAMs are still removed.")
    }
    if (params.compression_level != null && (!(params.compression_level instanceof Integer) || params.compression_level < 0 || params.compression_level > 9)) {
        error("--compression_level must be an integer from 0 to 9, not '${params.compression_level}'. Leave it unset to keep every tool at its own default.")
    }

    // The negative-control filter. Its settings are validated in
    // controlSettings(); what belongs here is the interaction with the rest of
    // the run, which that function cannot see.
    if (params.negative_controls && params.skip_kraken2) {
        error("--negative_controls has nothing to filter with --skip_kraken2: the control levels are read off the combined Kraken2 and Bracken tables.")
    }
    if (params.negative_controls && (params.control_ratio as double) <= 0) {
        error("--control_ratio must be greater than 0. A ratio of 0 keeps every cell that has any reads at all, which is the same as not running the filter.")
    }
    // A control level is an average of the blanks, and an average of one blank
    // is that blank. Salter et al. (BMC Biol 2014) is the standard warning:
    // contamination varies substantially between extractions, so a single
    // control measures one draw from that variation and calling it the
    // background sets the threshold wherever that draw happened to land.
    if (params.negative_controls && resolveNegativeControls(params.negative_controls).split(',').size() < 3) {
        log.warn("--negative_controls names fewer than 3 libraries. The control level is then an estimate from almost no data, and contamination varies enough between extractions that one blank can be several-fold off in either direction. Consider --control_statistic max, which at least fails in the conservative direction.")
    }
    if ((params.prevalence_filter as double) > 0 && (params.prevalence_filter as double) < 0.5) {
        log.warn("--prevalence_filter ${params.prevalence_filter} removes any taxon present in ${(100 * (params.prevalence_filter as double)) as int}% of libraries. That is a statement that nothing real is shared by that many samples, which is false for most cohorts. It is meant to be used near 1.0.")
    }
}

//
// Replace the LLM endpoint and API key with a placeholder everywhere in
// pipeline_info/. The model name is deliberately left alone - it is provenance,
// not a secret. No-op when no endpoint was given.
//
// The parameter dump (params_*.json) is the file this is really for, but every
// text artifact in the directory is scrubbed rather than just that one, so a
// future template change cannot silently reintroduce a leak.
//
// Two places stay out of reach. Nextflow's work directory keeps every
// `.command.sh`, as it does for any parameter. And execution_report_*.html -
// which quotes both the run command line and each task script - is written by
// Nextflow's own observer AFTER this hook, so it still shows whatever was typed
// on the command line. Passing the credentials through `-params-file` or a
// private `-c` config keeps them out of that too.
//
def redactPipelineInfo(outdir) {
    if (!params.llm_endpoint) {
        return
    }
    def secrets = [params.llm_endpoint, params.llm_api_key]
        .findAll { secret -> secret && secret != 'dummy' && secret.toString().size() >= 4 }
        .collect { secret -> secret.toString() }
        .sort { secret -> -secret.size() }

    def info_dir = file("${outdir}/pipeline_info")
    if (secrets.isEmpty() || !info_dir.exists()) {
        return
    }

    def scrubbed = []
    info_dir.list().findAll { name -> name ==~ /(?i).*\.(html|json|txt|tsv|csv|log|yml|yaml)$/ }.each { name ->
        def target = info_dir.resolve(name)
        try {
            def text = target.text
            def redacted = secrets.inject(text) { acc, secret -> acc.replace(secret, '[redacted]') }
            if (redacted != text) {
                target.text = redacted
                scrubbed << name
            }
        }
        catch (Exception e) {
            log.warn("Could not scrub the LLM endpoint from ${outdir}/pipeline_info/${name}: ${e.message}")
        }
    }
    if (scrubbed) {
        log.info("Scrubbed the LLM endpoint/API key from ${scrubbed.size()} file(s) in ${outdir}/pipeline_info/.")
    }
}

//
// The host references to deplete against, in the order they will be applied.
//
// Each of --hisat2_index / --fasta / --host_accession / --host_taxid accepts a
// comma-separated list, and they can be mixed: `--fasta hg38.fa --host_accession
// GCA_009914755.4` is one local assembly plus one fetched from NCBI, which is
// exactly the GRCh38 + T2T-CHM13 pairing that Monteleone et al. use to stop
// centromeric, satellite and structurally variant reads from being mistaken for
// microbes.
//
// The FIRST reference is the primary one: it is the one --gtf describes, so it
// is the one indexed splice-aware, QC'd with Qualimap and counted for host
// expression.
//
//
// The four minimizer-evidence thresholds, gathered so the taxonomy subworkflows
// take one argument rather than four positional numbers it would be easy to
// transpose.
//
def minimizerThresholds() {
    return [
        min_reads: params.minimizer_min_reads,
        min_distinct: params.minimizer_min_distinct,
        max_duplication: params.minimizer_max_duplication,
        min_coverage: params.minimizer_min_coverage,
        distinct_scale: params.minimizer_distinct_scale,
    ]
}

//
// The host-k-mer filter's settings. The taxid defaults to whatever the report
// already treats as host leakage, so the two agree unless deliberately split -
// --host_kmer_taxid also takes a comma-separated list, to add the host's genus
// when Kraken2 can only place its k-mers that high.
//
def hostKmerSettings() {
    return [
        taxid: params.host_kmer_taxid ?: params.host_carryover_taxid,
        max_fraction: params.host_kmer_max_fraction,
        min_reads: params.host_kmer_min_reads,
    ]
}

//
// decontam's settings, or null when it is off. Returned as one map so the
// subworkflow takes a single argument and "is decontam on" is the same test as
// "are its settings present".
//
def decontamSettings() {
    if (!params.decontam) {
        return null
    }
    return [
        metadata: params.decontam_metadata ?: params.da_metadata,
        neg_column: params.decontam_neg_column,
        neg_value: params.decontam_neg_value,
        conc_column: params.decontam_conc_column,
        method: params.decontam_method,
        threshold: params.decontam_threshold,
        batch_column: params.decontam_batch_column,
        batch_combine: params.decontam_batch_combine,
    ]
}

//
// The shuffled-read control's settings, or null when it is off.
//
def shuffleSettings() {
    if (!params.shuffle_control) {
        return null
    }
    def known = ['dinuc', 'mono', 'reverse']
    if (!known.contains(params.shuffle_method)) {
        error("--shuffle_method: '${params.shuffle_method}' is not one of ${known.join(', ')}.")
    }
    return [
        method: params.shuffle_method,
        seed: params.shuffle_seed,
        max_reads: params.shuffle_reads,
        max_ratio: params.shuffle_max_ratio,
        min_reads: params.shuffle_min_reads,
    ]
}

//
// The negative-control filter's settings, or null when neither half is on.
//
// --negative_controls is resolved to a plain comma-separated list of sample IDs
// here rather than in the module, because the four accepted forms are a
// usability affordance and not something a container should have to know about:
//
//   SRX1,SRX2          the IDs themselves
//   controls.txt       a file of one ID per line
//   meta.tsv:column    every row whose `column` is truthy - true/yes/1/control
//   meta.tsv:column:x  every row whose `column` equals x exactly
//
// The metadata forms exist so a cohort that already declares its blanks for
// --decontam does not have to declare them twice. The ID column is the first
// one, matching --da_metadata and bin/decontam_filter.R.
//
def controlSettings() {
    def prevalence = params.prevalence_filter as double
    if (!params.negative_controls && prevalence <= 0) {
        return null
    }
    def known = ['mean', 'median', 'max']
    if (!known.contains(params.control_statistic)) {
        error("--control_statistic: '${params.control_statistic}' is not one of ${known.join(', ')}.")
    }
    if (prevalence < 0 || prevalence > 1) {
        error("--prevalence_filter is a FRACTION of libraries: ${params.prevalence_filter} is outside 0-1. 1.0 means 'present in every library'.")
    }
    return [
        controls: resolveNegativeControls(params.negative_controls),
        ratio: params.control_ratio,
        statistic: params.control_statistic,
        min_reads: params.control_min_reads,
        floor_reads: params.control_floor_reads,
        prevalence: prevalence,
        prevalence_min_reads: params.prevalence_min_reads,
    ]
}

//
// The four forms above, collapsed to one comma-separated string.
//
def resolveNegativeControls(spec) {
    if (!spec) {
        return ''
    }
    def text = spec.toString().trim()

    // The metadata forms, recognised by the part before the first ':' being a
    // file that exists. Tested that way round so a path is never mistaken for
    // a column separator, and a list of IDs never for a path.
    def parts = text.split(':') as List
    if (parts.size() >= 2 && file(parts[0]).exists()) {
        def wanted = parts.size() >= 3 ? parts[2] : null
        def truthy = ['true', 'yes', 'y', '1', 'control', 'blank', 'negative']
        def lines = file(parts[0]).readLines().findAll { line -> line.trim() && !line.startsWith('#') }
        if (lines.size() < 2) {
            error("--negative_controls: ${parts[0]} has no data rows.")
        }
        def header = lines[0].split('\t') as List
        def column = header.findIndexOf { field -> field.trim() == parts[1] }
        if (column < 0) {
            error("--negative_controls: ${parts[0]} has no column '${parts[1]}'. It has: ${header.join(', ')}")
        }
        def ids = lines.drop(1).collect { line -> line.split('\t') as List }
            .findAll { fields ->
                def value = column < fields.size() ? fields[column].trim() : ''
                wanted != null ? value == wanted : truthy.contains(value.toLowerCase())
            }
            .collect { fields -> fields[0].trim() }
        if (!ids) {
            error("--negative_controls: no row of ${parts[0]} has ${parts[1]}" + (wanted != null ? " = ${wanted}" : " set") + ".")
        }
        return ids.join(',')
    }

    if (file(text).exists()) {
        def ids = file(text).readLines().collect { line -> line.trim() }.findAll { line -> line && !line.startsWith('#') }
        if (!ids) {
            error("--negative_controls: ${text} is empty.")
        }
        return ids.join(',')
    }

    def ids = text.split(',').collect { entry -> entry.trim() }.findAll { entry -> entry }
    if (!ids) {
        error("--negative_controls: '${text}' named no samples, and is neither a file nor 'metadata.tsv:column'.")
    }
    return ids.join(',')
}

//
// Barcode geometry per chemistry. STARsolo needs the barcode length, the UMI
// start and the UMI length explicitly; getting them wrong produces a run that
// completes with almost every barcode unmatched, which is easy to miss.
//
// The whitelists themselves ship with Cell Ranger rather than with STAR, so
// --sc_whitelist has to be supplied. Passing 'None' disables correction: every
// observed barcode is then taken at face value, and a sequencing error in a
// barcode manufactures a new cell.
//
def scChemistry() {
    // Inlined rather than held in a script-level constant: Nextflow allows no
    // statements outside a process or workflow, so a top-level `def` map does
    // not parse.
    def presets = [
        '10xv2': [cb_len: 16, umi_start: 17, umi_len: 10, whitelist: '737K-august-2016.txt'],
        '10xv3': [cb_len: 16, umi_start: 17, umi_len: 12, whitelist: '3M-february-2018.txt'],
        '10xv4': [cb_len: 16, umi_start: 17, umi_len: 12, whitelist: '3M-3pgex-may-2023.txt'],
    ]
    def preset = presets[params.sc_chemistry] ?: [cb_len: null, umi_start: null, umi_len: null, whitelist: null]
    return [
        name: params.sc_chemistry,
        cb_len: params.sc_cb_len ?: preset.cb_len,
        umi_start: params.sc_umi_start ?: preset.umi_start,
        umi_len: params.sc_umi_len ?: preset.umi_len,
        expected_whitelist: preset.whitelist,
    ]
}

def commaList(value) {
    return value
        ? value.toString().split(',').collect { entry -> entry.trim() }.findAll { entry -> entry }
        : []
}

//
// Work out what a single --host entry is, so one parameter can take an NCBI
// accession, a taxid, a genome FASTA or a prebuilt HISAT2 index.
//
def classifyHostReference(entry) {
    def value = entry.toString().trim()
    def tarball = value.endsWith('.tar.gz') || value.endsWith('.tgz')

    // Accessions and taxids are identifiers, never paths.
    if (value ==~ /(?i)^GC[AF]_[0-9]+(\.[0-9]+)?$/) {
        return [kind: 'accession', value: value]
    }
    if (value ==~ /^[0-9]+$/) {
        return [kind: 'taxid', value: value]
    }

    // A remote URL cannot be probed without fetching it, so it is classified by
    // name alone; everything else is on disk and can be inspected.
    if (value.contains('://')) {
        return [kind: tarball ? 'index' : 'fasta', value: value]
    }

    def path = file(value)
    if (!path.exists()) {
        error("--host entry '${value}' is not an NCBI accession (GCF_*/GCA_*), not a taxid, and not an existing path.")
    }
    if (path.isDirectory()) {
        if (!path.list().any { name -> name.endsWith('.ht2') || name.endsWith('.ht2l') }) {
            error("--host entry '${value}' is a directory but holds no HISAT2 index (*.ht2). Point it at the index directory, or pass the genome FASTA instead.")
        }
        return [kind: 'index', value: value]
    }
    return [kind: tarball ? 'index' : 'fasta', value: value]
}

//
// Resolve --da_method into the set of differential-abundance methods to run.
// Each has its own container, so 'ancombc,aldex2' runs two processes.
//
def daMethods(da_method) {
    def known = ['aldex2', 'ancombc']
    def selection = commaList((da_method ?: 'ancombc').toString().toLowerCase())
        .collect { entry -> entry == 'ancombc2' ? 'ancombc' : entry }
    def unknown = selection.findAll { entry -> !known.contains(entry) }
    if (unknown) {
        error("--da_method: unknown method(s) ${unknown.join(', ')}. Expected any of ${known.join(', ')}.")
    }
    return selection.unique()
}

//
// Resolve --host_de_method into the set of gene-level DE methods to run. Same
// shape as daMethods(): each has its own container, so both means two processes.
//
def hostDeMethods(host_de_method) {
    def known = ['deseq2', 'edger']
    def selection = commaList((host_de_method ?: '').toString().toLowerCase())
    def unknown = selection.findAll { entry -> !known.contains(entry) }
    if (unknown) {
        error("--host_de_method: unknown method(s) ${unknown.join(', ')}. Expected any of ${known.join(', ')}.")
    }
    return selection.unique()
}

//
// The target genome, classified the same way a host reference is: accession,
// taxid, FASTA or prebuilt index, detected from the value.
//
def targetReference() {
    return params.target_reference ? classifyHostReference(params.target_reference) : null
}

def hostReferences() {
    def refs = commaList(params.host).collect { value -> classifyHostReference(value) }

    // The superseded parameters still work, but they are grouped by kind rather
    // than kept in the order they were written, so they cannot express "this
    // accession first, that local file second". --host can, because one list
    // preserves its own order - and the first entry is the primary reference.
    commaList(params.hisat2_index).each { value -> refs << [kind: 'index', value: value] }
    commaList(params.fasta).each { value -> refs << [kind: 'fasta', value: value] }
    commaList(params.host_accession).each { value -> refs << [kind: 'accession', value: value] }
    commaList(params.host_taxid).each { value -> refs << [kind: 'taxid', value: value] }

    // No truncation here. Silently dropping a host reference produces a
    // non-host fraction that is still full of host, and nothing downstream can
    // tell that from a real microbial signal - so the count is validated in
    // validateInputParameters() and refused there instead.
    return refs
}

//
// How many depletion passes the workflow declares. See the comment on the
// check in validateInputParameters() for what raising it involves.
//
def maxHostPasses() {
    return 4
}

//
// Mean read length after trimming, from fastp's JSON. Used to pick Bracken's
// -r, which has to match the read length the database was built for.
//
def meanReadLength(json_file) {
    try {
        def parsed = new groovy.json.JsonSlurper().parse(json_file.toFile())
        def after = parsed?.summary?.after_filtering
        def lengths = [after?.read1_mean_length, after?.read2_mean_length].findAll { value -> value }
        return lengths ? (lengths.sum() / lengths.size()) as Integer : null
    }
    catch (Exception e) {
        log.warn("Could not read the trimmed read length from ${json_file}: ${e.message}")
        return null
    }
}

//
// Bracken's -r must match one of the k-mer distributions in the database, not
// merely the reads: `-r 100` against a database that only has 50/150 silently
// gives you the wrong estimates. So take the observed length and snap it to the
// nearest distribution the database actually ships.
//
def brackenDistributions(db) {
    if (!db) {
        return []
    }
    def dir = file(db)
    if (!dir.exists() || !dir.isDirectory()) {
        return []                                  // a tarball; nothing to inspect yet
    }
    return dir
        .list()
        .collect { name -> (name =~ /^database(\d+)mers\.kmer_distrib$/) }
        .findAll { matcher -> matcher.matches() }
        .collect { matcher -> matcher.group(1) as Integer }
        .sort()
}

def resolveBrackenReadLength(observed, available, configured) {
    if (configured?.toString()?.toLowerCase() != 'auto') {
        return configured as Integer
    }
    if (!observed) {
        log.warn("--bracken_read_length auto could not determine the trimmed read length (no fastp JSON); falling back to 100.")
        return 100
    }
    if (!available) {
        return observed
    }
    def nearest = available.min { candidate -> Math.abs(candidate - observed) }
    if (Math.abs(nearest - observed) > 25) {
        log.warn("Trimmed reads are ~${observed} bp but the Bracken database only provides distributions for ${available.join(', ')} bp; using ${nearest}. Abundance estimates will be approximate.")
    }
    return nearest
}

//
// Resolve --ai_insights into the set of AI annotations to actually run.
//
// Accepts 'all'/'yes' (everything), 'no'/'none' (nothing) or a comma-separated
// subset. Everything is off unless --llm_endpoint is set, so the AI features are
// opt-in and the default parameter set never contacts an external service.
//
def aiInsightOptions(ai_insights, llm_endpoint) {
    def known = ['multiqc', 'qualimap', 'taxonomy'] as Set

    if (!llm_endpoint) {
        return [] as Set
    }

    def selection = (ai_insights ?: 'all').toString().toLowerCase().replaceAll(/\s/, '')

    if (['no', 'none', 'false', ''].contains(selection)) {
        return [] as Set
    }
    if (['all', 'yes', 'true'].contains(selection)) {
        return known
    }

    def chosen = selection.tokenize(',') as Set
    def unknown = chosen - known
    if (unknown) {
        error("Unknown --ai_insights value(s): ${unknown.sort().join(', ')}. Choose from ${known.sort().join(', ')}, or use 'all' / 'no'.")
    }
    return chosen
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
            params.skip_host_removal || !params.save_host_bam || params.skip_qualimap ? "" : "Qualimap (Okonechnikov et al. 2016),",
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
            params.skip_host_removal || !params.save_host_bam || params.skip_qualimap ? "" : "<li>Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2), 292-294. doi: 10.1093/bioinformatics/btv566</li>",
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


//
// Stop the resident `k2 classify --use-daemon` process at the end of a run.
//
// Two things this is careful about. It kills only a PID that is currently a
// `classify` process AND owned by whoever owns the pid file, because /tmp is
// shared node-wide and a stale pid file can name a PID the kernel has since
// handed to something else - that exact confusion is what wedged earlier runs,
// where /tmp/classify.pid named a live kernel thread. And it never throws:
// failing to tidy up must not turn a successful run into a failed one.
//
def stopKraken2Daemon() {
    if (!params.kraken2_use_daemon) {
        return
    }
    try {
        def pidFile = new File('/tmp/classify.pid')
        if (!pidFile.exists()) {
            return
        }
        def pid = pidFile.text.replaceAll(/\D/, '')
        if (!pid) {
            return
        }
        def comm = new File("/proc/${pid}/comm")
        if (!comm.exists() || comm.text.trim() != 'classify') {
            log.debug("Kraken2 daemon: /tmp/classify.pid names PID ${pid}, which is not a running classify process. Leaving it alone.")
            return
        }
        def us = fileOwner('/tmp/classify.pid')
        if (us == null || fileOwner("/proc/${pid}") != us) {
            log.warn("Kraken2 daemon: PID ${pid} is not owned by ${us}, so it belongs to another user's run. Not stopping it.")
            return
        }
        log.info("Stopping the Kraken2 daemon (PID ${pid}) so it does not hold the index after this run.")
        ["kill", "-TERM", pid].execute().waitFor()
        // A flat wait rather than a poll loop: Nextflow's strict syntax has
        // removed `while`, and the increment operator crashes 25.10.4's parser
        // outright with an internal index error rather than a syntax message.
        // Five seconds is ample for a process whose only job on SIGTERM is to
        // close two FIFOs and exit.
        Thread.sleep(5000)
        if (new File("/proc/${pid}").exists()) {
            log.warn("Kraken2 daemon (PID ${pid}) ignored SIGTERM; sending SIGKILL.")
            ["kill", "-KILL", pid].execute().waitFor()
            Thread.sleep(2000)
        }
        // File.delete() returns false for a path that is already gone, so
        // nothing here needs guarding.
        ['/tmp/classify.pid', '/tmp/classify_stdin', '/tmp/classify_stdout'].each { f ->
            new File(f).delete()
        }
    }
    catch (Exception e) {
        log.warn("Could not stop the Kraken2 daemon: ${e.message}. Stop it by hand with `k2 clean --stop-daemon` on the execution node.")
    }
}

//
// Owner of a path, or null if it cannot be read. Used to check that a daemon
// belongs to us before signalling it.
//
def fileOwner(String path) {
    try {
        return java.nio.file.Files.getOwner(java.nio.file.Paths.get(path)).getName()
    }
    catch (Exception ignored) {
        return null
    }
}

// PRISM (sjdlabgroup/PRISM), run as its authors wrote it.
//
// PRISM is not a classifier. It takes a library, classifies it with Kraken2,
// pulls the candidate microbial reads out, and then tries to CONFIRM each
// candidate taxon with independent evidence: full-length BLAST alignment
// against nt, the GenBank annotation the alignments land on, host mapping with
// STAR and minimap2, and the k-mer composition of the reads across every
// taxonomic rank. Forty features per taxon go into an XGBoost model
// (`prismxg.RDS`) and come out as one number - the probability that the taxon
// is genuinely present rather than a contaminant.
//
// TRIMMED READS, NOT NON-HOST ONES. Every other classifier in this pipeline
// takes the fraction that survived host depletion; PRISM must not. Its host
// mapping is where several of its forty features come from, so handing it
// reads that have already had the host removed empties those columns and
// changes what the model is scoring. It does its own depletion.
//
// The reads are decompressed first. PRISM builds its input paths as
// `<data_path>/<sample><fq1_end>`, and its Kraken2 call does not pass
// --gzip-compressed; rather than rely on autodetection, the ambiguity is
// removed at the cost of scratch space.
//
// NOT PACKAGED. --prism_path must point at a prepared clone holding the
// repository's own files plus two things that are not in it: `genbank/`
// (the processed GenBank RDS files, distributed separately by the authors) and
// `sorted_accession_map.txt` (built once with blastdbcmd). See docs/usage.md.
process PRISM_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${params.prism_container ?: ''}"

    input:
    tuple val(meta), path(reads)
    path prism_dir
    val kraken_db
    val blast_db
    val star_genome_dir
    val minimap2_index
    val min_qcovs
    val min_read_per
    val min_uniq_frac
    val max_sample
    val barcode_only

    output:
    tuple val(meta), path("prism_out/*-counts.csv"), emit: counts
    tuple val(meta), path("prism_out/*-results.csv"), emit: results
    tuple val(meta), path("prism_out/*.fa"), emit: fasta, optional: true
    tuple val(meta), path("prism_out/data/*-xgmat.csv"), emit: features, optional: true
    tuple val("${task.process}"), val('prism'), eval("git -C ${prism_dir} rev-parse --short HEAD 2>/dev/null || echo unversioned"), emit: versions_prism, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // With --prism_barcode_only the barcode read is DISCARDED here rather than
    // passed to PRISM with its own --barcode_only flag. PRISM's flag reaches
    // only the BLAST step - prism_blast() skips read 1 when it is set - while
    // its Kraken2 pass still runs `--paired` over both mates, so a 28 bp
    // barcode+UMI read with no biological sequence in it is still concatenated
    // to the cDNA read and classified. Dropping it here is strictly cleaner:
    // every step downstream, Kraken2 included, sees only the cDNA read.
    //
    // To get PRISM's literal behaviour instead, leave this false and pass
    // `--prism_args '--barcode_only TRUE'`; the module never sets that flag
    // itself, so there is nothing to collide with.
    def cdna_only = barcode_only && !meta.single_end
    def paired = (meta.single_end || cdna_only) ? 'FALSE' : 'TRUE'
    def decompress = meta.single_end
        ? "gzip -cdf ${reads} > prism_in/${prefix}_1.fastq"
        : cdna_only
            ? "gzip -cdf ${reads[1]} > prism_in/${prefix}_1.fastq"
            : "gzip -cdf ${reads[0]} > prism_in/${prefix}_1.fastq\n    gzip -cdf ${reads[1]} > prism_in/${prefix}_2.fastq"
    """
    mkdir -p prism_in prism_out
    ${decompress}

    # PRISM takes PATHS, not names on \$PATH - and they are not all the same
    # kind of path. kraken2, seqkit, STAR and minimap2 are validated with
    # file.exists() and must be the EXECUTABLE; blast is validated with
    # dir.exists() and used as `paste0(blast_path, "/blastn")`, so it must be
    # the DIRECTORY holding blastn, blastdbcmd and makeblastdb. PRISM's README
    # describes --blast_path as "path to blastn", which is not what the code
    # accepts.
    Rscript ${prism_dir}/PRISM.R \\
        ${args} \\
        --sample ${prefix} \\
        --data_path prism_in \\
        --out_path prism_out \\
        --prism_path ${prism_dir} \\
        --kraken_path \$(command -v kraken2) \\
        --kraken_db_path ${kraken_db} \\
        --seqkit_path \$(command -v seqkit) \\
        --star_path \$(command -v STAR) \\
        --star_genome_dir ${star_genome_dir} \\
        --blast_path \$(dirname \$(command -v blastn)) \\
        --blast_db_path ${blast_db} \\
        --minimap2_path \$(command -v minimap2) \\
        --minimap2_index ${minimap2_index} \\
        --model_org_taxids ${prism_dir}/model_org_taxids.txt \\
        --paired ${paired} \\
        --fq1_end _1.fastq \\
        --fq2_end _2.fastq \\
        --min_qcovs ${min_qcovs} \\
        --min_read_per ${min_read_per} \\
        --min_uniq_frac ${min_uniq_frac} \\
        --max_sample ${max_sample} \\
        --threads ${task.cpus}

    # The decompressed copies are the size of the library and have no use once
    # PRISM has read them.
    rm -rf prism_in
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p prism_out/data
    touch prism_out/${prefix}-counts.csv prism_out/${prefix}-results.csv \\
        prism_out/data/${prefix}-xgmat.csv
    """
}

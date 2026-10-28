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

    // The four databases arrive as values, not staged paths - they run to
    // hundreds of gigabytes - so nothing mounts them into a container and
    // PRISM would find none of them under --prism_container. Bound read-only:
    // the Kraken2 and STAR directories themselves, and the directory holding
    // the BLAST database prefix and the minimap2 index file.
    containerOptions {
        def dirs = [
            kraken_db,
            star_genome_dir,
            blast_db ? file(blast_db).parent : null,
            minimap2_index ? file(minimap2_index).parent : null,
        ].findAll { dir -> dir }.collect { dir -> dir.toString() }.unique()
        if (workflow.containerEngine in ['singularity', 'apptainer']) {
            dirs.collect { dir -> "-B ${dir}:${dir}:ro" }.join(' ')
        }
        else if (workflow.containerEngine == 'docker') {
            dirs.collect { dir -> "-v ${dir}:${dir}:ro" }.join(' ')
        }
        else {
            ''
        }
    }

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
    val bin_dir
    val decoy_taxids

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
    // Headers keep only the read ID. PRISM decides whether its read IDs carry
    // /1 and /2 by searching the WHOLE header line, so an ENA header with the
    // mate number in its comment (`@ERR2935854.1 HWI-...:2249/1`) makes it
    // append suffixes the IDs do not have; seqkit then matches no read, and
    // PRISM reports no microbial species at all - a silent empty result. A
    // header like `@ID/1` keeps its suffix, which PRISM handles correctly.
    def ids_only = "awk 'NR % 4 == 1 { sub(/[ \\t].*/, \"\") } { print }'"
    def decompress = meta.single_end
        ? "gzip -cdf ${reads} | ${ids_only} > prism_in/${prefix}_1.fastq"
        : cdna_only
            ? "gzip -cdf ${reads[1]} | ${ids_only} > prism_in/${prefix}_1.fastq"
            : "gzip -cdf ${reads[0]} | ${ids_only} > prism_in/${prefix}_1.fastq\n    gzip -cdf ${reads[1]} | ${ids_only} > prism_in/${prefix}_2.fastq"
    def dirs = bin_dir.tokenize(',').collect { dir -> dir.trim() }.findAll { dir -> dir }
    def prepend = dirs ? "export PATH=${dirs.collect { dir -> dir + '/bin' }.join(':')}:\$PATH" : ''
    // PRISM's model organisms are decoys: a read whose best BLAST hit is one
    // of them is not counted for the microbe Kraken2 called. Its bundled list
    // (mouse, rat, fly, worm, zebrafish, Arabidopsis) suits the human data its
    // model was trained on, where the host itself is the decoy; for any other
    // host the decoy should be that host. PRISM adds human (9606) itself.
    // BLAST expands taxids to their descendants only when taxonomy4blast.sqlite3
    // sits beside the database, so they are expanded here from the Kraken2
    // database's nodes.dmp, and passed as they are where it has none.
    def seeds = decoy_taxids.toString().tokenize(',').collect { taxid -> taxid.trim() }.findAll { taxid -> taxid }
    def decoys = seeds == ['prism']
        ? "cp ${prism_dir}/model_org_taxids.txt prism_decoys.txt"
        : """
    if [ -s ${kraken_db}/nodes.dmp ]; then
        awk -F '\\t[|]\\t' -v seeds="${seeds.join(',')}" '
            { kids[\$2] = kids[\$2] " " \$1 }
            END {
                n = split(seeds, queue, ",")
                for (i = 1; i <= n; i++) seen[queue[i]] = 1
                for (i = 1; i <= n; i++) {
                    m = split(kids[queue[i]], c, " ")
                    for (j = 1; j <= m; j++) if (!(c[j] in seen) && c[j] != queue[i]) { seen[c[j]] = 1; queue[++n] = c[j] }
                }
                for (t in seen) print t
            }' ${kraken_db}/nodes.dmp | sort -n > prism_decoys.txt
    else
        printf '%s\\n' ${seeds.join(' ')} > prism_decoys.txt
    fi"""
    """
    ${prepend}
    mkdir -p prism_in prism_out prism_bin
    ${decompress}
    ${decoys}

    # PRISM's own Kraken2 pass loads the whole database into its process, so
    # every PRISM task held a private copy (~330 GB for core_nt) and a few at
    # once exhausted a 1 TB node. With --memory-mapping the tasks share one
    # copy in the page cache; reading it through once first spares the mapped
    # pass from faulting it in page by page over a network filesystem.
    printf '#!/usr/bin/env bash\\nexec %s --memory-mapping "\$@"\\n' "\$(command -v kraken2)" > prism_bin/kraken2
    chmod +x prism_bin/kraken2
    cat ${kraken_db}/*.k2d > /dev/null

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
        --data_path \$PWD/prism_in \\
        --out_path \$PWD/prism_out \\
        --prism_path \$PWD/${prism_dir} \\
        --kraken_path \$PWD/prism_bin/kraken2 \\
        --kraken_db_path ${kraken_db} \\
        --seqkit_path \$(command -v seqkit) \\
        --star_path \$(command -v STAR) \\
        --star_genome_dir ${star_genome_dir} \\
        --blast_path \$(dirname \$(command -v blastn)) \\
        --blast_db_path ${blast_db} \\
        --minimap2_path \$(command -v minimap2) \\
        --minimap2_index ${minimap2_index} \\
        --model_org_taxids \$PWD/prism_decoys.txt \\
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

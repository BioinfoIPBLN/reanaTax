// Metax: species calls from reference alignment plus genome-wide coverage
// evidence (Deng, Safaei & McHardy 2026, Cell), as an alternative or a
// companion to Kraken2/Bracken.
//
// Reads are aligned to a taxonomically annotated genome collection with the MA
// aligner (maCMD) and multi-mapping reads are shared out by EM. Every candidate
// taxon is then asked whether its reads cover the genome the way random
// sampling at that depth would: observed over expected breadth (OEBR), the same
// at chunk level (COEBR), and a probability on each. A taxon whose reads pile
// onto one region - a contaminated contig in the reference, a kitome fragment,
// a conserved gene it shares with the organism really present - fails and is
// dropped. That is the check this pipeline otherwise approximates with host
// k-mers and minimizer counts, made directly on coverage. It reports NCBI
// taxids, so its output lines up with Kraken2's without a name mapping.
//
// Its authors are explicit about the price: true taxa can be over-filtered
// where genome-wide coverage is not expected, "such as in RNA-derived data". A
// transcriptome covers what is expressed, so on RNA a missing Metax call is
// weak evidence of absence. RNA viruses are the exception worth knowing: there
// the transcript is the genome.
//
// There is no public container. Metax is published on its author's conda
// channel, not on bioconda, and that package ships the metax binary without
// the maCMD aligner it calls, so the environment below supplies both.
// conf/modules.config sets an image only when --metax_container names one;
// leaving it empty is what lets `-profile wave` build one from this
// environment.
//
// The SAM that maCMD writes is removed once the profile exists. It holds every
// read's alignment - under --skip_host_removal, every host read, against the
// human reference where the database includes one - and nothing downstream
// reads it.
// The per-read classify.txt is kept, gzipped, only when asked for.
process METAX_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(reads)
    path db, stageAs: 'metax_db'
    val db_json
    path taxonomy, stageAs: 'taxonomy'
    val save_classify

    output:
    tuple val(meta), path("*.profile.txt")    , emit: profile
    tuple val(meta), path("*.classify.txt.gz"), emit: classify, optional: true
    tuple val(meta), path("*.log")            , emit: log
    tuple val("${task.process}"), val('metax'), eval("metax --version 2>&1 | sed 's/^metax //'"), emit: versions_metax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files = [reads].flatten()
    def input = meta.single_end ? "-i ${files.join(',')}" : "-i ${files[0]},${files[1]} -p"
    def classify = save_classify
        ? "if [ -f ${prefix}.classify.txt ]; then gzip -n ${prefix}.classify.txt; fi"
        : "rm -f ${prefix}.classify.txt"
    """
    metax profile \\
        --db metax_db/${db_json} \\
        --dmp-dir taxonomy \\
        ${input} \\
        -o ${prefix} \\
        -t ${task.cpus} \\
        ${args}

    rm -f ${prefix}.sam ${prefix}.sam.gz
    ${classify}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.profile.txt ${prefix}.log
    """
}

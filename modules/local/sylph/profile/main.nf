// sylph: species-level profiling by containment ANI (Shaw & Yu 2024, Nat
// Biotechnol 43:1348), as an alternative or a companion to Kraken2/Bracken.
//
// It measures something different from Kraken2. Kraken2 assigns every read by
// its k-mers; sylph never assigns a read at all. It sketches the sample (about
// one k-mer in c, c = 200 by default), asks how much of each reference
// genome's sketch the sample contains, corrects that containment for low
// coverage with a zero-inflated Poisson model, and reports a genome only when
// the resulting ANI clears 95%. A call is therefore a statement about the whole
// genome, which is why it makes far fewer false positives than per-read
// classification - and why it needs reads that sample the genome rather than a
// handful of loci. Its documented floor is 0.01-0.05x coverage, and a
// transcriptome samples only what is expressed: uneven coverage lowers the
// containment, the ANI with it, and a present organism can land under 95% and
// go unreported. On RNA, a missing sylph call is weak evidence of absence.
//
// Reads are sketched inside `sylph profile` rather than in a separate step:
// sylph takes FASTQ directly, sketching is the cheap part, and one task per
// sample keeps every library resumable on its own.
//
// Sample_file is rewritten to the sample id. sylph records the read file name
// there - for a pair, the R1 path - and sylph-tax names its output files and
// its merged columns after that field, so leaving it would put
// `SRX123_1.fastq.gz` where every other table in the run says `SRX123`.
process SYLPH_PROFILE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph:1.0.0--hb42e459_0'
        : 'quay.io/biocontainers/sylph:1.0.0--hb42e459_0'}"

    input:
    tuple val(meta), path(reads)
    path databases, stageAs: 'db/*'

    output:
    tuple val(meta), path("*.sylph.tsv"), emit: profile
    tuple val("${task.process}"), val('sylph'), eval("sylph -V | sed 's/sylph //'"), emit: versions_sylph, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    sylph profile \\
        -t ${task.cpus} \\
        ${args} \\
        db/* \\
        ${input} \\
        -o sylph_raw.tsv

    awk -v id="${prefix}" 'BEGIN { FS = OFS = "\\t" } NR == 1 { print; next } { \$1 = id; print }' \\
        sylph_raw.tsv > ${prefix}.sylph.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sylph.tsv
    """
}

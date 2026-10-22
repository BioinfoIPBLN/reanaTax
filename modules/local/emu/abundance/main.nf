// Emu: species-level profiling of full-length 16S reads by expectation
// maximisation over minimap2 alignments (Curry et al., Nat Methods 2022).
//
// Emu aligns with `minimap2 -ax <type> -N 50 -p .9` straight against the
// database's species_taxid.fasta, which minimap2 indexes in memory in seconds,
// and writes only to its output directory; unlike MetaPhlAn it needs no shared
// index built beforehand. The minimap2 preset follows the platform of the
// library (map-ont for ONT, map-hifi for PacBio) unless --emu_args names one
// with --type, e.g. `--type lr:hq` for ONT reads of Q20 and above or
// `--type map-pb` for PacBio CLR.
process EMU_ABUNDANCE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/emu:3.6.2--hdfd78af_0'
        : 'quay.io/biocontainers/emu:3.6.2--hdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    path db

    output:
    tuple val(meta), path("*_rel-abundance.tsv"), emit: profile
    tuple val(meta), path("*_rel-abundance-threshold-*.tsv"), optional: true, emit: profile_threshold
    tuple val(meta), path("*_{unmapped,unclassified_mapped,filtered_mapped}*.fast{a,q}"), optional: true, emit: unclassified
    tuple val(meta), path("*.emu.log"), emit: log
    tuple val("${task.process}"), val('emu'), eval("emu --version 2>&1 | sed 's/^emu v//'"), emit: versions_emu, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def given = args.tokenize(' ')
    def type = given.contains('--type') || given.contains('-x') ? '' : "--type ${meta.long_reads == 'pacbio' ? 'map-hifi' : 'map-ont'}"
    """
    emu abundance \\
        ${reads} \\
        --db ${db} \\
        ${type} \\
        --threads ${task.cpus} \\
        --output-dir . \\
        --output-basename ${prefix} \\
        ${args} \\
        2>&1 | tee ${prefix}.emu.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rel-abundance.tsv ${prefix}.emu.log
    """
}

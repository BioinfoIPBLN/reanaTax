// Taxpasta: one profiler's per-sample outputs as a single standardised table.
//
// Every profiler writes its own format - a kreport, Bracken's abundance table,
// KrakenUniq's report, a MetaPhlAn profile - so comparing two of them on the
// same library means writing a parser for each. Taxpasta reads all of those
// and writes one shape: a taxonomy_id column and a count column per sample.
//
// Sample names come from a samplesheet rather than from file names, which
// Taxpasta would otherwise use whole (`SRX123.kraken2.report` as a column).
// `taxpasta merge` needs at least two profiles, so a single-sample run is
// written with `taxpasta standardise` instead, in Taxpasta's long layout.
//
// MetaPhlAn is trimmed to its first three columns first. Taxpasta insists on
// MetaPhlAn's four-column default output, and this pipeline runs MetaPhlAn with
// `-t rel_ab_w_read_stats` by default, which writes five. The three kept -
// clade, taxid, relative abundance - are the only ones Taxpasta reads. The
// trimmed copies are written beside the staged profiles, never over them: those
// are links into other tasks' directories.
process TAXPASTA_MERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/taxpasta:0.7.0--pyhdfd78af_1'
        : 'biocontainers/taxpasta:0.7.0--pyhdfd78af_1'}"

    input:
    tuple val(meta), val(samples), path(profiles, stageAs: 'profiles/*')
    val profiler
    path taxonomy

    output:
    tuple val(meta), path("*.taxpasta.*"), emit: table
    tuple val("${task.process}"), val('taxpasta'), eval("taxpasta --version | sed 's/.* //'"), emit: versions_taxpasta, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'tsv'
    def files = [profiles].flatten()
    def source = profiler == 'metaphlan' ? 'trimmed' : 'profiles'
    def pairs = [samples, files].transpose()
    def sheet_args = pairs.collect { sample, profile -> "'${sample}' '${source}/${profile.name}'" }.join(' ')
    def tax = taxonomy ? "--taxonomy ${taxonomy}" : ''
    def trim = profiler == 'metaphlan'
        ? """mkdir -p trimmed
    for f in profiles/*; do
        awk 'BEGIN { FS = OFS = "\\t" } /^#/ { print; next } { print \$1, \$2, \$3, "" }' "\$f" > "trimmed/\$(basename "\$f")"
    done"""
        : ''
    """
    ${trim}
    printf 'sample\\tprofile\\n' > samplesheet.tsv
    printf '%s\\t%s\\n' ${sheet_args} >> samplesheet.tsv

    if [ ${pairs.size()} -gt 1 ]; then
        taxpasta merge \\
            --profiler ${profiler} \\
            --samplesheet samplesheet.tsv \\
            --output ${prefix}.taxpasta.${format} \\
            ${tax} \\
            ${args}
    else
        taxpasta standardise \\
            --profiler ${profiler} \\
            --output ${prefix}.taxpasta.${format} \\
            ${tax} \\
            ${args} \\
            ${source}/${files[0].name}
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def format = task.ext.format ?: 'tsv'
    """
    touch ${prefix}.taxpasta.${format}
    """
}

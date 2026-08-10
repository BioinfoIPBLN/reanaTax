// Convert a Kraken2 report into a MetaPhlAn-style lineage profile.
//
// This exists purely to feed HUMAnN, which selects the pangenomes to align
// against from a taxonomic profile and only speaks MetaPhlAn's format. Rather
// than run MetaPhlAn as a second classifier - a second database, a second set
// of abundances that would disagree with the Bracken tables the report is built
// from - the Kraken2 report we already have is translated into that format.
//
// The translation is lossy in one direction that matters: HUMAnN keys on
// `s__` species lines, so anything Kraken2 could only place above species is
// invisible to the functional step. See docs/usage.md.
process KRAKENTOOLS_KREPORT2MPA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakentools:1.2.1--pyh7e72e81_0'
        : 'biocontainers/krakentools:1.2.1--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.mpa.txt"), emit: mpa
    tuple val("${task.process}"), val('krakentools'), eval("echo 1.2.1"), emit: versions_krakentools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    kreport2mpa.py \\
        ${args} \\
        --report ${report} \\
        --output ${prefix}.mpa.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.mpa.txt
    """
}

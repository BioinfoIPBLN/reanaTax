//
// Everything that is not a microbe, listed for removal from the final tables.
//
// Host reads that escape depletion and miss the host's own clade still land
// somewhere: on another vertebrate, an insect or a grass whose assembly holds
// a conserved gene or contaminating host sequence. On the mouse stroke
// single-cell cohort those calls outnumbered the microbes nine to one per
// library, and --drop_host_clade phylum cannot reach the invertebrates and
// plants among them.
//
// Like HOSTCLADE_EXPAND it removes nothing itself: the taxids join the shared
// drop list, so ABUNDANCE_FILTER owns every removal and the single-cell matrix
// loses the same rows. The anchors (Metazoa, Viridiplantae, other and
// unclassified sequences) are read from the run's own reports. Aliased with
// `--environmental` (task.ext.kind 'environmental') it lists uncultured,
// unidentified and environmental-sample taxa instead.
//
process NONMICROBIAL_EXPAND {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')

    output:
    tuple val(meta), path("*_drop.txt")    , emit: drop_list
    tuple val(meta), path("*_evidence.tsv"), emit: evidence
    tuple val(meta), path("*_mqc.tsv")                   , emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kind = task.ext.kind ?: 'non_microbial'
    """
    non_microbial.py \\
        --reports reports/* \\
        --drop-list ${prefix}.${kind}_drop.txt \\
        --evidence ${prefix}.${kind}_evidence.tsv \\
        --mqc ${prefix}.${kind}_mqc.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def kind = task.ext.kind ?: 'non_microbial'
    """
    touch ${prefix}.${kind}_drop.txt ${prefix}.${kind}_evidence.tsv
    """
}

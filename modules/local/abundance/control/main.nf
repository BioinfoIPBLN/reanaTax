// Judge a taxon against the level that arrives WITHOUT a sample.
//
// The only filter here whose verdict is per LIBRARY rather than per taxon, and
// that is the whole point of it. Index hopping, well-to-well carryover and
// ambient template put genuine reads of a genuine organism into the wrong
// library; nothing about those reads is wrong, so no evidence filter can see
// them, and no single verdict per taxon can express "real here, carryover
// there". A control library measures that background directly.
//
// Unlike every other filter in this pipeline it therefore REWRITES the table
// rather than only naming taxids: cells that fail are zeroed in place, and it
// runs before ABUNDANCE_FILTER so the surviving abundances are recomputed
// against what is left. Taxa that fail in every real library, and taxa the
// prevalence rule condemns outright, still leave by the usual route - a taxid
// list ABUNDANCE_FILTER removes - so those removals stay in one place.
//
// Complementary to DECONTAM_FILTER, not a replacement: decontam is given a
// measurement of the KIT and asks whether a taxon is reagent; this is given
// control LIBRARIES and asks whether there is more of a taxon here than turns
// up on its own. See bin/control_filter.py for why decontam's prevalence
// method cannot settle the plate case, with the numbers.
process ABUNDANCE_CONTROL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(table, stageAs: 'input/*')
    val controls
    val ratio
    val statistic
    val min_reads
    val floor_reads
    val prevalence
    val prevalence_min_reads

    output:
    tuple val(meta), path("*.control_filtered.tsv"), emit: filtered
    tuple val(meta), path("*.control_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.control_drop.txt"), emit: drop_list
    tuple val(meta), path("*.control_cells.tsv"), emit: drop_cells
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // A comma-separated list of sample IDs, already resolved by the pipeline
    // from whichever form --negative_controls took. Empty when only
    // --prevalence_filter is in play, which the script accepts.
    def controls_arg = controls ? "--controls '${controls}'" : ''
    """
    control_filter.py \\
        ${args} \\
        --input ${table} \\
        ${controls_arg} \\
        --ratio ${ratio} \\
        --statistic ${statistic} \\
        --min-reads ${min_reads} \\
        --floor-reads ${floor_reads} \\
        --prevalence ${prevalence} \\
        --prevalence-min-reads ${prevalence_min_reads} \\
        --output ${prefix}.control_filtered.tsv \\
        --evidence ${prefix}.control_evidence.tsv \\
        --drop-list ${prefix}.control_drop.txt \\
        --drop-cells ${prefix}.control_cells.tsv \\
        --mqc ${prefix}_control_filter_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.control_filtered.tsv ${prefix}.control_evidence.tsv \\
          ${prefix}.control_drop.txt ${prefix}.control_cells.tsv \\
          ${prefix}_control_filter_mqc.tsv
    """
}

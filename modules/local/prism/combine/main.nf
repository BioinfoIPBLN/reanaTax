// PRISM's per-sample verdicts into one table per cohort.
//
// Reads and scores are kept in separate matrices because a score has no
// meaning where there were no reads to score: a taxon absent from a sample is
// blank in the score matrix, not zero. Zero would read as "PRISM was confident
// this is a contaminant", which is the opposite of "PRISM never saw it".
//
// The drop list is aggregated on the MEDIAN score across the samples that saw
// the taxon - not the mean, which one confidently-scored library drags either
// way, and not "any sample", which on a large cohort condemns everything.
process PRISM_COMBINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(counts, stageAs: 'counts/*')
    val score_threshold
    val min_reads

    output:
    tuple val(meta), path("*.prism_reads.tsv"), emit: reads
    tuple val(meta), path("*.prism_score.tsv"), emit: score
    tuple val(meta), path("*.prism_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.prism_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_prism.py \\
        ${args} \\
        --counts counts/* \\
        --score-threshold ${score_threshold} \\
        --min-reads ${min_reads} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prism_reads.tsv ${prefix}.prism_score.tsv \\
        ${prefix}.prism_evidence.tsv ${prefix}.prism_drop.txt ${prefix}_prism_mqc.tsv
    """
}

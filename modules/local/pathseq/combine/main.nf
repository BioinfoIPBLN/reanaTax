// PathSeq's per-sample scores into one table per quantity.
//
// Kept in PathSeq's own units rather than converted into a Kraken report: the
// two tools disagree about what "a read assigned to a taxon" means - PathSeq
// divides an ambiguous read between the taxa it aligns to, Kraken2 pushes it up
// to their common ancestor - and a conversion would have to pick one convention
// and hide the other. See bin/merge_pathseq.py.
process PATHSEQ_COMBINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(scores, stageAs: 'scores/*')
    val rank

    output:
    tuple val(meta), path("*_reads.tsv"), emit: reads
    tuple val(meta), path("*_unambiguous.tsv"), emit: unambiguous
    tuple val(meta), path("*_score.tsv"), emit: score
    tuple val(meta), path("*_score_normalized.tsv"), emit: score_normalized
    tuple val(meta), path("*_lineage.tsv"), emit: lineage
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rank_arg = rank ? "--rank ${rank}" : ''
    """
    merge_pathseq.py \\
        ${args} \\
        --scores scores/* \\
        ${rank_arg} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_reads.tsv ${prefix}_unambiguous.tsv ${prefix}_score.tsv \\
        ${prefix}_score_normalized.tsv ${prefix}_lineage.tsv ${prefix}_pathseq_mqc.tsv
    """
}

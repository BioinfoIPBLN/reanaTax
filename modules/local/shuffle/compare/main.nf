// Compare the real classification against the shuffled one.
//
// Cohort-level, like the minimizer and host-k-mer filters, and for the same
// reason: how much of a taxon's signal is composition is a property of the
// taxon, and chance matches are far too rare per sample for a per-sample ratio
// to mean anything. The taxids it condemns go to ABUNDANCE_FILTER, so one step
// still owns every removal from the combined tables.
//
// The two report sets are staged into separate directories under the SAME file
// names, because bin/shuffle_compare.py pairs them by basename - staging them
// side by side would collide.
process SHUFFLE_COMPARE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(real, stageAs: 'real/*'), path(shuffled, stageAs: 'shuffled/*')
    val max_ratio
    val min_reads

    output:
    tuple val(meta), path("*.shuffle_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.shuffle_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    shuffle_compare.py \\
        ${args} \\
        --real real/* \\
        --shuffled shuffled/* \\
        --max-ratio ${max_ratio} \\
        --min-reads ${min_reads} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.shuffle_evidence.tsv ${prefix}.shuffle_drop.txt ${prefix}_shuffle_mqc.tsv
    """
}

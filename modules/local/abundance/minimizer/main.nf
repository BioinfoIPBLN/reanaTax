// Judge taxa on the BREADTH of their evidence, not its volume.
//
// One task for the cohort: a taxon is taken at its best showing across
// libraries, so breadth demonstrated once is enough and a shallow library
// cannot condemn an organism a deep one resolved. The output is a taxid list
// that ABUNDANCE_FILTER removes from the combined tables, plus the evidence
// table the decision was made from - the same audit trail the abundance filter
// keeps.
//
// Works from either classifier's report; see bin/minimizer_filter.py for the
// two layouts and what is derived from each.
process ABUNDANCE_MINIMIZER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')
    path inspect
    val min_reads
    val min_distinct
    val max_duplication
    val min_coverage
    val distinct_scale

    output:
    tuple val(meta), path("*.minimizer_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.minimizer_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inspect_arg = inspect ? "--inspect ${inspect}" : ''
    """
    minimizer_filter.py \\
        ${args} \\
        --reports reports/* \\
        ${inspect_arg} \\
        --min-reads ${min_reads} \\
        --min-distinct ${min_distinct} \\
        --max-duplication ${max_duplication} \\
        --min-coverage ${min_coverage} \\
        --distinct-scale ${distinct_scale} \\
        --evidence ${prefix}.minimizer_evidence.tsv \\
        --drop-list ${prefix}.minimizer_drop.txt \\
        --mqc ${prefix}_minimizer_filter_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.minimizer_evidence.tsv ${prefix}.minimizer_drop.txt ${prefix}_minimizer_filter_mqc.tsv
    """
}

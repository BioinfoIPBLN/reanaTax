// Ask the LLM to describe the combined taxonomic profile in a few sentences.
//
// Two outputs on purpose: a markdown file next to the tables, and a MultiQC
// custom-content section so the summary shows up at the top of the report
// rather than in a file nobody opens. That is also why this runs BEFORE
// MULTIQC: its `*_mqc.html` is one of MultiQC's inputs.
//
// Every failure mode here is non-fatal (no endpoint configured, unreachable
// server, empty answer): the script exits 0 having written nothing, and the
// outputs are optional, so an unreachable LLM costs the report a section and
// nothing else.
process LLM_INSIGHT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(tables)
    val llm_endpoint
    val llm_model
    val llm_api_key

    output:
    tuple val(meta), path("*.ai_insight.md"), emit: insight, optional: true
    tuple val(meta), path("*_mqc.html"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    llm_insight.py \\
        ${args} \\
        --input ${tables} \\
        --out ${prefix}.ai_insight.md \\
        --mqc-out ${prefix}_ai_insight_mqc.html \\
        --endpoint '${llm_endpoint}' \\
        --model '${llm_model}' \\
        --api-key '${llm_api_key}'
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.ai_insight.md
    touch ${prefix}_ai_insight_mqc.html
    """
}

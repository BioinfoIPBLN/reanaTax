// Add per-section AI summaries to the MultiQC report, and scrub the LLM
// endpoint from everything MultiQC wrote.
//
// The report and the data directory are COPIED before being touched: a staged
// input is a symlink into the MULTIQC task's work directory, so annotating it
// in place would rewrite another task's published output and break -resume.
//
// This process is also the redaction point for MultiQC's own AI feature
// (--multiqc_ai_builtin), which bakes the endpoint URL into the HTML. That is
// why the script redacts unconditionally, even when the LLM pass is skipped.
process MULTIQC_AI {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(report, stageAs: 'input/*'), path(data_dir, stageAs: 'input/*')
    val llm_endpoint
    val llm_model
    val llm_api_key
    val annotate

    output:
    tuple val(meta), path("*.html"), emit: report
    tuple val(meta), path("*_data"), emit: data
    // Deliberately NOT pushed to the `versions` topic. MULTIQC's input channel
    // waits for that topic to close, and this process runs after MULTIQC, so
    // publishing to it would make the two wait for each other forever - the
    // same trap the nf-core MULTIQC module documents. The interpreter is the
    // one LLM_INSIGHT already reports.
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def mode = annotate ? '' : '--redact-only'
    // NOT `report.name`: for a `stageAs: 'input/*'` input that returns the whole
    // staged name, `input/multiqc_report.html`, which would point the annotator
    // back at the symlink it is supposed to be copying away from.
    def report_name = report.toString().tokenize('/').last()
    def data_name = data_dir.toString().tokenize('/').last()
    """
    cp -L ${report} ${report_name}
    cp -RL ${data_dir} ${data_name}
    chmod -R u+w ${report_name} ${data_name}

    multiqc_ai.py \\
        ${mode} \\
        ${args} \\
        --report ${report_name} \\
        --data-dir ${data_name} \\
        --llm-endpoint '${llm_endpoint}' \\
        --llm-model '${llm_model}' \\
        --llm-api-key '${llm_api_key}'
    """

    stub:
    def report_name = report.toString().tokenize('/').last()
    def data_name = data_dir.toString().tokenize('/').last()
    """
    cp -L ${report} ${report_name}
    cp -RL ${data_dir} ${data_name}
    chmod -R u+w ${report_name} ${data_name}
    """
}

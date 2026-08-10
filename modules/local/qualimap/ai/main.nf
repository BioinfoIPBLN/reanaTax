// Add an AI summary box to one sample's Qualimap BamQC report.
//
// The report directory is COPIED before being touched, for the same reason as
// in MULTIQC_AI: a staged input is a symlink into the QUALIMAP_BAMQC task's
// work directory.
//
// `maxForks 1` (set in conf/base.config) is load-bearing, not a resource hint:
// it is what stops N samples from issuing N concurrent requests to a single
// LLM server. See docs/usage.md.
process QUALIMAP_AI {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(qualimap_dir, stageAs: 'input/*')
    val llm_endpoint
    val llm_model
    val llm_api_key

    output:
    tuple val(meta), path("${prefix}"), emit: results
    // Not pushed to the `versions` topic: this runs after MULTIQC, whose input
    // waits for that topic to close. See the note in modules/local/multiqcai.
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -RL ${qualimap_dir} ${prefix}
    chmod -R u+w ${prefix}

    qualimap_ai.py \\
        ${args} \\
        --analysis-dir ${prefix} \\
        --llm-endpoint '${llm_endpoint}' \\
        --llm-model '${llm_model}' \\
        --llm-api-key '${llm_api_key}'
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp -RL ${qualimap_dir} ${prefix}
    chmod -R u+w ${prefix}
    """
}

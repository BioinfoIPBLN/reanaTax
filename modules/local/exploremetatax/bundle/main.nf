// Package the taxonomic tables for the exploreMetaTax Shiny app.
//
// One task for the whole cohort: the inputs are tables, not reads, and a bundle
// holding one sample would defeat the point of an app built for comparing them.
//
// Each category is staged into its own subdirectory because the renaming rules
// differ per category and several of reanaTax's published names are ambiguous
// on their own - `<id>.tsv` from the MetaPhlAn merge could be anything. The
// script does the renaming and refuses to write an archive whose filenames the
// app would load under the wrong format; see bin/exploremetatax_bundle.py.
process EXPLOREMETATAX_BUNDLE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    val meta
    path kraken2_reports, stageAs: 'kraken2/*'
    path krakenuniq_reports, stageAs: 'krakenuniq/*'
    path bracken_tables, stageAs: 'bracken/*'
    path bracken_combined, stageAs: 'bracken_combined/*'
    path metaphlan_profiles, stageAs: 'metaphlan/*'
    path metaphlan_merged, stageAs: 'metaphlan_merged/*'
    path humann_tables, stageAs: 'humann/*'
    path metadata, stageAs: 'metadata/*'

    output:
    tuple val(meta), path("*.exploreMetaTax.tar.gz"), emit: bundle
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_arg = metadata ? "--metadata ${metadata}" : ''
    """
    exploremetatax_bundle.py \\
        ${args} \\
        --kraken2 kraken2 \\
        --krakenuniq krakenuniq \\
        --bracken bracken \\
        --bracken-combined bracken_combined \\
        --metaphlan metaphlan \\
        --metaphlan-merged metaphlan_merged \\
        --humann humann \\
        ${metadata_arg} \\
        --outdir ${prefix}_exploreMetaTax \\
        --output ${prefix}.exploreMetaTax.tar.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.exploreMetaTax.tar.gz
    """
}

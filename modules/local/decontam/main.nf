// Identify reagent contaminants against the negative controls.
//
// The only filter in this pipeline given an external measurement of the kit,
// and therefore the only one that can separate a reagent contaminant from a
// genuinely rare organism. Every other filter judges a taxon by its own
// evidence and is blind to that distinction by construction.
//
// Entirely opt-in: most public datasets carry neither blanks nor DNA
// concentrations, and a reanalysis pipeline that required them could not run on
// the archives it exists to mine. bin/decontam_filter.R fails loudly when the
// inputs for the requested method are absent rather than reporting a clean run.
//
// Like the other evidence filters it removes nothing; the taxids go to
// ABUNDANCE_FILTER so one step owns every removal from the combined tables.
process DECONTAM_FILTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-decontam:1.26.0--r44hdfd78af_0'
        : 'quay.io/biocontainers/bioconductor-decontam:1.26.0--r44hdfd78af_0'}"

    input:
    tuple val(meta), path(counts)
    path metadata
    val neg_column
    val neg_value
    val conc_column
    val method
    val threshold
    val batch_column
    val batch_combine

    output:
    tuple val(meta), path("*.decontam_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.decontam_drop.txt"), emit: drop_list
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('decontam'), eval("Rscript -e 'cat(as.character(packageVersion(\"decontam\")))'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def neg_arg = neg_column ? "--neg-column '${neg_column}' --neg-value '${neg_value}'" : ''
    def conc_arg = conc_column ? "--conc-column '${conc_column}'" : ''
    def batch_arg = batch_column ? "--batch-column '${batch_column}' --batch-combine ${batch_combine}" : ''
    """
    decontam_filter.R \\
        ${args} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        ${neg_arg} \\
        ${conc_arg} \\
        ${batch_arg} \\
        --method ${method} \\
        --threshold ${threshold} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.decontam_evidence.tsv ${prefix}.decontam_drop.txt ${prefix}_decontam_mqc.tsv
    """
}

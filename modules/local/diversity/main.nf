// Alpha diversity, and how much of the between-sample variation each metadata
// variable explains.
//
// Runs on the FILTERED combined Bracken table, unlike differential abundance,
// and the difference is deliberate. ALDEx2 and ANCOM-BC2 apply their own
// prevalence filters, so feeding them a pre-thinned table would filter twice.
// Diversity indices apply none: observed richness on an unfiltered Kraken2 or
// Bracken table is very largely a count of database artefacts and single-read
// tails, and Shannon inherits that. Whatever the abundance and evidence filters
// removed should stay removed here.
//
// Metadata is optional. Without it there is still alpha diversity and an
// unconstrained ordination; see bin/diversity.R.
process DIVERSITY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'
        : 'quay.io/biocontainers/bioconductor-phyloseq:1.50.0--r44hdfd78af_0'}"

    input:
    tuple val(meta), path(counts)
    path metadata
    val permutations

    output:
    tuple val(meta), path("*.alpha_diversity.tsv"), emit: alpha
    tuple val(meta), path("*.beta_variance.tsv"), optional: true, emit: beta
    tuple val(meta), path("*.ordination.tsv"), optional: true, emit: ordination
    tuple val(meta), path("*.aitchison_distance.tsv"), emit: distance
    tuple val(meta), path("*.png"), optional: true, emit: plots
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('vegan'), eval("Rscript -e 'cat(as.character(packageVersion(\"vegan\")))'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def metadata_arg = metadata ? "--metadata ${metadata}" : ''
    """
    diversity.R \\
        ${args} \\
        --counts ${counts} \\
        ${metadata_arg} \\
        --permutations ${permutations} \\
        --prefix ${prefix}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alpha_diversity.tsv ${prefix}.aitchison_distance.tsv
    """
}

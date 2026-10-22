//
// A minimap2 index of a host genome, for the long-read libraries. Built with the
// preset the libraries are aligned with (meta.preset): k-mer and window size
// are fixed in the .mmi, and minimap2 ignores the indexing half of -x when it
// loads a prebuilt index, so an index made for map-ont would quietly align
// spliced reads with the wrong seeds. One index per preset the run needs.
//
process MINIMAP2_INDEX {
    tag "${meta.id}:${meta.preset}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.mmi"), emit: index
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), emit: versions_minimap2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${fasta.baseName}.${meta.preset.replaceAll(/[^A-Za-z0-9_-]/, '_')}"
    """
    minimap2 -x ${meta.preset} -t ${task.cpus} ${args} -d ${prefix}.mmi ${fasta}
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}.${meta.preset.replaceAll(/[^A-Za-z0-9_-]/, '_')}"
    """
    touch ${prefix}.mmi
    """
}

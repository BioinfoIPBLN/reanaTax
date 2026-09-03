//
// kallisto index over a transcriptome FASTA.
//
process KALLISTO_INDEX {
    tag "${fasta.name}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--ha4fb952_1'
        : 'quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1'}"

    input:
    path fasta

    output:
    path "kallisto.idx", emit: index
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*version //'"), emit: versions_kallisto, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    kallisto index -i kallisto.idx ${args} ${fasta}
    """

    stub:
    """
    touch kallisto.idx
    """
}

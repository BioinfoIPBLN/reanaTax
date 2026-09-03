process HISAT2_EXTRACTSPLICESITES {
    tag "${gtf}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hisat2:2.2.2--h503566f_0' :
        'quay.io/biocontainers/hisat2:2.2.2--h503566f_0'}"

    input:
    tuple val(meta), path(gtf)

    output:
    tuple val(meta), path("*.splice_sites.txt"), emit: txt
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n 's/.*version \\([^ ]*\\).*/\\1/p'"), emit: versions_hisat2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    hisat2_extract_splice_sites.py ${gtf} > ${gtf.baseName}.splice_sites.txt
    """

    stub:
    """
    touch ${gtf.baseName}.splice_sites.txt
    """
}

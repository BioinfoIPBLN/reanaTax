// The minimap2 index MetaPhlAn's --long_reads mode maps against, built once
// for the run.
//
// Left to itself, MetaPhlAn builds this index inside --db_dir the first time it
// runs: every sample of a cohort that starts together finds it missing and
// writes the same file at once, and a shared read-only database cannot take it
// at all. It also only adopts an index it finds already there when the file is
// over 30 GB; below that it warns and hands minimap2 the directory. So the
// index is built here, with MetaPhlAn's own command (`-v 0 -I 20G` and the
// presets of --metaphlan_minimap2_ps), and METAPHLAN_PROFILE maps against it
// itself. The name is MetaPhlAn's, so a published copy can also sit in a
// database directory for MetaPhlAn run by hand.
process METAPHLAN_MM2INDEX {
    tag "${index}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metaphlan:4.2.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/metaphlan:4.2.6--pyhdfd78af_0'}"

    input:
    path db
    val index
    val minimap2_ps

    output:
    path "*.mmi", emit: index
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), emit: versions_minimap2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def suffix = minimap2_ps.trim().replace(' ', '').replace(',', '_').replace('-', '_')
    """
    minimap2 ${minimap2_ps} -v 0 -I 20G -t ${task.cpus} -d ${index}${suffix}.mmi ${db}/${index}.fna
    """

    stub:
    def suffix = minimap2_ps.trim().replace(' ', '').replace(',', '_').replace('-', '_')
    """
    touch ${index}${suffix}.mmi
    """
}

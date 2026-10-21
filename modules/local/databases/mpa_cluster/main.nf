// The proteins of every genome clustered at 90% identity over 80% of both
// sequences, the UniRef90/Uniclust90 criteria MetaPhlAn 4's gene families
// follow. A cluster is one gene family; which genomes of a species it reaches
// is that family's coreness. linclust, as HROM used: linear in the number of
// proteins, which is what makes a catalogue of thousands of genomes feasible.
process DB_MPA_CLUSTER {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/mmseqs2:18.8cc5c--hd6d6fdc_0'
        : 'quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0'}"

    input:
    tuple val(meta), path(proteins, stageAs: 'proteins/*')

    output:
    tuple val(meta), path("${meta.id}_cluster.tsv"), emit: clusters
    tuple val("${task.process}"), val('mmseqs2'), eval('mmseqs version'), emit: versions_mmseqs2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat proteins/* > all_proteins.faa
    mmseqs easy-linclust \\
        all_proteins.faa \\
        ${meta.id} \\
        tmp \\
        ${args} \\
        --threads ${task.cpus}
    rm -rf tmp all_proteins.faa ${meta.id}_all_seqs.fasta ${meta.id}_rep_seq.fasta
    """

    stub:
    """
    touch ${meta.id}_cluster.tsv
    """
}

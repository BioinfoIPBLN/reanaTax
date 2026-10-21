// Every marker candidate cut into 150-nt fragments and mapped with Bowtie2
// (-a --very-sensitive, as MetaPhlAn 4 did) against the mapping targets of
// every species and the --db_metaphlan_screen sequences; written down is which
// genomes each candidate reached. The MetaPhlAn container, for its Bowtie2.
// See bin/mpa_uniqueness.py.
process DB_MPA_UNIQUENESS {
    tag "${meta.id}"
    label 'process_high'
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metaphlan:4.2.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/metaphlan:4.2.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(candidates), path(targets), path(genomes, stageAs: 'genomes/*')
    path screen, stageAs: 'screen/*'

    output:
    tuple val(meta), path("hits.tsv"), emit: hits
    tuple val("${task.process}"), val('bowtie2'), eval("bowtie2 --version 2>&1 | head -1 | sed 's/.*version //'"), emit: versions_bowtie2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def screen_arg = screen ? "--screen ${screen}" : ''
    """
    mpa_uniqueness.py prepare \\
        --candidates ${candidates} \\
        --targets ${targets} \\
        --genomes-dir genomes \\
        ${screen_arg} \\
        --fragments-out fragments.fa \\
        --targets-out targets.fa

    bowtie2-build --threads ${task.cpus} targets.fa targets
    bowtie2 \\
        -p ${task.cpus} \\
        -f \\
        -a \\
        --very-sensitive \\
        --no-unal \\
        --no-hd \\
        ${args} \\
        -x targets \\
        -U fragments.fa \\
        | mpa_uniqueness.py hits --output hits.tsv

    rm -f targets.fa targets.*.bt2 targets.*.bt2l fragments.fa
    """

    stub:
    """
    touch hits.tsv
    """
}

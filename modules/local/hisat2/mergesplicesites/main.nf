//
// The splice sites pass 1 of --hisat2_two_pass found, from every library
// aligned against one host reference, as one file for the second pass's
// --known-splicesite-infile.
//
// HISAT2 reads that file with `>>` in a loop that ends only at end of file, so
// a line that is not `name left right strand` would stop the stream short of
// it and hang every second-pass task. Only well-formed lines are kept. A run
// whose libraries span no splice site gives an empty file, which HISAT2 reads
// as no sites.
//
process HISAT2_MERGESPLICESITES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(splicesites, stageAs: 'pass1/*')

    output:
    tuple val(meta), path("*.splicesites.tsv"), emit: splicesites
    tuple val("${task.process}"), val('coreutils'), eval("sort --version | sed -n '1s/.* //p'"), emit: versions_coreutils, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat pass1/* \\
        | { grep -E '^[^[:space:]]+[[:space:]][0-9]+[[:space:]][0-9]+[[:space:]][-+.]\$' || [ \$? -eq 1 ]; } \\
        | LC_ALL=C sort -T . --parallel=${task.cpus} -k1,1 -k2,2n -k3,3n -k4,4 -u \\
        > ${prefix}.splicesites.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.splicesites.tsv
    """
}

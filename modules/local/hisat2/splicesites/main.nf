//
// Pass 1 of --hisat2_two_pass: align a library, or one chunk of it, only to
// collect the splice sites its reads span. The alignments go to /dev/null.
// HISAT2_MERGESPLICESITES joins the sites of every library aligned against the
// same reference, and the real alignment reads them back with
// --known-splicesite-infile.
//
// The sites come from --novel-splicesite-outfile, which HISAT2 fills whether or
// not it may reuse a site within the run, so pass 1 keeps the speed of
// --no-temp-splicesite (added by conf/modules.config when --hisat2_args lacks
// it). HISAT2 writes the file when it finishes; a task stopped early leaves
// none, and fails on the missing output rather than passing on a partial list.
//
process HISAT2_SPLICESITES {
    tag "${meta.chunk ? "${meta.id} ${meta.chunk}" : meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.splicesites.tsv"), emit: splicesites
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), emit: versions_hisat2, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }
    def input = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    """
    INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
    hisat2 \\
        -x \$INDEX \\
        ${input} \\
        ${strandedness} \\
        --threads ${task.cpus} \\
        --novel-splicesite-outfile ${prefix}.splicesites.tsv \\
        ${args} \\
        > /dev/null
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.splicesites.tsv
    """
}

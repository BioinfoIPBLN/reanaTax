// Index the host FASTA, once, for the CRAM writer.
//
// CRAM stores differences from a reference rather than the sequence, so every
// task that writes or reads one needs random access to that reference. htslib
// will build the .fai itself if it is missing - beside whatever path it was
// handed, which under Nextflow is the symlink in the task work directory - so
// omitting this step is correct but pays the build again in every task, on a
// file that is several gigabytes for a human assembly and larger still for the
// two-genome host sets this pipeline is usually run with.
//
// Runs only with --alignment_output_format cram. A BAM run has no use for it.
process SAMTOOLS_FAIDX {
    tag "${fasta}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fai"), emit: fai
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // The staged FASTA is a symlink; samtools names the index after the path it
    // was given, not after the link target, so the .fai lands in this task's
    // directory and the shared reference tree is never written to.
    """
    samtools faidx ${args} ${fasta}
    """

    stub:
    """
    touch ${fasta}.fai
    """
}

//
// The FASTA a HISAT2 index was built from, recovered from the index itself.
//
// A host given as a prebuilt HISAT2 index has no FASTA beside it, and Bowtie2
// cannot read a HISAT2 index. hisat2-inspect prints the reference sequences the
// index holds, so a DNA library can still be depleted against the same host
// without the user having to supply the genome a second time.
//
process HISAT2_INSPECT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(index)

    output:
    tuple val(meta), path("*.fa"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
    hisat2-inspect \$INDEX > ${prefix}.fa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf '>stub\\nACGT\\n' > ${prefix}.fa
    """
}

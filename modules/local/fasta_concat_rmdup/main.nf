//
// Concatenate several host genome FASTA files (each optionally gzipped) into
// one stream and collapse exact duplicate sequences with seqkit rmdup, giving
// a single deduplicated reference that every index-build falls out of.
//
// This is the "combined host" half of --host_bundle: multiple genomes the user
// wants depleted against ONCE, as one genome, rather than as separate passes.
// Duplicates that survive the concatenation are almost always the same
// contig/haplotype present in more than one of the inputs, so collapsing by
// exact sequence (case-insensitively) removes index-tar-like redundancy without
// touching anything a mapper would separately use. Members that are gzipped are
// decompressed into the stream so the result is a plain, concatenable FASTA.
//
process FASTA_CONCAT_RMDUP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.8.2--h9ee0642_0'
        : 'quay.io/biocontainers/seqkit:2.8.2--h9ee0642_0'}"

    input:
    tuple val(meta), path(genomes)

    output:
    tuple val(meta), path("collapsed.fasta"), emit: fasta
    tuple val("${task.process}"), val('seqkit'), eval('seqkit version'), emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    : > combined.fa
    for g in ${genomes.join(' ')}; do
        case "\$g" in
            *.gz) gzip -cdf "\$g" ;;
            *)    cat "\$g" ;;
        esac
    done >> combined.fa
    seqkit rmdup --by-seq --ignore-case ${args} combined.fa -o collapsed.fasta

    # The concatenated stream is only intermediate; the collapsed reference is
    # what every index is built from.
    rm -f combined.fa
    """

    stub:
    """
    touch collapsed.fasta
    """
}

// A host genome with every region that a protected genome also holds replaced
// with N, so that reads of the protected organisms have no host sequence to
// align to and survive depletion. Hostile's masking, step for step: the
// protected genomes cut into 150-mers every 10 bp (bin/host_mask_kmers.py),
// aligned to the host with minimap2 (-N 10000 -p 0.0001, every placement kept),
// the hit intervals merged with bedtools and masked with `bedtools maskfasta`.
// Coordinates are unchanged, so an annotation of the host still applies. Run in
// Hostile's container, which carries exactly these tools.
process HOST_MASK {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(fasta)
    path targets, stageAs: 'target?/*'
    val kmer_length
    val kmer_step

    output:
    tuple val(meta), path("${prefix}.masked.fa"), emit: fasta
    tuple val(meta), path("${prefix}.mask.bed"), emit: bed
    tuple val(meta), path("${prefix}.mask_summary.tsv"), emit: summary
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), emit: versions_minimap2, topic: versions
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed 's/bedtools v//'"), emit: versions_bedtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    host_mask_kmers.py \\
        --targets target*/* \\
        --length ${kmer_length} \\
        --step ${kmer_step} \\
        --output - \\
        | minimap2 -N 10000 -p 0.0001 -t ${task.cpus} ${args} ${fasta} - \\
        | awk -v OFS='\\t' '{ print \$6, \$8, \$9 }' \\
        | sort -k1,1 -k2,2n \\
        | bedtools merge -i stdin \\
        > ${prefix}.mask.bed

    bedtools maskfasta -fi ${fasta} -bed ${prefix}.mask.bed -fo ${prefix}.masked.fa

    awk -v OFS='\\t' -v bed=${prefix}.mask.bed '
        !/^>/ { length_total += length(\$0) }
        END {
            while ((getline line < bed) > 0) { split(line, field, "\\t"); masked += field[3] - field[2]; regions += 1 }
            print "reference", "length", "masked_regions", "masked_positions", "masked_fraction"
            print "${meta.id}", length_total, regions + 0, masked + 0, (length_total ? masked / length_total : 0)
        }' ${fasta} > ${prefix}.mask_summary.tsv
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.masked.fa ${prefix}.mask.bed ${prefix}.mask_summary.tsv
    """
}

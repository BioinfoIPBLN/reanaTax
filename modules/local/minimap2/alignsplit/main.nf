// minimap2 host depletion for long-read libraries, with the contract of
// BOWTIE2_ALIGN_SPLIT: non-host reads as FASTQ, host reads as an optional BAM,
// and a summary in Bowtie2's format, so read accounting and MultiQC read it
// unchanged. Hostile's long-read route, aligned with the preset the index was
// built with (meta2.preset) and --secondary=no.
//
// A read is host when its primary alignment is mapped. Supplementary records
// (the other pieces of a chimeric or split read) are left out of every count
// and of the host BAM, which then holds one record per aligned read, as the
// short-read modules do. In the summary "aligned >1 times" is a primary
// alignment with MAPQ 0, minimap2's sign of equally good placements.
//
// --split-prefix makes a host larger than one index batch (-I, 8 Gb by default)
// come out as one alignment per read: without it each batch reports the read
// separately, and a read unmapped in one batch would be released as non-host
// while aligning in the next.
process MINIMAP2_ALIGN_SPLIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/hostile:2.0.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/hostile:2.0.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val save_unaligned
    val save_bam

    output:
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.summary.log"), emit: summary
    tuple val(meta), path("*fastq.gz"), optional: true, emit: fastq
    tuple val("${task.process}"), val('minimap2'), eval('minimap2 --version'), emit: versions_minimap2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_level = task.ext.compression != null ? "-l ${task.ext.compression} " : ''
    def fq_level = task.ext.compression != null ? "-c ${task.ext.compression} " : ''
    def rg = args.contains('-R ') ? '' : "-R '@RG\\tID:${prefix}\\tSM:${prefix}'"
    def host_bam = save_bam ? "samtools view -b ${bam_level}-F 0x904 tmp/${prefix}.all.bam > ${prefix}.bam" : ''
    def unaligned = save_unaligned ? "samtools fastq ${fq_level}-f 4 -F 0x900 -0 ${prefix}.unmapped.fastq.gz -n tmp/${prefix}.all.bam" : ''
    """
    mkdir -p tmp

    minimap2 \\
        -a \\
        -x ${meta2.preset} \\
        -t ${task.cpus} \\
        --secondary=no \\
        --split-prefix tmp/${prefix}.split \\
        ${rg} \\
        ${args} \\
        ${index} \\
        ${reads} \\
        | samtools view -b ${bam_level}-F 0x100 - > tmp/${prefix}.all.bam

    total=\$(samtools view -c -F 0x900 tmp/${prefix}.all.bam)
    unaligned=\$(samtools view -c -f 4 -F 0x900 tmp/${prefix}.all.bam)
    unique=\$(samtools view -c -F 0x904 -q 1 tmp/${prefix}.all.bam)
    awk -v total=\$total -v unaligned=\$unaligned -v unique=\$unique 'BEGIN {
        aligned = total - unaligned
        scale = total ? 100 / total : 0
        printf "%d reads; of these:\\n", total
        printf "  %d (100.00%%) were unpaired; of these:\\n", total
        printf "    %d (%.2f%%) aligned 0 times\\n", unaligned, unaligned * scale
        printf "    %d (%.2f%%) aligned exactly 1 time\\n", unique, unique * scale
        printf "    %d (%.2f%%) aligned >1 times\\n", aligned - unique, (aligned - unique) * scale
        printf "%.2f%% overall alignment rate\\n", aligned * scale
    }' > ${prefix}.minimap2.summary.log

    ${host_bam}
    ${unaligned}
    rm -f tmp/${prefix}.all.bam
    cat ${prefix}.minimap2.summary.log >&2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz" : ''
    def host_bam = save_bam ? "touch ${prefix}.bam" : ''
    """
    ${unaligned}
    touch ${prefix}.minimap2.summary.log
    ${host_bam}
    """
}

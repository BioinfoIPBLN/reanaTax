// Bowtie2 host depletion for DNA libraries, with the same contract as
// HISAT2_ALIGN_SPLIT: a pair is released only when NEITHER mate aligned.
//
// HISAT2 is splice-aware, which is right for RNA and wrong for DNA: it will
// split a genomic read across an "intron" that is not one, and has no local
// mode. Bowtie2 aligns the read as one stretch, and takes the full preset
// range - --very-sensitive by default, --sensitive-local and the rest through
// --bowtie2_args.
//
// The split is made on SAM flags exactly as in the HISAT2 module, and Bowtie2
// emits both mates of a pair adjacently as HISAT2 does, so `samtools fastq`
// stays in sync. Its summary goes to stderr in the same format HISAT2 writes
// with --summary-file, so the chunk merge, read accounting and MultiQC read it
// unchanged; it is kept as <prefix>.bowtie2.summary.log.
//
// Without a host BAM to keep (`save_bam` false) the stream goes straight into
// `samtools fastq`, as in the HISAT2 module.
process BOWTIE2_ALIGN_SPLIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b41b403e81883126c3227fc45840015538e8e2212f13abc9ae84e4b98891d51c/data' :
        'community.wave.seqera.io/library/bowtie2_htslib_samtools_pigz:edeb13799090a2a6' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val save_unaligned
    val save_bam

    output:
    tuple val(meta), path("*.bam")    , optional:true, emit: bam
    tuple val(meta), path("*.summary.log")           , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val("${task.process}"), val('bowtie2'), eval("bowtie2 --version 2>&1 | sed -n 's/.*bowtie2-align-s version //p'"), emit: versions_bowtie2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_level = task.ext.compression != null ? "-l ${task.ext.compression} " : ''
    def fq_level = task.ext.compression != null ? "-c ${task.ext.compression} " : ''
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"
    def inputs = meta.single_end ? "-U ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    // Unmapped is flag 4 for a single read, 4 and 8 (both mates) for a pair;
    // the host BAM is everything aligned, -F 4 -F 8 for a pair.
    def unmapped = meta.single_end ? '-f 4' : '-f 12'
    def mapped = meta.single_end ? '-F 4' : '-F 4 -F 8'
    def to_fastq = meta.single_end
        ? "samtools fastq ${fq_level}-0 ${prefix}.unmapped.fastq.gz -n -"
        : "samtools fastq ${fq_level}-1 ${prefix}.unmapped_1.fastq.gz -2 ${prefix}.unmapped_2.fastq.gz -0 /dev/null -s /dev/null -n -"

    def sink = "| samtools view -bS ${bam_level}-F 256 - > tmp/${prefix}.all.bam"
    def host_bam = "samtools view -b ${bam_level}${mapped} tmp/${prefix}.all.bam > ${prefix}.bam"
    def unaligned = save_unaligned ? "samtools view -u ${unmapped} tmp/${prefix}.all.bam | ${to_fastq}" : ''
    if (!save_bam) {
        sink = save_unaligned ? "| samtools view -u ${unmapped} - | ${to_fastq}" : '> /dev/null'
        host_bam = ''
        unaligned = ''
    }
    """
    INDEX=`find -L ./ -name "*.rev.1.bt2" | sed "s/\\.rev.1.bt2\$//"`
    [ -z "\$INDEX" ] && INDEX=`find -L ./ -name "*.rev.1.bt2l" | sed "s/\\.rev.1.bt2l\$//"`
    [ -z "\$INDEX" ] && echo "Bowtie2 index files not found" 1>&2 && exit 1
    mkdir -p tmp

    bowtie2 \\
        -x \$INDEX \\
        ${inputs} \\
        --threads $task.cpus \\
        $rg \\
        $args \\
        2> ${prefix}.bowtie2.summary.log \\
        ${sink}

    ${host_bam}
    ${unaligned}
    cat ${prefix}.bowtie2.summary.log >&2
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned
        ? (meta.single_end
            ? "echo '' | gzip > ${prefix}.unmapped.fastq.gz"
            : "echo '' | gzip > ${prefix}.unmapped_1.fastq.gz\n    echo '' | gzip > ${prefix}.unmapped_2.fastq.gz")
        : ''
    def host_bam = save_bam ? "touch ${prefix}.bam" : ''
    """
    ${unaligned}
    touch ${prefix}.bowtie2.summary.log
    ${host_bam}
    """
}

// HISAT2 host depletion that releases a pair only when NEITHER mate aligns.
//
// The nf-core module hardcodes `--no-mixed --no-discordant` and takes its
// non-host reads from `--un-conc-gz`, which is defined as "pairs that failed to
// align CONCORDANTLY". A pair whose R1 is unambiguously host and whose R2 is not
// therefore leaves as non-host, and per-mate alignment status is never even
// computed - `--no-mixed` suppresses it, and the BAM's `-F 4 -F 8` discards the
// unaligned side. Measured on a vulvar-swab RNA library, that leak put 13.5% of
// the surviving pairs back into Kraken2 as Homo sapiens.
//
// So here the two flags are dropped, the whole alignment stream is kept
// (`-F 256` only), and the split is made on explicit flags: `-F 4 -F 8` is the
// host BAM, `-f 12` (both mates unmapped) is the non-host FASTQ. HISAT2 emits
// both mates of a pair adjacently and `samtools view` preserves that order, so
// `samtools fastq` stays in sync without a collate or re-pairing pass.
//
// The whole alignment stream is parked in `tmp/` rather than beside the outputs
// because it has to be read twice (once per side of the split) and the `bam`
// output is declared as a glob. A `*.all.bam` sitting next to `*.bam` matches
// that glob too, and the nf-core SAMTOOLS_SORT `samtools cat`s whatever it is
// handed - which silently published a host BAM holding every unaligned record
// plus a second copy of every aligned one.
process HISAT2_ALIGN_SPLIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/hisat2_samtools:5a258fe6e30b2c20' :
        'community.wave.seqera.io/library/hisat2_samtools:6ca0ef72b662d5c8' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    tuple val(meta3), path(splicesites)
    val save_unaligned

    output:
    tuple val(meta), path("*.bam")                   , emit: bam
    tuple val(meta), path("*.log")                   , emit: summary
    tuple val(meta), path("*fastq.gz"), optional:true, emit: fastq
    tuple val("${task.process}"), val('hisat2'), eval("hisat2 --version | sed -n '1s/.*version //p'"), emit: versions_hisat2, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --compression_level, for every BGZF/gzip stream written below. Applied to
    // the throwaway tmp/*.all.bam too, so the setting means one thing across
    // the module rather than being quietly ignored on the largest file it
    // writes. `samtools view` spells it -l, `samtools fastq` spells it -c.
    def level = task.ext.compression != null ? task.ext.compression : 6

    def strandedness = ''
    if (meta.strandedness == 'forward') {
        strandedness = meta.single_end ? '--rna-strandness F' : '--rna-strandness FR'
    } else if (meta.strandedness == 'reverse') {
        strandedness = meta.single_end ? '--rna-strandness R' : '--rna-strandness RF'
    }
    def ss = "$splicesites" ? "--known-splicesite-infile $splicesites" : ''
    def rg = args.contains("--rg-id") ? "" : "--rg-id ${prefix} --rg SM:${prefix}"

    if (meta.single_end) {
        def unaligned = save_unaligned
            ? "samtools view -u -f 4 tmp/${prefix}.all.bam | samtools fastq -c ${level} -0 ${prefix}.unmapped.fastq.gz -n -"
            : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        mkdir -p tmp
        hisat2 \\
            -x \$INDEX \\
            -U $reads \\
            $strandedness \\
            $ss \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $rg \\
            $args \\
            | samtools view -bS -l ${level} -F 256 - > tmp/${prefix}.all.bam

        samtools view -b -l ${level} -F 4 tmp/${prefix}.all.bam > ${prefix}.bam
        ${unaligned}
        """
    } else {
        def unaligned = save_unaligned
            ? """samtools view -u -f 12 tmp/${prefix}.all.bam \\
            | samtools fastq -c ${level} -1 ${prefix}.unmapped_1.fastq.gz -2 ${prefix}.unmapped_2.fastq.gz -0 /dev/null -s /dev/null -n -"""
            : ''
        """
        INDEX=`find -L ./ -name "*.1.ht2*" | sed 's/\\.1.ht2.*\$//'`
        mkdir -p tmp
        hisat2 \\
            -x \$INDEX \\
            -1 ${reads[0]} \\
            -2 ${reads[1]} \\
            $strandedness \\
            $ss \\
            --summary-file ${prefix}.hisat2.summary.log \\
            --threads $task.cpus \\
            $rg \\
            $args \\
            | samtools view -bS -l ${level} -F 256 - > tmp/${prefix}.all.bam

        samtools view -b -l ${level} -F 4 -F 8 tmp/${prefix}.all.bam > ${prefix}.bam
        ${unaligned}
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def unaligned = save_unaligned ? "echo '' | gzip >  ${prefix}.unmapped_1.fastq.gz \n echo '' | gzip >  ${prefix}.unmapped_2.fastq.gz" : ''
    """
    ${unaligned}

    touch ${prefix}.hisat2.summary.log
    touch ${prefix}.bam
    """
}

// The unmapped (non-host) reads, with their cell barcode carried in the name.
//
// Kraken2 truncates a read name at the first whitespace, and `samtools fastq -T`
// appends its tags after a TAB. Left as-is the barcode would therefore be
// dropped on the way into the classifier, and the cell-by-taxon matrix could
// never be built. So the tags are folded into the name itself, with no
// whitespace, and recovered from Kraken2's read-level output afterwards.
//
// Reads whose barcode STARsolo could not correct to the whitelist carry CB:NA.
// They are kept rather than discarded: they cannot be placed in a cell, but
// they are still part of the sample and belong in the pseudobulk profile. The
// matrix builder is what skips them.
process STARSOLO_UNMAPPED {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path('*.nonhost.fastq.gz'), emit: reads
    tuple val(meta), path('*.barcode_stats.tsv'), emit: stats
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // --compression_level, for the gzip stream below.
    def level = task.ext.compression != null ? "-${task.ext.compression} " : ''
    // -f 4 keeps the unmapped records. STARsolo writes one record per cDNA
    // read, so this is single-end from here on: the barcode read carries no
    // biological sequence and nothing would be gained by classifying it.
    """
    mkdir -p tmp

    samtools view -b -f 4 ${args} ${bam} > tmp/${prefix}.unmapped.bam

    samtools fastq -T CB,UB -n tmp/${prefix}.unmapped.bam 2> tmp/${prefix}.fastq.log \\
        | gawk -v stats=tmp/${prefix}.counts '
            NR % 4 == 1 {
                n = split(\$0, f, "\\t")
                cb = "NA"; ub = "NA"
                for (i = 2; i <= n; i++) {
                    if (f[i] ~ /^CB:Z:/)      cb = substr(f[i], 6)
                    else if (f[i] ~ /^UB:Z:/) ub = substr(f[i], 6)
                }
                total++
                if (cb == "NA" || cb == "-") { nocb++; cb = "NA" } else { seen[cb] = 1 }
                print f[1] "|CB:" cb "|UB:" ub
                next
            }
            { print }
            END {
                printf "reads\\t%d\\nreads_without_barcode\\t%d\\nbarcodes\\t%d\\n", total, nocb, length(seen) > stats
            }' \\
        | gzip ${level}-c > ${prefix}.nonhost.fastq.gz

    {
        echo -e "sample\\t${prefix}"
        cat tmp/${prefix}.counts
    } > ${prefix}.barcode_stats.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip -c > ${prefix}.nonhost.fastq.gz
    touch ${prefix}.barcode_stats.tsv
    """
}

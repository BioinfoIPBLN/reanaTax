// SortMeRNA in annotate-only mode: which of a library's non-host reads are
// rRNA, as a list of read IDs. Nothing is removed - the reads Kraken2
// classified stay exactly as they were, and the IDs are only joined to
// Kraken2's per-read calls afterwards (RRNA_SHARE). Pairs are flagged together
// (`--paired-in`), as Kraken2 classifies a pair as one unit.
process RRNA_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47c387ea4cc97404886b95ece4bd7aecc61aa63ddf93b71a31fabf61beb9b6ae/data' :
        'community.wave.seqera.io/library/sortmerna_gzip:659d4ea003ce0d12' }"

    input:
    tuple val(meta), path(reads)
    path fastas
    path index

    output:
    tuple val(meta), path("*.rrna_ids.txt.gz"), emit: ids
    tuple val(meta), path("*.rrna_annotation.log"), emit: log
    tuple val("${task.process}"), val('sortmerna'), eval('sortmerna --version 2>&1 | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+" | head -1'), topic: versions, emit: versions_sortmerna

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def files = [reads].flatten()
    def paired = files.size() > 1 ? '--paired-in --out2' : ''
    """
    sortmerna \\
        ${fastas.collect { fasta -> "--ref ${fasta}" }.join(' ')} \\
        --idx-dir ${index} \\
        --index 0 \\
        ${files.collect { fastq -> "--reads ${fastq}" }.join(' ')} \\
        --threads ${task.cpus} \\
        --workdir . \\
        --aligned rRNA_reads \\
        --fastx \\
        ${paired} \\
        ${args}

    # Mate 1 carries the pair's ID; --paired-in put both mates of every
    # flagged pair in the aligned output.
    for fastq in rRNA_reads.f*q* rRNA_reads_fwd.f*q*; do
        if [ -e "\$fastq" ]; then gzip -cdf "\$fastq"; fi
    done | awk 'NR % 4 == 1 { sub(/^@/, "", \$1); sub(/\\/[12]\$/, "", \$1); print \$1 }' | gzip -c > ${prefix}.rrna_ids.txt.gz
    rm -f rRNA_reads*.f*q*
    mv rRNA_reads.log ${prefix}.rrna_annotation.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo -n | gzip -c > ${prefix}.rrna_ids.txt.gz
    touch ${prefix}.rrna_annotation.log
    """
}

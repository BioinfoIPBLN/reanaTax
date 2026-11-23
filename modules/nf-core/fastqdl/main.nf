// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process FASTQDL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ed/ed693e48ea3639b0a0f6a0ad191ad31f4e8fa15674a42c4053df6e06eff30cdc/data'
:         'community.wave.seqera.io/library/fastq-dl:4.0.1--0758300a069ba945' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fastq.gz")       , emit: fastq
    tuple val(meta), path("*-run-info.tsv")   , emit: runinfo
    tuple val(meta), path("*-run-mergers.tsv"), emit: runmergers, optional: true
    tuple val(meta), path("*.sra_fallback.txt"), emit: sra_fallback, optional: true
    tuple val("${task.process}"), val('fastq-dl'), eval('fastq-dl --version |& sed "s/.* //"'), emit: versions_fastqdl, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-dl \\
        ${args} \\
        --prefix ${prefix} \\
        --accession ${accession} \\
        --cpus ${task.cpus} \\
        --outdir .

    # fastq-dl before 3.1 exited 0 even when every attempt at a file failed,
    # leaving a zero-byte .fastq.gz that failed several steps later inside
    # whichever tool read it first. 4.x exits non-zero instead; the check stays
    # because a file can still arrive empty or truncated, and failing here makes
    # Nextflow retry the download rather than the step that tripped over it.
    # gzip -t reads the whole stream because a truncated transfer has a valid
    # header and only fails at the end.
    found=0
    for fastq in *.fastq.gz; do
        [ -e "\$fastq" ] || continue
        found=1
        if [ ! -s "\$fastq" ]; then
            echo "ERROR: ${accession}: '\$fastq' is empty - the archive most likely lists a file it cannot serve. Check the other provider (--fastqdl_provider)." >&2
            exit 1
        fi
        if ! gzip -t "\$fastq" 2>/dev/null; then
            echo "ERROR: ${accession}: '\$fastq' is not a complete gzip stream - the download was truncated." >&2
            exit 1
        fi
    done
    if [ "\$found" -eq 0 ]; then
        echo "ERROR: ${accession}: fastq-dl produced no FASTQ at all." >&2
        exit 1
    fi
    ${task.ext.sra_fallback ? "echo '${accession}: fetched from SRA on attempt ${task.attempt}, after ENA failed' > ${prefix}.sra_fallback.txt" : ''}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${accession}.fastq.gz
    echo "" | gzip > ${accession}_1.fastq.gz
    echo "" | gzip > ${accession}_2.fastq.gz
    touch ${prefix}-run-info.tsv
    """
}

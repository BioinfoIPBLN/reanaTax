process FASTQDL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/25de473366050bf652634865e5cf8450e682852a03dc34f77243264794cb989f/data':
        'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3' }"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.fastq.gz")       , emit: fastq
    tuple val(meta), path("*-run-info.tsv")   , emit: runinfo
    tuple val(meta), path("*-run-mergers.tsv"), emit: runmergers, optional: true
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

    # fastq-dl exits 0 even when every attempt at a file failed, leaving a
    # zero-byte .fastq.gz behind. Nothing downstream can tell that apart from a
    # real download, so the run continues and dies several steps later inside
    # whichever tool reads the file first - reporting a corrupt FASTQ rather
    # than a failed download, and never retrying the thing that actually broke.
    #
    # Checking here makes the task fail where the fault is, so Nextflow's retry
    # retries the download. gzip -t reads the whole stream because a truncated
    # transfer has a valid header and only fails at the end.
    found=0
    for fastq in *.fastq.gz; do
        [ -e "\$fastq" ] || continue
        found=1
        if [ ! -s "\$fastq" ]; then
            echo "ERROR: ${accession}: '\$fastq' is empty. fastq-dl exhausted its retries and still exited 0 - the archive most likely lists a file it cannot serve. Check the other provider (--fastqdl_provider)." >&2
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

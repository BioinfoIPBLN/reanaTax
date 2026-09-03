// Delete an intermediate once the pipeline knows nothing else reads it.
//
// Nextflow has no way to say "this file is finished with": every output lives
// until the run ends or `nextflow clean` removes it. For a cohort that is
// downloaded rather than supplied, the raw FASTQ is the largest thing on disk
// and is dead the moment trimming and raw FastQC have both read it - measured
// on one CSI-Microbes run, 413 GB of it.
//
// Ordering is enforced by data, not by hope: the sentinels are the OUTPUTS of
// the tasks that had to finish first, so this task cannot be scheduled until
// they exist. They are staged and never looked at.
//
// The safety rule is the whole design. A staged name is a symlink, so it is the
// LINK TARGET that would have to go, and a target is removed only when it lies
// inside the work directory. A FASTQ named in a samplesheet is the user's file,
// living wherever they put it; it is skipped, reported, and left alone. That
// check is what makes this safe to point at any channel.
process CLEANUP_INTERMEDIATES {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(targets, stageAs: 'targets/*'), path(sentinels, stageAs: 'sentinels/*')

    output:
    tuple val(meta), path("*.cleanup.log"), emit: log

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    : > ${prefix}.cleanup.log
    for staged in targets/*; do
        [ -e "\$staged" ] || continue
        target=\$(readlink -f "\$staged" || true)
        size=\$(stat -Lc%s "\$staged" 2>/dev/null || echo 0)
        case "\$target" in
            ${workflow.workDir}/*)
                rm -f "\$target"
                echo -e "removed\\t\${size}\\t\$target" >> ${prefix}.cleanup.log
                ;;
            *)
                echo -e "kept\\t\${size}\\t\$target\\t(outside the work directory)" >> ${prefix}.cleanup.log
                ;;
        esac
    done
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cleanup.log
    """
}

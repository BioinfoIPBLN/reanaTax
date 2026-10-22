// The cohort's Emu profiles as one table per rank, with `emu combine-outputs`:
// relative abundances, and estimated counts too when EMU_ABUNDANCE kept them
// (--keep-counts). combine-outputs reads every file in a directory whose name
// holds `rel-abundance`, which the -threshold- tables EMU_ABUNDANCE can also
// write would match as extra samples, so only the full profiles are staged.
process EMU_COMBINE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/emu:3.6.2--hdfd78af_0'
        : 'quay.io/biocontainers/emu:3.6.2--hdfd78af_0'}"

    input:
    tuple val(meta), path(profiles, stageAs: 'profiles/*')
    val ranks
    val counts

    output:
    tuple val(meta), path("emu-combined-*.tsv"), emit: combined
    tuple val("${task.process}"), val('emu'), eval("emu --version 2>&1 | sed 's/^emu v//'"), emit: versions_emu, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def count_tables = counts ? 'emu combine-outputs profiles "$rank" --counts' : ''
    """
    for rank in ${ranks.join(' ')}; do
        emu combine-outputs profiles "\$rank"
        ${count_tables}
    done
    mv profiles/emu-combined-*.tsv .
    """

    stub:
    """
    touch emu-combined-${ranks[0]}.tsv
    """
}

// sylph-tax: turn sylph's genome-level rows into a taxonomic profile.
//
// sylph reports genome files, not taxa - a database is a bag of FASTA
// sketches - so the names come from a metadata table that maps each genome
// file to a lineage. It has to be the table built for the database that was
// profiled against: sylph-tax matches on the genome file name.
//
// Files, not sylph-tax's built-in identifiers (GTDB_r232, IMGVR_4.1, ...).
// Those resolve through a folder that `sylph-tax download` records in a config
// under $HOME, which a container does not have; the same tables are published
// with the databases and are passed here directly.
//
// SYLPH_TAXONOMY_CONFIG points into the task directory rather than /tmp. /tmp is
// shared by every task on the node, and a fixed path there is exactly how the
// Kraken2 daemon wedged.
//
// A sample with nothing detected is a legitimate outcome - on a transcriptome,
// a common one - so it gets an empty profile instead of a call into sylph-tax
// with no rows, whose behaviour on that input nothing downstream should depend
// on. Anything else that fails to produce exactly one profile is an error.
process SYLPH_TAXPROF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sylph-tax:1.9.2--pyhdfd78af_0'
        : 'quay.io/biocontainers/sylph-tax:1.9.2--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(profile)
    path taxonomy, stageAs: 'taxonomy/*'

    output:
    tuple val(meta), path("*.sylphmpa"), emit: profile
    tuple val("${task.process}"), val('sylph-tax'), eval("sylph-tax --version 2>&1 | tail -1"), emit: versions_sylphtax, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export SYLPH_TAXONOMY_CONFIG="\$PWD/sylph-tax-config.json"

    if [ "\$(awk 'NR > 1' ${profile} | wc -l)" -eq 0 ]; then
        printf '#SampleID\\t%s\\n' "${prefix}" > ${prefix}.sylphmpa
        printf 'clade_name\\trelative_abundance\\tsequence_abundance\\tANI (if strain-level)\\tCoverage (if strain-level)\\n' >> ${prefix}.sylphmpa
    else
        sylph-tax taxprof ${profile} ${args} -t taxonomy/*
        shopt -s nullglob
        made=( *.sylphmpa )
        if [ "\${#made[@]}" -ne 1 ]; then
            echo "ERROR: ${meta.id}: expected one .sylphmpa from sylph-tax, found \${#made[@]}." >&2
            exit 1
        fi
        [ "\${made[0]}" = "${prefix}.sylphmpa" ] || mv "\${made[0]}" ${prefix}.sylphmpa
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sylphmpa
    """
}

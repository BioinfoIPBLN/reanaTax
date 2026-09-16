//
// Score every taxon against the assembly.
//
// The one line of evidence this pipeline otherwise has none of. The minimizer,
// host-k-mer, shuffle, decontam and negative-control filters all judge a taxon
// by counting things about its reads; not one of them ever produces a longer
// sequence. A contig is orthogonal to all of them: reads placed by index
// hopping, or piled on one conserved locus, do not assemble into anything.
//
// Like the other evidence filters it removes nothing itself. The taxids go to
// ABUNDANCE_FILTER so one step owns every removal from the combined tables, and
// the read gate means a taxon too rare to have assembled is recorded rather
// than condemned.
//
process CONTIG_EVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(table, stageAs: 'table/*')
    path assignments, stageAs: 'assignments/*'
    path reports, stageAs: 'reports/*'
    val min_reads
    val min_length
    val min_contigs

    output:
    tuple val(meta), path("*.contig_evidence.tsv"), emit: evidence
    tuple val(meta), path("*.contig_drop.txt")    , emit: drop_list
    tuple val(meta), path("*_mqc.tsv")            , emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    contig_evidence.py \\
        --table ${table} \\
        --assignments assignments/* \\
        --reports reports/* \\
        --output ${prefix}.contig_evidence.tsv \\
        --drop-list ${prefix}.contig_drop.txt \\
        --mqc ${prefix}.contig_evidence_mqc.tsv \\
        --min-reads ${min_reads} \\
        --min-length ${min_length} \\
        --min-contigs ${min_contigs} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.contig_evidence.tsv
    touch ${prefix}.contig_drop.txt
    """
}

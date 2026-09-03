// Which host cell types carry which microbes.
//
// The question the single-cell branch exists to answer, and the one step whose
// statistics come from CSI-Microbes rather than SAHMI: one-sided Fisher exact
// per sample, combined across samples by Stouffer's Z weighted by the expected
// number of infected cells, with log2(summed observed / summed expected) as the
// effect size. See bin/sc_enrichment.py for why Fisher rather than chi-square,
// presence rather than abundance, and why pooling samples is refused outright.
//
// One task for the cohort, because combining across samples IS the analysis.
// It needs cell-type annotations, which the pipeline cannot produce - STARsolo
// gives a count matrix, not clusters - so --sc_cell_metadata is required and
// the step is skipped without it.
process SCTAXA_ENRICHMENT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/python:3.12'
        : 'biocontainers/python:3.12'}"

    input:
    tuple val(meta), path(counts, stageAs: 'counts/*')
    path cell_metadata
    val min_cells
    val p_threshold

    output:
    tuple val(meta), path("*.cell_type_enrichment.tsv"), emit: enrichment
    tuple val(meta), path("*.cooccurrence.tsv"), emit: cooccurrence, optional: true
    tuple val(meta), path("*.doublet_check.tsv"), emit: doublet, optional: true
    tuple val(meta), path("*_mqc.tsv"), emit: mqc, optional: true
    tuple val("${task.process}"), val('python'), eval('python3 --version | sed "s/Python //"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sc_enrichment.py \\
        ${args} \\
        --counts counts/* \\
        --cells ${cell_metadata} \\
        --min-cells ${min_cells} \\
        --p-threshold ${p_threshold} \\
        --enrichment ${prefix}.cell_type_enrichment.tsv \\
        --cooccurrence ${prefix}.cooccurrence.tsv \\
        --doublet-test ${prefix}.doublet_check.tsv \\
        --mqc ${prefix}_sc_enrichment_mqc.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cell_type_enrichment.tsv ${prefix}.cooccurrence.tsv ${prefix}.doublet_check.tsv
    """
}

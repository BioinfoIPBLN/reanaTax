//
// Differential expression of host gene counts with DESeq2.
//
// One process per method because each needs its own container, exactly as the
// microbial side does. See bin/gene_de.R for why gene counts get a negative
// binomial model rather than the compositional one ALDEx2/ANCOM-BC2 apply.
//
// The container is pinned to an r44 build deliberately: the newer r45 builds of
// both DESeq2 and edgeR exist on quay.io but are NOT on the Galaxy singularity
// depot, so a singularity run would fail to pull them.
//
process GENE_DE_DESEQ2 {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-deseq2:1.46.0--r44he5774e6_1'
        : 'quay.io/biocontainers/bioconductor-deseq2:1.46.0--r44he5774e6_1'}"

    input:
    tuple val(meta), path(counts)
    path metadata

    output:
    tuple val(meta), path("*.results.tsv"), emit: results
    tuple val(meta), path("*_mqc.tsv")    , emit: mqc, optional: true
    tuple val("${task.process}"), val('DESeq2'), eval("Rscript -e 'cat(as.character(packageVersion(\"DESeq2\")))'"), emit: versions_r, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gene_de.R \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --method deseq2 \\
        --output ${prefix}.deseq2.results.tsv \\
        --mqc ${prefix}.deseq2_mqc.tsv \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.deseq2.results.tsv ${prefix}.deseq2_mqc.tsv
    """
}

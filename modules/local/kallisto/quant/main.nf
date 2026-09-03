//
// kallisto quant: pseudoalignment and abundance estimation against a
// transcriptome, without ever producing an alignment.
//
// For a single-end library kallisto cannot infer the fragment-length
// distribution from the data - there is only one end - so -l and -s must be
// supplied, and the abundances inherit whatever error those estimates carry.
// Paired-end input needs neither.
//
// The outputs are renamed out of kallisto's fixed directory layout and given
// the sample id, because KALLISTO_MERGE recovers the sample from the filename
// when it stages many of them side by side.
//
process KALLISTO_QUANT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kallisto:0.51.1--ha4fb952_1'
        : 'quay.io/biocontainers/kallisto:0.51.1--ha4fb952_1'}"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.abundance.tsv"), emit: abundance
    tuple val(meta), path("*.run_info.json"), emit: run_info
    tuple val("${task.process}"), val('kallisto'), eval("kallisto version | sed 's/.*version //'"), emit: versions_kallisto, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "${reads}" : "${reads[0]} ${reads[1]}"
    """
    kallisto quant \\
        -i ${index} \\
        -o quant \\
        --threads ${task.cpus} \\
        ${args} \\
        ${input_reads}

    mv quant/abundance.tsv ${prefix}.abundance.tsv
    mv quant/run_info.json ${prefix}.run_info.json
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'target_id\\tlength\\teff_length\\test_counts\\ttpm\\n' > ${prefix}.abundance.tsv
    echo '{}' > ${prefix}.run_info.json
    """
}

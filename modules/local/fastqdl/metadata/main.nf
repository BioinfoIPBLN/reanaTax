// Resolve an ENA/SRA/BioProject accession into the list of sequencing runs it
// contains, WITHOUT downloading any reads.
//
// This exists so that a single umbrella accession (e.g. a BioProject with 200
// runs) can be fanned out into 200 independent, parallel `FASTQDL` download
// tasks instead of being downloaded serially inside one long-running job.
// It deliberately runs the exact same `fastq-dl` build as the nf-core
// `FASTQDL` module so that the run list and the downloads can never disagree.
process FASTQDL_METADATA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/25/25de473366050bf652634865e5cf8450e682852a03dc34f77243264794cb989f/data'
        : 'community.wave.seqera.io/library/fastq-dl:3.0.1--fa446f61dfc85bc3'}"

    input:
    tuple val(meta), val(accession)

    output:
    tuple val(meta), path("*.runsheet.csv"), emit: runsheet
    tuple val(meta), path("*-run-info.tsv"), emit: runinfo
    tuple val("${task.process}"), val('fastq-dl'), eval('fastq-dl --version |& sed "s/.* //"'), emit: versions_fastqdl, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fastq-dl \\
        ${args} \\
        --only-download-metadata \\
        --prefix ${prefix} \\
        --accession ${accession} \\
        --outdir .

    runinfo_to_runsheet.py \\
        ${prefix}-run-info.tsv \\
        ${prefix}.runsheet.csv \\
        --query ${accession}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf 'run_accession\\texperiment_accession\\tsample_accession\\tlibrary_layout\\n' > ${prefix}-run-info.tsv
    printf '%s\\t%s\\t%s\\tPAIRED\\n' SRR0000001 SRX0000001 SAMN0000001 >> ${prefix}-run-info.tsv

    runinfo_to_runsheet.py \\
        ${prefix}-run-info.tsv \\
        ${prefix}.runsheet.csv \\
        --query ${accession}
    """
}

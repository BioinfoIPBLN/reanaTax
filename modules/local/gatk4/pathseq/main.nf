// GATK PathSeq: a third classification route, and a different kind of one.
//
// Kraken2 and KrakenUniq classify by exact k-mer match against a database of
// whole genomes; MetaPhlAn matches against clade-specific markers. PathSeq
// ALIGNS the surviving reads to its microbe reference with BWA and scores taxa
// from those alignments, distributing an ambiguously-mapping read's weight over
// the taxa it hits instead of pushing it up to their common ancestor.
//
// That difference is the reason to run it. A read that Kraken2 can only place
// at a genus - because its k-mers are shared - contributes nothing at species
// level; PathSeq gives it fractional weight at every species it aligns to, and
// reports `unambiguous` separately so the two kinds of evidence stay
// distinguishable. It is also alignment-based, so it is not fooled by the
// composition-driven chance matches the shuffled-read control exists to
// measure. Where PathSeq and Kraken2 agree on a taxon, two quite different
// failure modes have both been avoided.
//
// The cost is the resource bundle: a BWA index image and a k-mer file for the
// host, a BWA index image, sequence dictionary and taxonomy database for the
// microbes. Broad distributes prebuilt ones. This pipeline does not build them
// - they are large, slow and specific to a reference choice the user has to
// make - so they are parameters.
//
// Runs on the reads that survived host depletion, like every other classifier
// here. PathSeq's own host filter then acts as a second, alignment-plus-k-mer
// pass over what HISAT2 left behind, which is the same layering as pairing
// GRCh38 with T2T-CHM13.
process GATK4_PATHSEQ {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0'
        : 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'}"

    input:
    tuple val(meta), path(bam)
    path microbe_bwa_image
    path microbe_dict
    path taxonomy_db
    path host_bwa_image
    path host_kmers
    val save_bam

    output:
    tuple val(meta), path("*.pathseq.scores.txt"), emit: scores
    tuple val(meta), path("*.pathseq.bam"), emit: bam, optional: true
    tuple val(meta), path("*.filter_metrics.txt"), emit: filter_metrics, optional: true
    tuple val(meta), path("*.score_metrics.txt"), emit: score_metrics, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version 2>&1 | sed -n 's/^The Genome Analysis Toolkit (GATK) v//p'"), emit: versions_gatk, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def memory = task.memory ? "-Xmx${task.memory.toGiga()}g" : ''
    // Both host resources or neither: the k-mer file is what makes the BWA
    // filter affordable, and PathSeq's filter stage expects them together.
    // Given neither, PathSeq skips host filtering entirely - which is a
    // reasonable thing to ask of it here, because HISAT2 has already done that
    // job and the reads reaching this point are the ones it could not place.
    def host_args = host_bwa_image && host_kmers
        ? "--filter-bwa-image ${host_bwa_image} --kmer-file ${host_kmers}"
        : ''
    // The output BAM holds every read PathSeq aligned, tagged with its call. It
    // is the size of the library and is only worth keeping when the individual
    // alignments are going to be looked at.
    def keep_bam = save_bam ? "--output ${prefix}.pathseq.bam" : ''
    """
    gatk --java-options "${memory}" PathSeqPipelineSpark \\
        ${args} \\
        --input ${bam} \\
        ${host_args} \\
        --microbe-bwa-image ${microbe_bwa_image} \\
        --microbe-dict ${microbe_dict} \\
        --taxonomy-file ${taxonomy_db} \\
        --scores-output ${prefix}.pathseq.scores.txt \\
        ${keep_bam} \\
        --filter-metrics ${prefix}.filter_metrics.txt \\
        --score-metrics ${prefix}.score_metrics.txt \\
        --spark-master local[${task.cpus}]
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pathseq.scores.txt ${prefix}.filter_metrics.txt ${prefix}.score_metrics.txt
    """
}

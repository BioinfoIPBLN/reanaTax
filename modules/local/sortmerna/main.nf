// nf-core's SORTMERNA module, kept here to run SortMeRNA 8 rather than nf-core's
// 4.3.7. Two things changed under it. Since 6.0, `--index 1` means "build the
// index if it is missing, then align" and index-only is `--task 5`, so that is
// what switches this module into index mode. And 8.0 spells the paired flag
// `--paired-in`; the underscore form is only an alias now. Bioconda stopped
// at 4.3.7, so the image is built from conda-forge's package.
process SORTMERNA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/47c387ea4cc97404886b95ece4bd7aecc61aa63ddf93b71a31fabf61beb9b6ae/data' :
        'community.wave.seqera.io/library/sortmerna_gzip:659d4ea003ce0d12' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(fastas)
    tuple val(meta3), path(index)

    output:
    tuple val(meta), path("*non_rRNA.fastq.gz"), emit: reads, optional: true
    tuple val(meta), path("*.log")             , emit: log, optional: true
    tuple val(meta2), path("idx")              , emit: index, optional: true
    tuple val("${task.process}"), val('sortmerna'), eval('sortmerna --version 2>&1 | grep -oE "[0-9]+\\.[0-9]+\\.[0-9]+" | head -1'), topic: versions, emit: versions_sortmerna

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args  ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    def index_only    = args.contains('--task 5')? true : false
    def skip_index    = args.contains('--index 0')? true : false
    def paired_end    = reads instanceof List
    def paired_cmd    = ''
    def reads_args    = ''
    def out2_cmd      = ''
    def mv_cmd        = ''
    def reads_input   = ''
    def refs_input    = (skip_index && index) ? "--idx-dir ${index}" : ''

    if (! index_only){
        reads_args = '--aligned rRNA_reads --fastx --other non_rRNA_reads'
        reads_input = paired_end ? reads.collect{ r -> "--reads $r"}.join(' ') : "--reads $reads"
        def n_fastq = paired_end ? reads.size() : 1
        if ( n_fastq == 1 ) {
            mv_cmd = """
            mv non_rRNA_reads.f*q.gz ${prefix}.non_rRNA.fastq.gz
            mv rRNA_reads.log ${prefix}.sortmerna.log
            """
        } else {
            mv_cmd = """
            mv non_rRNA_reads_fwd.f*q.gz ${prefix}_1.non_rRNA.fastq.gz
            mv non_rRNA_reads_rev.f*q.gz ${prefix}_2.non_rRNA.fastq.gz
            mv rRNA_reads.log ${prefix}.sortmerna.log
            """
            paired_cmd = "--paired-in"
            out2_cmd   = "--out2"
        }
    }
    """
    sortmerna \\
        ${'--ref '+fastas.join(' --ref ')} \\
        $refs_input \\
        $reads_input \\
        --threads $task.cpus \\
        --workdir . \\
        $reads_args \\
        $paired_cmd \\
        $out2_cmd \\
        $args

    $mv_cmd
    """

    stub:
    def args          = task.ext.args  ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"

    def index_only    = args.contains('--task 5')? true : false
    def paired_end    = reads instanceof List
    def mv_cmd        = ''

    if (! index_only){
        def n_fastq = paired_end ? reads.size() : 1
        if ( n_fastq == 1 ) {
            mv_cmd = "echo | gzip > ${prefix}.non_rRNA.fastq.gz"
        } else {
            mv_cmd = """
            echo "" | gzip > ${prefix}_1.non_rRNA.fastq.gz
            echo "" | gzip > ${prefix}_2.non_rRNA.fastq.gz
            """
        }
    }
    """
    $mv_cmd
    mkdir -p idx
    touch ${prefix}.sortmerna.log
    """
}

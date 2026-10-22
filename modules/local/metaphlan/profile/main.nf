// MetaPhlAn marker-gene profiling, as an alternative to Kraken2/Bracken.
//
// The two answer the same question differently: Kraken2 assigns every read by
// k-mer and Bracken redistributes what it could only place above species, so
// the output is read counts over whatever the database contains. MetaPhlAn maps
// against clade-specific marker genes only, so it reports relative abundance
// over a curated set and never sees most reads. Running both and comparing is
// the point - a taxon only one of them calls is worth a second look.
//
// The invocation mirrors ezRun's app-metaphlan.R (FGCZ), which is what the
// databases under /srv/GT/databases/metaphlan_databases are laid out for:
// --db_dir plus --index naming a bowtie2 index basename with a sibling .pkl.
//
// A long-read library (meta.long_reads) is mapped here rather than by
// MetaPhlAn, against the index METAPHLAN_MM2INDEX built once for the run (see
// there for why), with the command MetaPhlAn's --long_reads mode runs:
// read_fastx.py, then minimap2 with --metaphlan_minimap2_ps, `-v 0
// --sam-hit-only -a`. The SAM goes to MetaPhlAn as `--input_type sam`, which
// applies the same hit filters as its FASTQ path, with the base count
// read_fastx.py reports as --nbases. The read-stats columns then count bases.
process METAPHLAN_PROFILE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/metaphlan:4.2.6--pyhdfd78af_0'
        : 'quay.io/biocontainers/metaphlan:4.2.6--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(reads)
    path db
    val index
    path mm2_index
    val minimap2_ps

    output:
    tuple val(meta), path("*_metaphlan.txt")     , emit: profile
    tuple val(meta), path("*.bowtie2.bz2")       , emit: mapout, optional: true
    tuple val(meta), path("*.metaphlan.log")     , emit: log
    tuple val("${task.process}"), val('metaphlan'), eval("metaphlan --version | sed -n 's/.*version //p'"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.long_reads) {
        def given = args.tokenize(' ')
        def at = given.indexOf('--read_min_len')
        def min_len = at >= 0 && at + 1 < given.size() ? given[at + 1] : '70'
        """
        read_fastx.py -l ${min_len} ${reads} 2> read_fastx.stats \\
            | minimap2 ${minimap2_ps} -v 0 --sam-hit-only -a -t ${task.cpus} ${mm2_index} - \\
            > ${prefix}.minimap2.sam

        nbases=\$(awk 'NR == 1 { printf "%d", \$3 }' read_fastx.stats)

        metaphlan \\
            ${prefix}.minimap2.sam \\
            --input_type sam \\
            --long_reads \\
            --nbases \$nbases \\
            --db_dir ${db} \\
            --index ${index} \\
            --nproc ${task.cpus} \\
            --tmp_dir . \\
            -o ${prefix}_metaphlan.txt \\
            ${args} \\
            1> ${prefix}.metaphlan.log 2>&1

        rm -f ${prefix}.minimap2.sam
        """
    }
    else {
        // MetaPhlAn takes both mates as one comma-separated argument; it treats
        // them as independent reads rather than as a pair.
        def read_arg = meta.single_end ? "${reads}" : "${reads[0]},${reads[1]}"
        """
        metaphlan \\
            ${read_arg} \\
            --db_dir ${db} \\
            --index ${index} \\
            --input_type fastq \\
            --nproc ${task.cpus} \\
            --tmp_dir . \\
            --mapout ${prefix}.bowtie2.bz2 \\
            -o ${prefix}_metaphlan.txt \\
            ${args} \\
            1> ${prefix}.metaphlan.log 2>&1
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_metaphlan.txt ${prefix}.metaphlan.log
    """
}

// Barcode-aware host alignment with STARsolo.
//
// This is the step that makes a single-cell branch possible at all, and the
// reason HOST_DEPLETION_HISAT2 cannot simply be reused. In bulk, host depletion
// DISCARDS the host reads. In single cell the host reads are the other half of
// the experiment - the cell-by-gene matrix - and the cell barcode is the only
// thing that joins the two halves. HISAT2 knows nothing about barcodes, so the
// aligner has to be one that emits the count matrix and the unmapped reads with
// their CB/UB tags still attached, in one pass. That is STARsolo, and it is
// what SAHMI uses for the same reason.
//
// Three settings are forced rather than left to ext.args, because the rest of
// the branch is meaningless without them:
//
//   --outSAMunmapped Within   keeps unmapped reads IN the BAM. Without it they
//                             are dropped and there is nothing to classify.
//   CB/UB in --outSAMattributes  the corrected barcode and deduplicated UMI.
//                             CR/UR (raw) are kept alongside for diagnostics.
//   --outSAMtype BAM SortedByCoordinate  NOT a preference. STAR refuses CB/UB
//                             on an unsorted BAM - "CB and/or UB attributes in
//                             --outSAMattributes can only be output in the
//                             sorted BAM file" - because the corrected barcode
//                             and the deduplicated UMI are only known once the
//                             solo pass has finished, which STAR does during
//                             the sort. Sorting reads that are wanted only as
//                             FASTQ does look like wasted work, and it is; it
//                             is also the only way to get the tags that make
//                             the branch work at all.
//
// --limitBAMsortRAM is therefore set too, and scaled from task.memory rather
// than left at STAR's default of 0 (which means "the genome index size", ~30 GB
// for human, on top of the ~30 GB the index already occupies - close enough to
// the 72 GB of process_high to matter, and it would not grow on retry).
//
// Taking CB rather than re-deriving the barcode from the read is a deliberate
// improvement on SAHMI, whose sckmer.r reads it positionally as
// substr(R1, 1, cb_len) with no whitelist and no error correction - so a single
// sequencing error in the barcode manufactures a new "cell".
process STARSOLO {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/26/268b4c9c6cbf8fa6606c9b7fd4fafce18bf2c931d1a809a0ce51b105ec06c89d/data'
        : 'community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4'}"

    input:
    tuple val(meta), path(reads, stageAs: 'input/*')
    tuple val(meta2), path(index)
    tuple val(meta3), path(gtf)
    path whitelist

    output:
    tuple val(meta), path('*.Aligned.sortedByCoord.out.bam'), emit: bam
    tuple val(meta), path('*Solo.out'), emit: solo
    tuple val(meta), path('*Log.final.out'), emit: log_final
    tuple val(meta), path('*Log.out'), emit: log_out
    tuple val("${task.process}"), val('star'), eval('STAR --version | sed -e "s/STAR_//g"'), emit: versions_star, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // STARsolo's --readFilesIn takes the cDNA read FIRST and the barcode read
    // SECOND, which is the opposite of the R1/R2 order the files are named in.
    // Getting this backwards produces a run that completes with almost nothing
    // aligned, so it is derived here rather than left to the caller.
    def files = reads instanceof List ? reads : [reads]
    if (files.size() != 2) {
        error("STARSOLO needs exactly two FASTQ files for '${meta.id}' (barcode+UMI read and cDNA read) but got ${files.size()}.")
    }
    def barcode_read = files[0]
    def cdna_read = files[1]
    def whitelist_arg = whitelist ? "${whitelist}" : 'None'
    def zipped = files.every { entry -> entry.name.endsWith('.gz') }
    def read_command = zipped ? '--readFilesCommand zcat' : ''
    def sort_ram = (task.memory.toBytes() / 2.5) as long
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${index} \\
        --sjdbGTFfile ${gtf} \\
        --readFilesIn ${cdna_read} ${barcode_read} \\
        ${read_command} \\
        --runThreadN ${task.cpus} \\
        --soloCBwhitelist ${whitelist_arg} \\
        --outSAMunmapped Within \\
        --outSAMtype BAM SortedByCoordinate \\
        --limitBAMsortRAM ${sort_ram} \\
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN \\
        --outFileNamePrefix ${prefix}. \\
        ${args}

    mv ${prefix}.Solo.out ${prefix}_Solo.out
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.Aligned.sortedByCoord.out.bam ${prefix}.Log.final.out ${prefix}.Log.out
    mkdir -p ${prefix}_Solo.out
    """
}

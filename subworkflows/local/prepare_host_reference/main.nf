//
// Make an index of the host genome available, building or downloading it only
// when the user has not supplied one.
//
// Three aligners, because the routes need different things from the host. The
// bulk route depletes RNA with HISAT2 and DNA with Bowtie2 (--host_aligner),
// and throws the host reads away. The single-cell route aligns with STARsolo,
// which has to produce the cell-by-gene matrix as well as the unmapped reads,
// and needs its own index.
//
// On the bulk route each index is built only once a library that needs it
// turns up (`ch_aligners` carries every library's aligner), so a DNA-only run
// never pays for a splice-aware HISAT2 build, nor an RNA-only one for Bowtie2. Everything
// before the index build - resolving --host to a FASTA, downloading from NCBI,
// gunzipping - is identical, which is why this is one subworkflow with a switch
// rather than two that would drift apart.
//

include { NCBIGENOMEDOWNLOAD } from '../../../modules/nf-core/ncbigenomedownload/main'
include { GUNZIP             } from '../../../modules/nf-core/gunzip/main'
include { UNTAR              } from '../../../modules/nf-core/untar/main'
include { HISAT2_BUILD       } from '../../../modules/nf-core/hisat2/build/main'
include { HISAT2_EXTRACTSPLICESITES } from '../../../modules/nf-core/hisat2/extractsplicesites/main'
include { STAR_GENOMEGENERATE       } from '../../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_INSPECT            } from '../../../modules/local/hisat2/inspect/main'
include { BOWTIE2_BUILD             } from '../../../modules/local/bowtie2/build/main'
include { UNTAR as UNTAR_BOWTIE2    } from '../../../modules/nf-core/untar/main'

workflow PREPARE_HOST_REFERENCE {

    take:
    fasta // string: path to a host genome FASTA (optionally gzipped), or null
    hisat2_index // string: path to a prebuilt HISAT2 index directory or tarball, or null
    host_accession // string: NCBI assembly accession (GCF_*/GCA_*), or null
    host_taxid // string: NCBI taxonomy ID, or null
    ncbi_group // string: ncbi-genome-download taxonomic group to search
    gtf // string: path to a GTF for splice-aware index building, or null
    aligner // string: 'hisat2' or 'star' - which index to make available
    star_index // string: path to a prebuilt STAR index directory or tarball, or null
    bowtie2_index // string: path to a prebuilt Bowtie2 index directory or tarball, or null
    ch_aligners // channel: 'hisat2' or 'bowtie2', one per library - which indexes the bulk route needs

    main:

    def prebuilt = aligner == 'star' ? star_index : hisat2_index

    def ch_gtf = gtf
        ? channel.value([[id: file(gtf).baseName], file(gtf, checkIfExists: true)])
        : channel.value([[:], []])

    def ch_index = channel.empty()
    def ch_fasta = channel.empty()

    if (prebuilt) {
        //
        // A prebuilt index short-circuits everything else.
        //
        if (prebuilt.endsWith('.tar.gz') || prebuilt.endsWith('.tgz')) {
            UNTAR(channel.value([[id: "${aligner}_index"], file(prebuilt, checkIfExists: true)]))
            ch_index = UNTAR.out.untar
        }
        else {
            ch_index = channel.value([[id: "${aligner}_index"], file(prebuilt, checkIfExists: true)])
        }
        ch_fasta = fasta
            ? channel.value([[id: file(fasta).baseName], file(fasta, checkIfExists: true)])
            : channel.value([[:], []])
    }
    else {
        //
        // Otherwise obtain a FASTA, either locally or from NCBI.
        //
        def ch_fasta_maybe_gz = channel.empty()

        if (fasta) {
            ch_fasta_maybe_gz = channel.value([[id: file(fasta).baseName], file(fasta, checkIfExists: true)])
        }
        else {
            // ncbi-genome-download reads its accessions/taxids from files, so
            // materialise the single value the user gave us into one.
            def query = host_accession ?: host_taxid
            def is_accession = host_accession as boolean
            def ch_query_file = channel
                .of(query.toString())
                .collectFile(name: is_accession ? 'accessions.txt' : 'taxids.txt', newLine: true)

            // GCA_* assemblies live in GenBank, GCF_* in RefSeq. Getting this
            // wrong is the single most common ncbi-genome-download failure.
            def section = is_accession && host_accession.toString().startsWith('GCA_') ? 'genbank' : 'refseq'
            def meta = [id: query.toString().replaceAll(/[^A-Za-z0-9._-]/, '_'), section: section]

            NCBIGENOMEDOWNLOAD(
                channel.value(meta),
                is_accession ? ch_query_file : channel.value([]),
                is_accession ? channel.value([]) : ch_query_file,
                ncbi_group,
            )

            ch_fasta_maybe_gz = NCBIGENOMEDOWNLOAD.out.fna.map { meta_fna, fna ->
                def files = fna instanceof List ? fna : [fna]
                if (files.size() > 1) {
                    error("ncbi-genome-download returned ${files.size()} genomes for '${query}' (${files*.name.join(', ')}). Narrow the query with --host_accession, or pass --fasta directly.")
                }
                [meta_fna, files.first()]
            }
        }

        def ch_branched = ch_fasta_maybe_gz.branch { _meta, genome ->
            gz: genome.name.endsWith('.gz')
            plain: true
        }

        // hisat2-build cannot read gzipped FASTA.
        GUNZIP(ch_branched.gz)
        ch_fasta = ch_branched.plain.mix(GUNZIP.out.gunzip)

        if (aligner == 'star') {
            // STARsolo needs the annotation in the index: it assigns reads to
            // genes at alignment time, so a --sjdbGTFfile given later cannot
            // recover a matrix built without one.
            STAR_GENOMEGENERATE(ch_fasta, ch_gtf)
            ch_index = STAR_GENOMEGENERATE.out.index
        }
        else {
            // hisat2-build rejects --exon without --ss ('Nongraph exception'), and
            // the build module derives --exon from the GTF but takes --ss as a
            // separate file it does not produce itself. So whenever a GTF is given,
            // the splice sites have to be extracted first.
            def ch_splicesites = channel.value([[:], []])
            if (gtf) {
                HISAT2_EXTRACTSPLICESITES(ch_gtf)
                ch_splicesites = HISAT2_EXTRACTSPLICESITES.out.txt
            }

            HISAT2_BUILD(
                ch_fasta.combine(ch_aligners.filter { name -> name == 'hisat2' }.first()).map { meta, genome, _name -> [meta, genome] },
                ch_gtf,
                ch_splicesites,
            )
            ch_index = HISAT2_BUILD.out.index
        }
    }

    //
    // Bowtie2, for the DNA libraries of the bulk route. A prebuilt index is
    // used as given. Otherwise it is built from the host FASTA - or, when the
    // host came as a prebuilt HISAT2 index with no FASTA beside it, from the
    // sequences that index holds.
    //
    def ch_bowtie2_index = channel.empty()
    if (aligner != 'star') {
        def ch_need_bowtie2 = ch_aligners.filter { name -> name == 'bowtie2' }.first()
        if (bowtie2_index) {
            if (bowtie2_index.endsWith('.tar.gz') || bowtie2_index.endsWith('.tgz')) {
                UNTAR_BOWTIE2(channel.value([[id: 'bowtie2_index'], file(bowtie2_index, checkIfExists: true)]))
                ch_bowtie2_index = UNTAR_BOWTIE2.out.untar
            }
            else {
                ch_bowtie2_index = channel.value([[id: 'bowtie2_index'], file(bowtie2_index, checkIfExists: true)])
            }
        }
        else {
            def ch_bowtie2_source = ch_fasta
                .combine(ch_need_bowtie2)
                .map { meta, genome, _name -> [meta, genome] }
                .branch { _meta, genome ->
                    fasta: !(genome instanceof List)
                    none: true
                }
            HISAT2_INSPECT(
                ch_index.combine(ch_bowtie2_source.none).map { meta, index, _meta_fasta, _genome -> [meta, index] }
            )
            BOWTIE2_BUILD(ch_bowtie2_source.fasta.mix(HISAT2_INSPECT.out.fasta))
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
        }
    }

    emit:
    index = ch_index // channel: [ val(meta), path(index_dir) ]
    bowtie2_index = ch_bowtie2_index // channel: [ val(meta), path(index_dir) ] - empty unless a library needs it
    fasta = ch_fasta // channel: [ val(meta), path(fasta) ]
    gtf = ch_gtf // channel: [ val(meta), path(gtf) ]
}

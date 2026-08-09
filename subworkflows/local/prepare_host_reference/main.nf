//
// Make a HISAT2 index of the host genome available, building or downloading it
// only when the user has not supplied one.
//

include { NCBIGENOMEDOWNLOAD } from '../../../modules/nf-core/ncbigenomedownload/main'
include { GUNZIP             } from '../../../modules/nf-core/gunzip/main'
include { UNTAR              } from '../../../modules/nf-core/untar/main'
include { HISAT2_BUILD       } from '../../../modules/nf-core/hisat2/build/main'

workflow PREPARE_HOST_REFERENCE {

    take:
    fasta // string: path to a host genome FASTA (optionally gzipped), or null
    hisat2_index // string: path to a prebuilt HISAT2 index directory or tarball, or null
    host_accession // string: NCBI assembly accession (GCF_*/GCA_*), or null
    host_taxid // string: NCBI taxonomy ID, or null
    ncbi_group // string: ncbi-genome-download taxonomic group to search
    gtf // string: path to a GTF for splice-aware index building, or null

    main:

    def ch_gtf = gtf
        ? channel.value([[id: file(gtf).baseName], file(gtf, checkIfExists: true)])
        : channel.value([[:], []])

    def ch_index = channel.empty()
    def ch_fasta = channel.empty()

    if (hisat2_index) {
        //
        // A prebuilt index short-circuits everything else.
        //
        if (hisat2_index.endsWith('.tar.gz') || hisat2_index.endsWith('.tgz')) {
            UNTAR(channel.value([[id: 'hisat2_index'], file(hisat2_index, checkIfExists: true)]))
            ch_index = UNTAR.out.untar
        }
        else {
            ch_index = channel.value([[id: 'hisat2_index'], file(hisat2_index, checkIfExists: true)])
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

        HISAT2_BUILD(ch_fasta, ch_gtf, channel.value([[:], []]))
        ch_index = HISAT2_BUILD.out.index
    }

    emit:
    index = ch_index // channel: [ val(meta), path(index_dir) ]
    fasta = ch_fasta // channel: [ val(meta), path(fasta) ]
    gtf = ch_gtf // channel: [ val(meta), path(gtf) ]
}

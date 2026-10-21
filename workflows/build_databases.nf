/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BUILD_DATABASES: the classifiers' databases, built from your own genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// A run of its own, switched on by --build_databases. Prebuilt databases are
// where most analyses start, and they are what leaves the gaps: a parasite only
// VEuPathDB holds, a symbiont described last year, a vector-specific virus. Here
// the genomes of a sheet (--db_input: id, taxid, and a FASTA or an NCBI assembly
// accession) are built into the databases the analysis routes read, in exactly
// the layout their parameters expect:
//
//   kraken2      <outdir>/databases/kraken2      --kraken2_db (inspect.txt included)
//   bracken      the same directory              its k-mer distributions, one per length
//   krakenuniq   <outdir>/databases/krakenuniq   --krakenuniq_db
//   sylph        <outdir>/databases/sylph        --sylph_db, with --sylph_taxonomy the
//                                                 <db_name>.sylph_taxonomy.tsv in databases/
//
// Every sequence is assigned to its genome's taxid through its header, so no
// accession-to-taxid table is needed, and every taxid is checked against the
// taxonomy (--db_taxonomy, or NCBI's taxdump downloaded for the run) before
// anything is built. What cannot be built here, and why, is in docs/usage.md.
//

include { samplesheetToList                     } from 'plugin/nf-schema'
include { NCBIGENOMEDOWNLOAD as DB_GENOMEDOWNLOAD } from '../modules/nf-core/ncbigenomedownload/main'
include { UNTAR as DB_TAXDUMP                   } from '../modules/nf-core/untar/main'
include { DB_PREPARE_GENOME                     } from '../modules/local/databases/prepare/main'
include { DB_LINEAGE                            } from '../modules/local/databases/lineage/main'
include { DB_KRAKEN2_BUILD                      } from '../modules/local/databases/kraken2/main'
include { DB_BRACKEN_BUILD                      } from '../modules/local/databases/bracken/main'
include { DB_KRAKENUNIQ_BUILD                   } from '../modules/local/databases/krakenuniq/main'
include { DB_SYLPH_SKETCH                       } from '../modules/local/databases/sylph/main'
include { buildDatabaseTools                    } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'

workflow BUILD_DATABASES {

    main:

    def tools = buildDatabaseTools(params.build_databases)
    def meta = [id: params.db_name]

    def next_run = [
        tools.contains('kraken2') ? "--kraken2_db ${params.outdir}/databases/kraken2" : null,
        tools.contains('krakenuniq') ? "--krakenuniq_db ${params.outdir}/databases/krakenuniq" : null,
        tools.contains('sylph') ? "--run_sylph --sylph_db ${params.outdir}/databases/sylph/${params.db_name}.syldb --sylph_taxonomy ${params.outdir}/databases/${params.db_name}.sylph_taxonomy.tsv" : null,
    ].findAll { arg -> arg }.join(' ')
    log.info("Building ${tools.join(', ')} database(s) from ${params.db_input} into ${params.outdir}/databases/. Afterwards: ${next_run}")

    //
    // The genome sheet: one row per genome, a FASTA or an accession each.
    //
    def ch_rows = channel
        .fromList(samplesheetToList(params.db_input, "${projectDir}/assets/schema_databases.json"))
        .map { row_meta, fasta, accession ->
            if (!fasta && !accession) {
                error("--db_input: genome '${row_meta.id}' has neither a fasta nor an accession.")
            }
            if (fasta && accession) {
                error("--db_input: genome '${row_meta.id}' has both a fasta and an accession; give one.")
            }
            [row_meta + [taxid: row_meta.taxid.toString()], fasta, accession]
        }
        .branch { _meta, fasta, _accession ->
            local: fasta
            remote: true
        }

    //
    // MODULE: assemblies by accession. ncbi-genome-download reads its query
    // from a file, one per genome; meta and file travel as one multiMap so the
    // module's two inputs stay paired. GCA_* is GenBank, GCF_* RefSeq.
    //
    def ch_queries = ch_rows.remote
        .collectFile { row_meta, _fasta, accession -> ["${row_meta.id}.accession.txt", "${accession}\n"] }
        .map { query -> [query.name.replaceAll(/\.accession\.txt$/, ''), query] }
    def ch_download = ch_rows.remote
        .map { row_meta, _fasta, accession -> [row_meta.id, row_meta + [section: accession.startsWith('GCA_') ? 'genbank' : 'refseq']] }
        .join(ch_queries)
        .multiMap { _id, row_meta, query ->
            meta: row_meta
            query: query
        }
    DB_GENOMEDOWNLOAD(ch_download.meta, ch_download.query, [], 'all')
    def ch_downloaded = DB_GENOMEDOWNLOAD.out.fna.map { row_meta, fna ->
        def files = fna instanceof List ? fna : [fna]
        if (files.size() != 1) {
            error("ncbi-genome-download returned ${files.size()} files for '${row_meta.id}'; expected one assembly.")
        }
        [row_meta.findAll { key, _value -> key != 'section' }, files.first()]
    }

    //
    // MODULE: headers carrying the taxid, and a seqid-to-taxid map beside each.
    //
    DB_PREPARE_GENOME(
        ch_rows.local.map { row_meta, fasta, _accession -> [row_meta, fasta] }.mix(ch_downloaded)
    )
    def ch_genomes = DB_PREPARE_GENOME.out.genome

    //
    // The taxonomy: given, or NCBI's taxdump fetched for this run.
    //
    def ch_taxonomy = channel.empty()
    if (params.db_taxonomy) {
        ch_taxonomy = channel.value(file(params.db_taxonomy, checkIfExists: true))
    }
    else {
        DB_TAXDUMP(channel.value([[id: 'taxdump'], file(params.db_taxdump_url)]))
        ch_taxonomy = DB_TAXDUMP.out.untar.map { _meta, dir -> dir }
    }

    //
    // MODULE: every taxid checked against that taxonomy, and sylph-tax's
    // lineages written. Runs whatever is being built: it is the check.
    //
    def ch_genome_taxids = ch_genomes
        .map { row_meta, fna, _map -> "${fna.name}\t${row_meta.taxid}" }
        .collectFile(name: 'genome_taxids.tsv', newLine: true, sort: true)
        .map { table -> [meta, table] }
    DB_LINEAGE(ch_genome_taxids, ch_taxonomy)

    // Builds wait for the check: a database is not built around an unknown taxid.
    def ch_checked = DB_LINEAGE.out.sylph_taxonomy.map { _meta, table -> table }.first()
    def ch_fastas = ch_genomes
        .map { _meta, fna, _map -> fna }
        .collect(sort: true)
        .combine(ch_checked)
        .map { entries -> [meta, entries.init()] }

    if (tools.contains('kraken2')) {
        DB_KRAKEN2_BUILD(ch_fastas, ch_taxonomy)
        if (tools.contains('bracken')) {
            DB_BRACKEN_BUILD(DB_KRAKEN2_BUILD.out.db, params.db_bracken_read_lengths.toString().tokenize(',').collect { length -> length.trim() as Integer })
        }
    }

    if (tools.contains('krakenuniq')) {
        DB_KRAKENUNIQ_BUILD(
            ch_genomes
                .map { _meta, fna, map_file -> [fna, map_file] }
                .collect(sort: true)
                .combine(ch_checked)
                .map { entries -> [meta, entries.init()] },
            ch_taxonomy,
        )
    }

    if (tools.contains('sylph')) {
        DB_SYLPH_SKETCH(ch_fastas)
    }
}

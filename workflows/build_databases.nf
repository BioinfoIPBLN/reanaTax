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
//   metaphlan    <outdir>/databases/metaphlan    --metaphlan_db, with --metaphlan_index
//                                                 <db_name>
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
include { DB_MPA_GENES                          } from '../modules/local/databases/mpa_genes/main'
include { DB_MPA_CLUSTER                        } from '../modules/local/databases/mpa_cluster/main'
include { DB_MPA_CANDIDATES                     } from '../modules/local/databases/mpa_candidates/main'
include { DB_MPA_UNIQUENESS                     } from '../modules/local/databases/mpa_uniqueness/main'
include { DB_METAPHLAN_BUILD                    } from '../modules/local/databases/metaphlan/main'
include { buildDatabaseTools                    } from '../subworkflows/local/utils_nfcore_reanatax_pipeline'

workflow BUILD_DATABASES {

    main:

    def tools = buildDatabaseTools(params.build_databases)
    def meta = [id: params.db_name]

    def next_run = [
        tools.contains('kraken2') ? "--kraken2_db ${params.outdir}/databases/kraken2" : null,
        tools.contains('krakenuniq') ? "--krakenuniq_db ${params.outdir}/databases/krakenuniq" : null,
        tools.contains('sylph') ? "--run_sylph --sylph_db ${params.outdir}/databases/sylph/${params.db_name}.syldb --sylph_taxonomy ${params.outdir}/databases/${params.db_name}.sylph_taxonomy.tsv" : null,
        tools.contains('metaphlan') ? "--run_metaphlan --metaphlan_db ${params.outdir}/databases/metaphlan --metaphlan_index ${params.db_name}" : null,
    ].findAll { arg -> arg }.join(' ')
    log.info("Building ${tools.join(', ')} database(s) from ${params.db_input} into ${params.outdir}/databases/. Afterwards: ${next_run}")

    //
    // The genome sheet: one row per genome, a FASTA or an accession each, and
    // optionally its coding sequences (the MetaPhlAn build's genes).
    //
    def ch_rows = channel
        .fromList(samplesheetToList(params.db_input, "${projectDir}/assets/schema_databases.json"))
        .map { row_meta, fasta, accession, cds ->
            if (!fasta && !accession) {
                error("--db_input: genome '${row_meta.id}' has neither a fasta nor an accession.")
            }
            if (fasta && accession) {
                error("--db_input: genome '${row_meta.id}' has both a fasta and an accession; give one.")
            }
            [row_meta + [taxid: row_meta.taxid.toString()], fasta, accession, cds ?: []]
        }
        .branch { _meta, fasta, _accession, _cds ->
            local: fasta
            remote: true
        }

    //
    // MODULE: assemblies by accession. ncbi-genome-download reads its query
    // from a file, one per genome; meta and file travel as one multiMap so the
    // module's two inputs stay paired. GCA_* is GenBank, GCF_* RefSeq.
    //
    def ch_queries = ch_rows.remote
        .collectFile { row_meta, _fasta, accession, _cds -> ["${row_meta.id}.accession.txt", "${accession}\n"] }
        .map { query -> [query.name.replaceAll(/\.accession\.txt$/, ''), query] }
    def ch_download = ch_rows.remote
        .map { row_meta, _fasta, accession, _cds -> [row_meta.id, row_meta + [section: accession.startsWith('GCA_') ? 'genbank' : 'refseq']] }
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
        [row_meta.id, row_meta.findAll { key, _value -> key != 'section' }, files.first()]
    }
    // For a MetaPhlAn build the assembly's own CDS file is fetched too
    // (cds-fasta, in conf/modules.config), and used unless the sheet gives one.
    // An assembly NCBI has not annotated has none, and falls back to pyrodigal.
    def ch_remote = ch_downloaded
        .join(ch_rows.remote.map { row_meta, _fasta, _accession, cds -> [row_meta.id, cds] })
        .join(DB_GENOMEDOWNLOAD.out.cds.map { row_meta, cds -> [row_meta.id, cds instanceof List ? cds.first() : cds] }, remainder: true)
        .filter { _id, row_meta, _fasta, _sheet_cds, _ncbi_cds -> row_meta }
        .map { _id, row_meta, fasta, sheet_cds, ncbi_cds -> [row_meta, fasta, sheet_cds ?: (ncbi_cds ?: [])] }
    def ch_sources = ch_rows.local
        .map { row_meta, fasta, _accession, cds -> [row_meta, fasta, cds] }
        .mix(ch_remote)
        .map { row_meta, fasta, cds -> [row_meta + [has_cds: cds as boolean], fasta, cds] }

    //
    // MODULE: headers carrying the taxid, and a seqid-to-taxid map beside each.
    //
    DB_PREPARE_GENOME(ch_sources.map { row_meta, fasta, _cds -> [row_meta, fasta] })
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
        .map { row_meta, fna, _map -> "${fna.name}\t${row_meta.taxid}\t${row_meta.has_cds ? 1 : 0}" }
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

    //
    // MetaPhlAn: markers chosen the way MetaPhlAn 4's were (Blanco-Miguez et
    // al., Nat Biotechnol 2023), as HROM rebuilt the recipe for its own
    // catalogue: genes, gene families, each species' core families, and of
    // those the ones whose 150-nt fragments reach no other species' genome and
    // no --db_metaphlan_screen sequence. Genomes sharing a taxid are one
    // species. Gene calling waits for the lineage step, which refuses a
    // eukaryote without coding sequences and leaves viruses out.
    //
    if (tools.contains('metaphlan')) {
        def ch_mpa_taxonomy = DB_LINEAGE.out.mpa_taxonomy.map { _meta, table -> table }.first()
        DB_MPA_GENES(
            ch_genomes
                .map { row_meta, fna, _map -> [row_meta.id, row_meta, fna] }
                .join(ch_sources.map { row_meta, _fasta, cds -> [row_meta.id, cds] })
                .combine(ch_mpa_taxonomy)
                .map { _id, row_meta, fna, cds, _taxonomy -> [row_meta, fna, cds] }
        )
        DB_MPA_CLUSTER(
            DB_MPA_GENES.out.genes
                .map { _meta, _fna, faa, _table, _genome -> faa }
                .collect(sort: true)
                .map { proteins -> [meta, proteins] }
        )
        DB_MPA_CANDIDATES(
            DB_MPA_CLUSTER.out.clusters
                .combine(
                    DB_MPA_GENES.out.genes
                        .map { _meta, fna, _faa, table, genome -> [fna, table, genome] }
                        .collect(sort: true)
                        .map { files -> [files] }
                )
                .combine(ch_mpa_taxonomy)
        )
        def screen = params.db_metaphlan_screen
            ? params.db_metaphlan_screen.toString().tokenize(',').collect { path -> file(path.trim(), checkIfExists: true) }
            : []
        DB_MPA_UNIQUENESS(
            DB_MPA_CANDIDATES.out.candidates
                .map { db_meta, fna, _table, targets, _species -> [db_meta, fna, targets] }
                .combine(ch_genomes.map { _meta, fna, _map -> fna }.collect(sort: true).map { files -> [files] }),
            screen,
        )
        DB_METAPHLAN_BUILD(
            DB_MPA_CANDIDATES.out.candidates
                .join(DB_MPA_UNIQUENESS.out.hits)
                .combine(ch_mpa_taxonomy)
        )
    }
}

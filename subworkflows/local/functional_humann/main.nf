//
// Functional profiling of the non-host fraction with HUMAnN 3.
//
// This is the third layer Monteleone et al. (Microbiome 2026) describe: on top
// of "who is there" (Kraken2/Bracken), "what are they doing" - gene families and
// pathway abundances, optionally regrouped to KEGG orthologs.
//
// HUMAnN picks which pangenomes to align against from a taxonomic profile, and
// only reads MetaPhlAn's format. Rather than run MetaPhlAn as a second
// classifier - a second database to install, and a second set of abundances
// that would quietly disagree with the Bracken tables in the same report - the
// report we already have is translated with KrakenTools' kreport2mpa.py.
//
// HUMAnN keys on `s__` species lines and ignores everything above them, which
// is why the caller hands this Bracken's kreport rather than Kraken2's when
// Bracken ran: a read Kraken2 could only place at a genus contributes nothing
// until Bracken redistributes it to species. What no re-estimation can rescue
// is a clade Bracken itself leaves above species, so a profile dominated by
// higher-rank assignments still produces a thin functional table - a property
// of the classification, not a bug.
//
// Two databases have to agree for the nucleotide search to do anything. The
// profile names species the NCBI way, so --humann_nucleotide_db must be the
// species-named ChocoPhlAn (v201901_v31); against the SGB-named releases the
// lookup matches nothing and HUMAnN falls back to translated search alone. See
// modules/local/krakentools/kreport2mpa for the version line that encodes this.
//

include { KRAKENTOOLS_KREPORT2MPA } from '../../../modules/local/krakentools/kreport2mpa/main'
include { HUMANN3_HUMANN          } from '../../../modules/nf-core/humann3/humann/main'
include { HUMANN3_REGROUP         } from '../../../modules/nf-core/humann3/regroup/main'
include { HUMANN3_RENORM          } from '../../../modules/nf-core/humann3/renorm/main'

workflow FUNCTIONAL_HUMANN {

    take:
    ch_reads // channel: [ val(meta), [ path(fastq) ] ] - the non-host fraction
    ch_taxonomic_report // channel: [ val(meta), path(report) ] - Bracken's kreport, or Kraken2's
    nucleotide_db // string: ChocoPhlAn directory
    protein_db // string: UniRef (translated search) directory
    utility_db // string: HUMAnN utility mapping directory, or null
    regroup_groups // string: regrouping to apply (e.g. 'uniref90_ko'), or null
    renorm_units // string: 'cpm' or 'relab', or null to skip renormalisation

    main:

    def ch_versions = channel.empty()

    //
    // MODULE: Kraken-style report -> MetaPhlAn-style lineage profile
    //
    KRAKENTOOLS_KREPORT2MPA(ch_taxonomic_report)

    //
    // HUMAnN takes a single sequence file. Paired-end input is concatenated
    // rather than aligned as pairs: HUMAnN treats reads independently anyway,
    // so pairing carries no information it can use.
    //
    def ch_humann_input = ch_reads.map { meta, reads ->
        [meta, reads instanceof List ? reads : [reads]]
    }

    //
    // MODULE: HUMAnN
    //
    HUMANN3_HUMANN(
        ch_humann_input,
        KRAKENTOOLS_KREPORT2MPA.out.mpa,
        file(nucleotide_db, checkIfExists: true),
        file(protein_db, checkIfExists: true),
        utility_db ? file(utility_db, checkIfExists: true) : [],
    )

    //
    // MODULE: Regroup gene families onto another namespace. `uniref90_ko` is
    // what reproduces the KEGG-ortholog view used in the paper.
    //
    def ch_regrouped = channel.empty()

    if (regroup_groups) {
        HUMANN3_REGROUP(
            HUMANN3_HUMANN.out.genefamilies,
            regroup_groups,
            utility_db ? file(utility_db, checkIfExists: true) : [],
        )
        ch_regrouped = HUMANN3_REGROUP.out.regroup
    }

    //
    // MODULE: Renormalise. Raw HUMAnN abundances are in RPK, which is not
    // comparable between samples of different depth - anything cross-sample
    // wants CPM or relative abundance.
    //
    def ch_renormed = channel.empty()

    if (renorm_units) {
        // The units are set through ext.args in conf/modules.config; the
        // module itself takes only the table.
        HUMANN3_RENORM(regroup_groups ? ch_regrouped : HUMANN3_HUMANN.out.genefamilies)
        ch_renormed = HUMANN3_RENORM.out.renorm
    }

    emit:
    genefamilies = HUMANN3_HUMANN.out.genefamilies // channel: [ val(meta), path(tsv.gz) ]
    pathabundance = HUMANN3_HUMANN.out.pathabundance // channel: [ val(meta), path(tsv.gz) ]
    regrouped = ch_regrouped // channel: [ val(meta), path(tsv.gz) ]
    renormed = ch_renormed // channel: [ val(meta), path(tsv.gz) ]
    mpa = KRAKENTOOLS_KREPORT2MPA.out.mpa // channel: [ val(meta), path(txt) ]
    versions = ch_versions // channel: [ path(versions.yml) ]
}

// Convert a Kraken-style report into a MetaPhlAn-style lineage profile.
//
// This exists purely to feed HUMAnN, which selects the pangenomes to align
// against from a taxonomic profile and only speaks MetaPhlAn's format. Rather
// than run MetaPhlAn as a second classifier - a second database, a second set
// of abundances that would disagree with the Bracken tables in the same report
// - the report we already have is translated into that format. It is Bracken's
// `-w` kreport whenever Bracken ran; see subworkflows/local/functional_humann.
//
// The translation is lossy in one direction that matters: HUMAnN keys on
// `s__` species lines, so anything the classifier could only place above
// species is invisible to the functional step. See docs/usage.md.
//
// The version line is not decoration. HUMAnN 3.6.1's prescreen step
// (humann/search/prescreen.py, create_custom_database) scans the profile for a
// comment line containing one of two database version strings - "v3" or
// "vJan21" in config.py - and calls sys.exit() if it finds neither. A profile
// carrying only KrakenTools' own `#Classification` header aborts the run.
//
// Of the two accepted strings only the v3 one is useful here. HUMAnN builds
// its per-sample nucleotide database by lowercasing each `g__`/`s__` pair and
// matching it against the ChocoPhlAn filenames, which are species-named in
// v201901_v31 and SGB-named (t__SGB…) from vJan21 on. NCBI species names, which
// is all a Kraken2 or Bracken report can offer, match the former and can never
// match the latter - so declaring vJan21 would pass the version check and then
// select no pangenomes at all. Pair this with the v201901_v31 ChocoPhlAn.
process KRAKENTOOLS_KREPORT2MPA {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/krakentools:1.2.1--pyh7e72e81_0'
        : 'biocontainers/krakentools:1.2.1--pyh7e72e81_0'}"

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.mpa.txt"), emit: mpa
    tuple val("${task.process}"), val('krakentools'), eval("echo 1.2.1"), emit: versions_krakentools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpa_version = task.ext.mpa_version ?: 'mpa_v30_CHOCOPhlAn_201901'
    // kreport2mpa.py cannot be told to write a version line, so it writes the
    // body to tmp/ and the header is prepended on the way out. tmp/, not the
    // task directory: an intermediate sitting next to the real output is one
    // loose glob away from being published or collected alongside it.
    """
    mkdir -p tmp

    kreport2mpa.py \\
        ${args} \\
        --report ${report} \\
        --output tmp/${prefix}.mpa.txt

    printf '#%s\\n' '${mpa_version}' > ${prefix}.mpa.txt
    cat tmp/${prefix}.mpa.txt >> ${prefix}.mpa.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpa_version = task.ext.mpa_version ?: 'mpa_v30_CHOCOPhlAn_201901'
    """
    printf '#%s\\n' '${mpa_version}' > ${prefix}.mpa.txt
    """
}

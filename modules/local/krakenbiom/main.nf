//
// A BIOM table of the classification, for tools that speak BIOM and not TSV.
//
// Entirely for interoperability. Nothing downstream in this pipeline reads it:
// the combined tables are already a taxa-by-sample matrix, which is why the
// differential-abundance step skips the kraken-biom round trip that DiffShot
// does. This exists so a run's output can be handed to QIIME 2, phyloseq,
// microbiome or anything else that wants BIOM without the reader having to
// reconstruct a taxonomy hierarchy from a TSV.
//
// Two properties of kraken-biom decide almost everything in the script below.
//
// It names each sample after the FILE, taking the basename minus one extension,
// so `SRR123.kraken2.report_bracken.txt` would become a sample called
// `SRR123.kraken2.report_bracken`. The reports are therefore re-linked under
// `<sample id>.txt` first, matched by name rather than by position so a sorting
// difference between Groovy and the shell cannot mislabel a column.
//
// And it parses the report with a fixed six-field reader. With
// --kraken2_report_minimizer_data the reports have eight columns, and the extra
// two would silently shift rank, taxid and name by two positions - a table full
// of plausible nonsense rather than an error. The awk drops them back to six.
//
process KRAKENBIOM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/kraken-biom:1.2.0--pyh5e36f6f_0'
        : 'quay.io/biocontainers/kraken-biom:1.2.0--pyh5e36f6f_0'}"

    input:
    tuple val(meta), path(reports, stageAs: 'reports/*')
    path keep_table
    path metadata

    output:
    tuple val(meta), path("*.biom")           , emit: biom
    tuple val(meta), path("*.unfiltered.biom"), emit: unfiltered, optional: true
    tuple val("${task.process}"), val('kraken-biom'), eval("kraken-biom --version 2>&1 | sed -n '1s/.*[[:space:]]//p'"), emit: versions_krakenbiom, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--max D --min S'
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def keep_arg = keep_table ? "--keep-table ${keep_table}" : ''
    def metadata_arg = metadata ? "--metadata ${metadata}" : ''
    def unfiltered_arg = keep_table ? "--unfiltered ${prefix}.unfiltered.biom" : ''
    """
    mkdir -p named

    # Matched by name rather than by position: `set -e -o pipefail` makes a
    # missing file an unexplained abort, so the loop looks and then says so.
    for name in ${meta.names}; do
        found=''
        for candidate in reports/"\${name}".*; do
            if [ -e "\${candidate}" ]; then
                found="\${candidate}"
                break
            fi
        done
        if [ -z "\${found}" ]; then
            echo "KRAKENBIOM: no report staged for sample \${name}" >&2
            exit 1
        fi
        awk -F'\\t' -v OFS='\\t' 'NF >= 8 { print \$1, \$2, \$3, \$6, \$7, \$8; next } { print }' \\
            "\${found}" > "named/\${name}.txt"
    done

    kraken-biom \\
        named/*.txt \\
        --fmt json \\
        -o raw.biom \\
        ${args}

    biom_export.py \\
        --biom raw.biom \\
        --output ${prefix}.biom \\
        ${keep_arg} \\
        ${metadata_arg} \\
        ${unfiltered_arg} \\
        ${args2}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.biom
    """
}

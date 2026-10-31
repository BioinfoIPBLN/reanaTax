// Copies a genome or index this run built into --reference_cache.
//
// The copy is made beside the entry and renamed into place, so another run
// sees an entry whole or not at all. When two runs build the same entry at
// once, the second rename fails and that copy is discarded; the entry the
// first one wrote is the same build.
//
// Runs on the host, without a container: it only copies files, and the cache
// directory is outside every path a container would be given.
process REFERENCE_CACHE_STORE {
    tag "${meta.id}:${entry.tokenize('/')[-2]}"
    label 'process_single'

    input:
    tuple val(meta), path(items, stageAs: 'store/*'), val(entry), val(subdir), val(key)

    output:
    tuple val(meta), val(entry), emit: entry

    when:
    task.ext.when == null || task.ext.when

    script:
    def key_lines = key.readLines().collect { line -> "'" + line.replace("'", "'\\''") + "'" }.join(' ')
    """
    if [ ! -d '${entry}' ]; then
        mkdir -p "\$(dirname '${entry}')"
        tmp='${entry}.tmp.'\$\$
        mkdir -p "\$tmp/${subdir}"
        cp -rL store/* "\$tmp/${subdir}/"
        printf '%s\\n' ${key_lines} > "\$tmp/key.txt"
        mv -T "\$tmp" '${entry}' 2> /dev/null || rm -rf "\$tmp"
    fi
    """
}

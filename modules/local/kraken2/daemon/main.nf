// Kraken2 classification through the `k2` wrapper's classifier daemon.
//
// `kraken2` (the legacy Perl wrapper) loads the whole hash into the process on
// every invocation. Against a 340 GB core_nt on network storage that is ~30
// minutes of pure I/O per sample, repeated for every sample in the cohort.
// `k2 classify --use-daemon` loads the index once into a background process
// that outlives the task; later invocations hand their arguments to it and
// start classifying immediately.
//
// Three properties of that daemon dictate how this process is written.
//
//   1. It is addressed through HARDCODED paths in /tmp - `classify.pid` plus
//      the `classify_stdin`/`classify_stdout` FIFOs. Singularity binds the
//      host /tmp into the container and shares the host PID namespace, so the
//      daemon a task starts is reachable from the next task. It also means the
//      daemon is shared with anything else on the node using k2, and that it
//      is NOT confined to this pipeline run.
//
//   2. The control FIFO carries one conversation at a time and takes no lock,
//      so two concurrent clients would interleave their handshakes. That costs
//      nothing to prevent: the point of the daemon is that the index is loaded
//      once, after which each sample classifies at full speed, and the
//      alternative was paying the load repeatedly anyway.
//
//      `maxForks 1` is NOT sufficient to prevent it. Nextflow counts maxForks
//      per process, and this module is included twice - once as KRAKEN2_DAEMON
//      and once aliased as KRAKEN2_DAEMON_SHUFFLED - which makes two processes
//      with one slot each, not one process with one slot. The shuffled-read
//      control therefore ran concurrently with the real classification, two
//      clients on one FIFO, and the observed result was a client blocked in
//      open(2) forever plus a SECOND daemon loading its own private copy of the
//      index: 2 x 314 GB resident for a 344 GB database.
//
//      The mutual exclusion is therefore taken where it is actually needed, on
//      a lock file the two process instances share. The kernel drops an flock
//      when the holder dies, so a killed task cannot leave the cohort wedged -
//      which a maxForks counter, or a lock made out of mkdir, would.
//
//   3. The daemon keys its resident index by the `--db` argument it is handed.
//      The database is therefore taken as a value and passed as the absolute
//      path the user gave, NOT staged into the work directory: a staged symlink
//      resolves to a different path in every task, which would make the daemon
//      treat each sample as a new database and load another copy alongside the
//      first.
//
//      That has a consequence worth stating, because getting it wrong looks
//      like a corrupt database rather than a missing mount: `autoMounts` only
//      binds what Nextflow knows about - the work directory and staged inputs.
//      An unstaged path is invisible to it, so the container would not see the
//      database at all and `k2` would report "is not a valid database" against
//      a directory that is perfectly intact on the host. Hence the explicit
//      bind in containerOptions below.
//
// The daemon does not stop when the run ends. It holds the index until
// `k2 clean --stop-daemon` is run on that node - see docs/usage.md.
//
// A FOURTH property is what the timeout and the retry-time reset below exist
// for. Nothing in the handshake has a deadline. A client that opens the
// control FIFO when no daemon is reading it blocks in open(2) - indefinitely,
// at 0% CPU, with an empty log and no error - and `maxForks 1` means the rest
// of the cohort queues behind it forever. `k2` decides whether a daemon is
// already running from /tmp/classify.pid, and it only asks whether that PID
// exists, not whether it is a k2 process; the PID it records is not always the
// daemon's (the wrapper logs an empty "Started background classifier process
// with PID:" and has been seen to leave a PID belonging to a kernel thread).
// Once that happens the check can never fail again and every later task hangs.
//
// So the wait is bounded rather than diagnosed: a timeout turns the hang into
// an ordinary failure, and only on the retry - once THIS run has established
// that the daemon is unreachable - is the /tmp state cleared so `k2` starts a
// fresh one. Clearing it up front would be wrong: the daemon is shared with
// anything else on the node, and a healthy daemon serving another run must not
// be torn down on suspicion.
//
// The reset VERIFIES the kill rather than assuming it, because a wedged daemon
// ignores SIGTERM and clearing /tmp around a survivor is worse than the hang.
// See the reset_daemon block for the mechanism.
process KRAKEN2_DAEMON {
    tag "$meta.id"
    label 'process_high'
    // See note 2: the control FIFO is not concurrency-safe.
    maxForks 1

    // See note 3: the database is passed by path, not staged, so it has to be
    // bound into the container by hand. Read-only under Docker because nothing
    // here writes to the index - though note that Docker gives every task its
    // own /tmp and PID namespace, so the daemon cannot outlive a task there and
    // --kraken2_use_daemon is really an Apptainer/Singularity feature.
    containerOptions {
        if (workflow.containerEngine in ['singularity', 'apptainer']) {
            "-B ${db}"
        }
        else if (workflow.containerEngine == 'docker') {
            "-v ${db}:${db}:ro"
        }
        else {
            ''
        }
    }

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0f/0f827dcea51be6b5c32255167caa2dfb65607caecdc8b067abd6b71c267e2e82/data' :
        'community.wave.seqera.io/library/kraken2_coreutils_pigz:920ecc6b96e2ba71' }"

    input:
    tuple val(meta), path(reads)
    val db
    val save_output_fastqs
    val save_reads_assignment

    output:
    tuple val(meta), path('*.classified{.,_}*'), optional: true, emit: classified_reads_fastq
    tuple val(meta), path('*.unclassified{.,_}*'), optional: true, emit: unclassified_reads_fastq
    tuple val(meta), path('*classifiedreads.txt'), optional: true, emit: classified_reads_assignment
    tuple val(meta), path('*report.txt'), emit: report
    tuple val(meta), path('*.k2.log'), emit: log
    tuple val("${task.process}"), val('kraken2'), eval('kraken2 --version 2>&1 | head -1 | sed "s/^.*Kraken version //; s/ .*//"'), topic: versions, emit: versions_kraken2
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1 | sed "s/pigz //g"'), topic: versions, emit: versions_pigz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def classify_timeout = task.ext.timeout ?: params.kraken2_daemon_timeout
    // See note 4. Only after an attempt of ours has already timed out, and only
    // inside the lock: holding it is what guarantees no other task of this run
    // is mid-conversation with the daemon we are about to stop.
    def reset_daemon = task.attempt > 1 ? """
        echo "[KRAKEN2_DAEMON] attempt ${task.attempt}: the previous attempt did not return within ${classify_timeout}, so the daemon is unreachable. Stopping it and clearing /tmp so a fresh one is started." >&2
        # Read the PID BEFORE k2 clean, which removes the file. k2 clean signals
        # with SIGTERM, and a daemon blocked in pipe_read does not act on it -
        # observed on the leprosy cohort, where one sat in that state for 18 hours
        # and only SIGKILL cleared it. Removing the pidfile around a survivor
        # ORPHANS it: still alive, still holding the hardcoded /tmp FIFOs, but no
        # longer named anywhere, so the next daemon shares those paths with it and
        # the two race for client messages. That is worse than the hang this reset
        # exists to fix, and it is why the kill is verified rather than assumed.
        stuck=\$(tr -dc '0-9' < /tmp/classify.pid 2>/dev/null || true)
        timeout 60 k2 clean --stop-daemon >/dev/null 2>&1 || true
        # Same safety rule as the pipeline's own cleanup: only ever signal a PID
        # that is genuinely a classify process and genuinely ours. k2 records a PID
        # it has not verified - a kernel thread's, in one observed case - so an
        # unchecked kill here could signal an unrelated process.
        if [ -n "\${stuck}" ] \\
           && [ "\$(cat /proc/\${stuck}/comm 2>/dev/null)" = "classify" ] \\
           && [ "\$(stat -c%u /proc/\${stuck} 2>/dev/null)" = "\$(id -u)" ]; then
            kill -TERM "\${stuck}" 2>/dev/null || true
            for _ in 1 2 3 4 5 6 7 8 9 10; do
                kill -0 "\${stuck}" 2>/dev/null || break
                sleep 1
            done
            if kill -0 "\${stuck}" 2>/dev/null; then
                echo "[KRAKEN2_DAEMON] PID \${stuck} ignored SIGTERM; sending SIGKILL." >&2
                kill -KILL "\${stuck}" 2>/dev/null || true
                sleep 2
            fi
            if kill -0 "\${stuck}" 2>/dev/null; then
                echo "[KRAKEN2_DAEMON] PID \${stuck} survived SIGKILL. NOT clearing /tmp: a fresh daemon would share its FIFOs and the two would race. Stop it by hand before resuming." >&2
                exit 1
            fi
        fi
        rm -f /tmp/classify.pid /tmp/classify_stdin /tmp/classify_stdout
""" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "" : "--paired"
    def classified = meta.single_end ? "${prefix}.classified.fastq" : "${prefix}.classified#.fastq"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq" : "${prefix}.unclassified#.fastq"
    // Every path handed to the daemon is absolute, and that is load-bearing.
    // The daemon is a process that is already running: its working directory is
    // the work dir of whichever task started it, not of the task sending the
    // command. With relative paths it writes sample N's report into sample 1's
    // directory, Nextflow looks in the right place, finds nothing, and fails the
    // task on a missing output. Observed exactly that: one task dir holding four
    // different samples' reports. This is invisible until the daemon is actually
    // reused, which is why it only surfaced once the PID namespace was shared.
    def classified_option = save_output_fastqs ? "--classified-out \$PWD/${classified}" : ""
    def unclassified_option = save_output_fastqs ? "--unclassified-out \$PWD/${unclassified}" : ""
    def readclassification_option = save_reads_assignment ? "--output \$PWD/${prefix}.kraken2.classifiedreads.txt" : "--output /dev/null"
    // Inputs for the same reason. The reports came out correct without this, so
    // k2 evidently resolves reads client-side, but relying on that distinction
    // is not worth a silently mis-paired classification.
    def reads_abs = [reads].flatten().collect { "\$PWD/${it}" }.join(' ')
    // --compression_level. These FASTQs are only written with
    // --kraken2_save_reads, and they are the size of the non-host library.
    def level = task.ext.compression != null ? "-${task.ext.compression} " : ''
    def compress_reads_command = save_output_fastqs ? "pigz ${level}-p ${task.cpus} *.fastq" : ""

    // No --gzip-compressed: `k2 classify` detects gzip/bz2/xz from the file
    // itself and rejects the flag the legacy wrapper needs.
    """
    # See note 2. One conversation with the daemon at a time, across every
    # process instance of this module - which is what maxForks cannot do. The
    # wait is deliberately unbounded: whoever holds the lock is classifying,
    # and the kernel releases it if they die. Only the classify call itself is
    # under the timeout, so a task queued behind a slow one is not charged for
    # the wait. The lock is per-user because the /tmp paths k2 uses are.
    # The redirect is `9>>` rather than `9>`: Nextflow runs the script under
    # `bash -C`, so a truncating redirect onto a lock file that already exists
    # is refused outright. Appending never writes a byte and never truncates.
    (
        flock 9
${reset_daemon}
        # Recorded, not acted on: if this task does time out, the log says
        # whether the PID k2 was trusting belonged to a k2 process at all.
        if [ -f /tmp/classify.pid ]; then
            daemon_pid=\$(tr -dc '0-9' < /tmp/classify.pid 2>/dev/null || true)
            daemon_cmd=\$(tr '\\0' ' ' < /proc/\${daemon_pid:-0}/cmdline 2>/dev/null || true)
            echo "[KRAKEN2_DAEMON] /tmp/classify.pid names PID \${daemon_pid:-<empty>}; its cmdline is \${daemon_cmd:-<no such process>}" >&2
        fi

        timeout --preserve-status --signal=TERM --kill-after=60 ${classify_timeout} \\
        k2 classify \\
            --use-daemon \\
            --db ${db} \\
            --threads ${task.cpus} \\
            --report \$PWD/${prefix}.kraken2.report.txt \\
            --log ${prefix}.k2.log \\
            ${unclassified_option} \\
            ${classified_option} \\
            ${readclassification_option} \\
            ${paired} \\
            ${args} \\
            ${reads_abs}
    ) 9>>/tmp/reanatax_k2_daemon.\$(id -u).lock

    ${compress_reads_command}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def classified = meta.single_end ? "${prefix}.classified.fastq.gz" : "${prefix}.classified_1.fastq.gz ${prefix}.classified_2.fastq.gz"
    def unclassified = meta.single_end ? "${prefix}.unclassified.fastq.gz" : "${prefix}.unclassified_1.fastq.gz ${prefix}.unclassified_2.fastq.gz"
    """
    touch ${prefix}.kraken2.report.txt
    touch ${prefix}.k2.log
    if [ "${save_output_fastqs}" == "true" ]; then
        touch ${classified}
        touch ${unclassified}
    fi
    if [ "${save_reads_assignment}" == "true" ]; then
        touch ${prefix}.kraken2.classifiedreads.txt
    fi
    """
}

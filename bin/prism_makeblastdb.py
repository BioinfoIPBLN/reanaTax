#!/usr/bin/env python3
"""
prism_makeblastdb.py -- PRISM's makeblastdb, building a database as several
in parallel.

makeblastdb is single-threaded: PRISM's per-library database (up to 115 GB of
FASTA on a mock community) took 4 to 15 min on one core while the task's other
cores idled. Here the FASTA is cut into consecutive byte ranges at record
boundaries, each range is streamed into a makeblastdb of its own, several at a
time, and an alias file lists the resulting databases in order. BLAST reads
the alias as one database whose sequences are in the FASTA's order, with the
same letters and sequence count, so scores and e-values are unchanged.

The hits are not byte-identical to a single database's: when more subjects
score equally than -max_target_seqs keeps, which of them are kept depends on
how the database is cut into volumes, as it already depends on blastn's thread
count. See assets/prism/README.md for what that changed on a mock library.

modules/local/prism/run gives PRISM a makeblastdb that runs

  prism_makeblastdb.py <real makeblastdb> <jobs> <makeblastdb arguments>

A call without -in and -out, from standard input, with jobs < 2, or on a FASTA
under 1 GB is passed to makeblastdb unchanged. Ranges are at most 4 GB, so
each database is one volume and the alias file is flat. Both sizes can be set
in bytes with PRISM_MAKEBLASTDB_MIN_BYTES and PRISM_MAKEBLASTDB_CHUNK_BYTES.
"""
import math
import os
import subprocess
import sys
import threading
import time

MIN_SIZE = int(os.environ.get("PRISM_MAKEBLASTDB_MIN_BYTES", 1 << 30))
CHUNK_SIZE = int(os.environ.get("PRISM_MAKEBLASTDB_CHUNK_BYTES", 4 << 30))
BLOCK = 16 << 20


def option(args, flag):
    return args[args.index(flag) + 1] if flag in args else None


def replace(args, flag, value):
    args = list(args)
    if flag in args:
        args[args.index(flag) + 1] = value
    else:
        args += [flag, value]
    return args


def boundaries(fasta, chunks):
    """Byte offsets where each of up to `chunks` ranges starts: 0, then the start of
    the first record at or after each equal share of the file."""
    size = os.path.getsize(fasta)
    starts = [0]
    with open(fasta, "rb") as handle:
        for k in range(1, chunks):
            handle.seek(max(size * k // chunks - 1, 0))
            handle.readline()
            while True:
                at = handle.tell()
                line = handle.readline()
                if not line or line.startswith(b">"):
                    break
            if line and at > starts[-1]:
                starts.append(at)
    return starts + [size]


def feed(fasta, start, end, sink):
    try:
        with open(fasta, "rb") as handle:
            handle.seek(start)
            left = end - start
            while left > 0:
                block = handle.read(min(BLOCK, left))
                if not block:
                    break
                sink.write(block)
                left -= len(block)
    except BrokenPipeError:
        pass
    finally:
        try:
            sink.close()
        except BrokenPipeError:
            pass


def main():
    makeblastdb, jobs, args = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    fasta, out = option(args, "-in"), option(args, "-out")
    if jobs < 2 or fasta in (None, "-") or out is None or not os.path.isfile(fasta) or os.path.getsize(fasta) < MIN_SIZE:
        os.execv(makeblastdb, [makeblastdb] + args)
    starts = boundaries(fasta, max(jobs, math.ceil(os.path.getsize(fasta) / CHUNK_SIZE)))
    ranges = list(zip(starts[:-1], starts[1:]))
    names = [f"{os.path.basename(out)}.p{k:03d}" for k in range(len(ranges))]
    title = option(args, "-title") or fasta

    pending = list(range(len(ranges)))
    running = {}
    failed = []
    while pending or running:
        while pending and len(running) < jobs:
            k = pending.pop(0)
            argv = replace(replace(replace(args, "-in", "-"), "-out", os.path.join(os.path.dirname(out), names[k])), "-title", title)
            job = subprocess.Popen([makeblastdb] + argv, stdin=subprocess.PIPE)
            feeder = threading.Thread(target=feed, args=(fasta, *ranges[k], job.stdin), daemon=True)
            feeder.start()
            running[k] = (job, feeder)
        done = [k for k, (job, _) in running.items() if job.poll() is not None]
        for k in done:
            job, feeder = running.pop(k)
            feeder.join()
            if job.returncode != 0:
                failed.append(k)
        if not done:
            time.sleep(0.2)
    if failed:
        print(f"[prism_makeblastdb] makeblastdb failed on {len(failed)} of {len(ranges)} range(s)", file=sys.stderr)
        sys.exit(1)

    suffix = ".pal" if option(args, "-dbtype") == "prot" else ".nal"
    with open(out + suffix, "w") as alias:
        alias.write(f"TITLE {title}\nDBLIST {' '.join(names)}\n")
    print(f"[prism_makeblastdb] {len(ranges)} range(s), {jobs} at a time, alias {os.path.basename(out + suffix)}", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
prism_blastn.py -- PRISM's blastn, run as several processes over slices of the query.

PRISM's megablast keeps about half of its threads busy: at -num_threads 16 it
used 8.5 cores on a 9.5M-read library, and -mt_mode 1 was slower still. Here
the query is cut into contiguous slices, each searched by its own blastn with
an equal share of -num_threads, and the outputs are joined in query order.

The output is not byte-identical to one blastn's, and neither is one blastn's
at another -num_threads: when more than ten subjects score equally, which ones
fill -max_target_seqs 10 depends on the threads and the query batch. On 200,000
reads of ERR2935854 (16 threads vs 4 slices of 4) the best bitscore, the
species hit, the species at the best score and whether a read is unique to one
species were the same for every read; 77% of reads listed other accessions
among their equal hits.

modules/local/prism/run gives PRISM a blastn that runs

  prism_blastn.py <real blastn> <slices> <blastn arguments>

A call without -query and -out, or with fewer reads than slices, is passed to
blastn unchanged.
"""
import os
import shutil
import subprocess
import sys
import tempfile


def replace(args, flag, value):
    args = list(args)
    if flag in args:
        args[args.index(flag) + 1] = value
    else:
        args += [flag, value]
    return args


def main():
    blastn, slices, args = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    query = args[args.index("-query") + 1] if "-query" in args else None
    out = args[args.index("-out") + 1] if "-out" in args else None
    threads = int(args[args.index("-num_threads") + 1]) if "-num_threads" in args else 1
    if slices < 2 or query is None or out is None:
        os.execv(blastn, [blastn] + args)
    with open(query) as handle:
        records = sum(1 for line in handle if line.startswith(">"))
    slices = min(slices, records)
    if slices < 2:
        os.execv(blastn, [blastn] + args)

    work = tempfile.mkdtemp(prefix=os.path.basename(out) + ".slices.", dir=os.path.dirname(os.path.abspath(out)))
    try:
        queries = [os.path.join(work, f"query{k}.fa") for k in range(slices)]
        outputs = [os.path.join(work, f"out{k}") for k in range(slices)]
        sinks = [open(path, "w") for path in queries]
        with open(query) as handle:
            record, k = -1, 0
            for line in handle:
                if line.startswith(">"):
                    record += 1
                    k = record * slices // records
                sinks[k].write(line)
        for sink in sinks:
            sink.close()

        share = str(max(1, threads // slices))
        jobs = []
        for k in range(slices):
            argv = replace(replace(replace(args, "-query", queries[k]), "-out", outputs[k]), "-num_threads", share)
            jobs.append(subprocess.Popen([blastn] + argv))
        codes = [job.wait() for job in jobs]
        failed = [code for code in codes if code != 0]
        if failed:
            print(f"[prism_blastn] {len(failed)} of {slices} slice(s) failed", file=sys.stderr)
            sys.exit(failed[0] if failed[0] > 0 else 1)

        with open(out, "wb") as sink:
            for path in outputs:
                with open(path, "rb") as part:
                    shutil.copyfileobj(part, sink)
        print(f"[prism_blastn] {records} read(s) in {slices} slices of {share} thread(s)", file=sys.stderr)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
prism_blastdbcmd.py -- PRISM's blastdbcmd, extracting from the volumes of a
database in parallel.

PRISM builds a BLAST database per library from the sequences of a few thousand
taxids, with one `blastdbcmd -taxidlist` over the whole database: on core_nt
(92 volumes, 300 GB) that took 10 min on one core for 115 GB of FASTA. Here
the same call runs once per volume, several at a time, and the outputs are
joined in the order the database's alias file lists its volumes, which is the
order one blastdbcmd walks them: the FASTA is byte-identical (checked on that
115 GB extraction: 363 s with 8 processes against 531-786 s for one call, on a
busy node).

modules/local/prism/run gives PRISM a blastdbcmd that runs

  prism_blastdbcmd.py <real blastdbcmd> <jobs> <blastdbcmd arguments>

A call that is not a -taxidlist/-taxids extraction in %f format to -out from a
database with an alias file, or with jobs < 2, is passed to blastdbcmd
unchanged. A volume with no sequence of the taxids is empty; when no volume has
one, the whole call fails as blastdbcmd's would. Each volume's part is deleted
once appended, so the disk holds about one copy of the FASTA, not two.
"""
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time


def option(args, flag):
    return args[args.index(flag) + 1] if flag in args else None


def replace(args, flag, value):
    args = list(args)
    args[args.index(flag) + 1] = value
    return args


def alias_file(db, suffix):
    candidates = [db] + [os.path.join(d, db) for d in os.environ.get("BLASTDB", "").split(":") if d]
    for path in candidates:
        if os.path.isfile(path + suffix):
            return path + suffix
    return None


def volumes(alias):
    with open(alias) as handle:
        for line in handle:
            if line.startswith("DBLIST"):
                names = shlex.split(line[len("DBLIST"):])
                base = os.path.dirname(alias)
                return [name if os.path.isabs(name) else os.path.join(base, name) for name in names]
    return []


def main():
    blastdbcmd, jobs, args = sys.argv[1], int(sys.argv[2]), sys.argv[3:]
    db, out = option(args, "-db"), option(args, "-out")
    by_taxid = "-taxidlist" in args or "-taxids" in args
    plain = option(args, "-outfmt") in (None, "%f") and "-entry" not in args and "-entry_batch" not in args
    suffix = ".pal" if option(args, "-dbtype") == "prot" else ".nal"
    alias = alias_file(db, suffix) if db else None
    if jobs < 2 or out is None or not by_taxid or not plain or alias is None:
        os.execv(blastdbcmd, [blastdbcmd] + args)
    parts = volumes(alias)
    if len(parts) < 2:
        os.execv(blastdbcmd, [blastdbcmd] + args)

    work = tempfile.mkdtemp(prefix=os.path.basename(out) + ".volumes.", dir=os.path.dirname(os.path.abspath(out)))
    try:
        outputs = [os.path.join(work, f"volume{k}") for k in range(len(parts))]
        errors = [os.path.join(work, f"volume{k}.err") for k in range(len(parts))]
        pending = list(range(len(parts)))
        running = {}
        codes = {}
        while pending or running:
            while pending and len(running) < jobs:
                k = pending.pop(0)
                argv = replace(replace(args, "-db", parts[k]), "-out", outputs[k])
                running[k] = subprocess.Popen([blastdbcmd] + argv, stderr=open(errors[k], "w"))
            done = [k for k, job in running.items() if job.poll() is not None]
            for k in done:
                codes[k] = running.pop(k).returncode
            if not done:
                time.sleep(0.2)
        found = 0
        for k in range(len(parts)):
            with open(errors[k]) as handle:
                message = handle.read()
            if codes[k] == 0:
                found += 1
            elif "not found" in message and not (os.path.exists(outputs[k]) and os.path.getsize(outputs[k])):
                open(outputs[k], "w").close()
            else:
                sys.stderr.write(message)
                print(f"[prism_blastdbcmd] blastdbcmd failed on {os.path.basename(parts[k])}", file=sys.stderr)
                sys.exit(codes[k] if codes[k] > 0 else 1)
        if found == 0:
            print(f"Error: [blastdbcmd] Taxonomy ID(s) not found in the {os.path.basename(db)} database.", file=sys.stderr)
            sys.exit(1)
        with open(out, "wb") as sink:
            for path in outputs:
                with open(path, "rb") as part:
                    shutil.copyfileobj(part, sink, 16 << 20)
                os.remove(path)
        print(f"[prism_blastdbcmd] {len(parts)} volumes, {found} with sequences, {jobs} at a time", file=sys.stderr)
    finally:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()

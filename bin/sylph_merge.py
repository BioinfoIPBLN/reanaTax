#!/usr/bin/env python3
"""
sylph_merge.py -- per-sample sylph and sylph-tax outputs into cohort tables.

Two kinds of input, two kinds of output.

  <prefix>.genomes.tsv
      Every sample's sylph rows under one header: one row per genome called,
      with its adjusted ANI, taxonomic and sequence abundance, and effective
      coverage. Sample_file already holds the sample id (SYLPH_PROFILE rewrites
      it), so the table needs no extra key.

  <prefix>.relative_abundance.tsv, <prefix>.sequence_abundance.tsv
      With sylph-tax profiles only: clades by samples, in the layout
      `sylph-tax merge` writes. relative_abundance is coverage-normalised, the
      same quantity MetaPhlAn reports; sequence_abundance is the share of reads
      assigned, the same quantity a Kraken2 percentage is.

`sylph-tax merge` itself is not used because it cannot tell an empty profile
from a missing one. A sample in which sylph detected nothing must stay in the
table as a column of zeros; dropping it would make an empty library look like
one that was never profiled.
"""
import argparse
import sys
from pathlib import Path

ABUNDANCES = ("relative_abundance", "sequence_abundance")


def sample_id(path, suffix):
    name = Path(path).name
    return name[: -len(suffix)] if name.endswith(suffix) else Path(path).stem


def stack_genomes(paths, out):
    header = None
    rows = 0
    with open(out, "w", encoding="utf-8") as sink:
        for path in sorted(paths, key=lambda p: sample_id(p, ".sylph.tsv")):
            with open(path, encoding="utf-8") as handle:
                first = handle.readline()
                if not first.strip():
                    continue
                if header is None:
                    header = first
                    sink.write(first)
                elif first != header:
                    sys.exit(f"ERROR: {path}: header differs from the first profile's. "
                             "Were these written by different sylph versions?")
                for line in handle:
                    if line.strip():
                        sink.write(line)
                        rows += 1
    return rows


def read_taxprof(path):
    header = None
    clades = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                continue
            clades[fields[0]] = dict(zip(header, fields))
    if clades:
        missing = [c for c in ABUNDANCES if c not in header]
        if missing:
            sys.exit(f"ERROR: {path}: no {', '.join(missing)} column. "
                     f"Header was: {chr(9).join(header)}")
    return clades


def merge_taxprof(paths, prefix):
    samples = {}
    for path in paths:
        sid = sample_id(path, ".sylphmpa")
        if sid in samples:
            sys.exit(f"ERROR: two sylph-tax profiles for sample {sid}.")
        samples[sid] = read_taxprof(path)
    order = sorted(samples)
    # Lexicographic order keeps every clade below its parent, since a child's
    # name is its parent's with `|<rank>__<name>` appended.
    clades = sorted({clade for profile in samples.values() for clade in profile})
    for column in ABUNDANCES:
        with open(f"{prefix}.{column}.tsv", "w", encoding="utf-8") as sink:
            sink.write("\t".join(["clade_name", *order]) + "\n")
            for clade in clades:
                values = [samples[s].get(clade, {}).get(column) or "0" for s in order]
                sink.write("\t".join([clade, *values]) + "\n")
    empty = [s for s in order if not samples[s]]
    return len(order), len(clades), empty


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--profiles", nargs="+", required=True,
                        help="per-sample sylph profiles (<sample>.sylph.tsv)")
    parser.add_argument("--taxprof", nargs="*", default=[],
                        help="per-sample sylph-tax profiles (<sample>.sylphmpa)")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    rows = stack_genomes(args.profiles, f"{args.prefix}.genomes.tsv")
    print(f"[sylph_merge] {len(args.profiles)} profile(s), {rows} genome call(s)", file=sys.stderr)

    if args.taxprof:
        n, clades, empty = merge_taxprof(args.taxprof, args.prefix)
        print(f"[sylph_merge] {n} taxonomic profile(s), {clades} clade(s)", file=sys.stderr)
        if empty:
            print(f"[sylph_merge] nothing detected in {len(empty)} sample(s): {', '.join(empty)}",
                  file=sys.stderr)


if __name__ == "__main__":
    main()

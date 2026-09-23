#!/usr/bin/env python3
"""
metax_merge.py -- per-sample Metax profiles into cohort tables.

Metax writes <sample>.profile.txt as tab-separated text with NO header row.
The columns, in the order the Metax manual gives them:

  taxon_name, taxid, rank
  reads              read count (scaled when the index was fractional)
  depth              depth of coverage after EM
  abundance          relative abundance, renormalised over the taxa that
                     passed Metax's coverage filters
  breadth            observed breadth of coverage
  expected_breadth   breadth expected under random sampling at that depth
  cov_prob           probability field on the breadth test
  fixed_chunk_breadth, flex_chunk_breadth, expected_flex_chunk_breadth
  chunk_prob         probability field on the chunk-breadth test

and, in pathogen mode (--host), host_names, host_taxids and diseases.
observed/expected breadth is Metax's OEBR, the ratio its filter is built on.

Three tables come out:

  <prefix>.long.tsv       every sample's rows under one named header
  <prefix>.abundance.tsv  taxa by samples, relative abundance
  <prefix>.reads.tsv      taxa by samples, reads

A sample in which nothing passed the filters is a column of zeros in the wide
tables, not a missing column: an empty library must read as empty.
"""
import argparse
import sys
from pathlib import Path

COLUMNS = [
    "taxon_name", "taxid", "rank", "reads", "depth", "abundance",
    "breadth", "expected_breadth", "cov_prob",
    "fixed_chunk_breadth", "flex_chunk_breadth", "expected_flex_chunk_breadth", "chunk_prob",
]
PATHOGEN_COLUMNS = ["host_names", "host_taxids", "diseases"]
RANK_ORDER = ["superkingdom", "domain", "kingdom", "phylum", "class", "order",
              "family", "genus", "species", "strain"]


def sample_id(path):
    name = Path(path).name
    for suffix in (".pathogen.profile.txt", ".profile.txt"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def read_profile(path):
    rows = []
    with open(path, encoding="utf-8") as handle:
        for n, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(COLUMNS):
                sys.exit(f"ERROR: {path}:{n}: {len(fields)} columns, Metax profiles have "
                         f"{len(COLUMNS)} (or {len(COLUMNS) + len(PATHOGEN_COLUMNS)} in pathogen mode).")
            rows.append(fields)
    return rows


def rank_key(rank):
    rank = rank.lower()
    return RANK_ORDER.index(rank) if rank in RANK_ORDER else len(RANK_ORDER)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("--profiles", nargs="+", required=True,
                        help="per-sample Metax profiles (<sample>.profile.txt)")
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()

    samples = {}
    for path in args.profiles:
        sid = sample_id(path)
        if sid in samples:
            sys.exit(f"ERROR: two Metax profiles for sample {sid}.")
        samples[sid] = read_profile(path)
    order = sorted(samples)

    width = max((len(r) for rows in samples.values() for r in rows), default=len(COLUMNS))
    header = COLUMNS + PATHOGEN_COLUMNS[: max(0, width - len(COLUMNS))]
    with open(f"{args.prefix}.long.tsv", "w", encoding="utf-8") as sink:
        sink.write("\t".join(["sample", *header]) + "\n")
        for sid in order:
            for row in samples[sid]:
                sink.write("\t".join([sid, *row, *[""] * (len(header) - len(row))]) + "\n")

    taxa = {}
    values = {"abundance": {}, "reads": {}}
    for sid in order:
        for row in samples[sid]:
            name, taxid, rank = row[0], row[1], row[2]
            taxa.setdefault(taxid, (name, rank))
            values["abundance"][(taxid, sid)] = row[COLUMNS.index("abundance")]
            values["reads"][(taxid, sid)] = row[COLUMNS.index("reads")]
    rows = sorted(taxa, key=lambda t: (rank_key(taxa[t][1]), taxa[t][0], t))
    for kind, table in values.items():
        with open(f"{args.prefix}.{kind}.tsv", "w", encoding="utf-8") as sink:
            sink.write("\t".join(["taxid", "taxon_name", "rank", *order]) + "\n")
            for taxid in rows:
                name, rank = taxa[taxid]
                sink.write("\t".join([taxid, name, rank, *[table.get((taxid, s), "0") for s in order]]) + "\n")

    empty = [s for s in order if not samples[s]]
    print(f"[metax_merge] {len(order)} profile(s), {len(taxa)} taxa", file=sys.stderr)
    if empty:
        print(f"[metax_merge] nothing passed the filters in {len(empty)} sample(s): {', '.join(empty)}",
              file=sys.stderr)


if __name__ == "__main__":
    main()

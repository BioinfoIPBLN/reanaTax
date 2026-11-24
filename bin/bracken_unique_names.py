#!/usr/bin/env python3
"""
bracken_unique_names.py -- one name per taxid across a cohort's Bracken outputs.

Bracken's combine_bracken_outputs.py keys its rows by the name column and stops
("Taxonomy IDs not matching for species ...") when one name belongs to two
taxids. NCBI-derived databases rarely do that; SILVA does at every rank: its
"uncultured" genera are one per family, so a 16S genus-level cohort carrying
two of them never got a combined table.

This reads every Bracken output of the cohort first, finds each name used by
more than one taxid, and writes a copy of every file with those rows renamed
`<name> (taxid <taxid>)`. Every other row is copied unchanged, so a cohort with
no such name gives the same combined table as before. The combined table is
then unique by name as well as by taxid, which the steps keyed on the name
(differential abundance, diversity, decontam, the host-microbe correlation)
rely on.
"""
import argparse
import os
import sys
from collections import defaultdict


def read_rows(path):
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    if not lines:
        sys.exit(f"[bracken_unique_names] '{path}' is empty.")
    header = lines[0].split("\t")
    if header[:2] != ["name", "taxonomy_id"]:
        sys.exit(f"[bracken_unique_names] '{path}' is not a Bracken output: it opens with "
                 f"'{lines[0][:80]}', not 'name<TAB>taxonomy_id'.")
    return lines[0], [line.split("\t") for line in lines[1:] if line.strip()]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", required=True, help="written here under each input's own file name")
    parser.add_argument("files", nargs="+", help="Bracken outputs (<sample>.bracken_<level>.tsv)")
    args = parser.parse_args()

    names = {os.path.basename(path) for path in args.files}
    if len(names) != len(args.files):
        sys.exit("[bracken_unique_names] two inputs share a file name, and the copies would overwrite each other.")

    tables = [(path, *read_rows(path)) for path in args.files]
    taxids = defaultdict(set)
    for _path, _header, rows in tables:
        for row in rows:
            taxids[row[0]].add(row[1])
    shared = {name: ids for name, ids in taxids.items() if len(ids) > 1}

    os.makedirs(args.outdir, exist_ok=True)
    for path, header, rows in tables:
        with open(os.path.join(args.outdir, os.path.basename(path)), "w", encoding="utf-8") as out:
            out.write(header + "\n")
            for row in rows:
                if row[0] in shared:
                    row = [f"{row[0]} (taxid {row[1]})"] + row[1:]
                out.write("\t".join(row) + "\n")

    for name in sorted(shared):
        ids = sorted(shared[name], key=lambda taxid: (len(taxid), taxid))
        print(f"[bracken_unique_names] '{name}' names {len(ids)} taxa ({', '.join(ids)}); "
              f"written as '{name} (taxid <taxid>)'.", file=sys.stderr)


if __name__ == "__main__":
    main()

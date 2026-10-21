#!/usr/bin/env python3
"""
db_lineage.py -- check every genome's taxid against the taxonomy, and write the
lineages sylph-tax needs.

Two jobs, both needing names.dmp and nodes.dmp:

1. A taxid the taxonomy does not know would be built into the database anyway:
   Kraken2 places such sequences nowhere useful and says little about it. So
   every taxid of the genome sheet is looked up here first (through merged.dmp
   when it has been merged into another), and the build stops on any unknown
   one, naming them.

2. sylph reports genomes, not taxa; sylph-tax attaches a lineage to each from a
   metadata file - two columns, genome file name and a GTDB-style string,
   `d__Bacteria;p__...;c__...;o__...;f__...;g__...;s__...` - matched on the file
   name sylph recorded. That file is written here from the NCBI taxonomy, one
   row per genome, so a custom sylph database profiles with taxonomy like a
   prebuilt one.

NCBI ranks map onto the seven prefixes; a rank a lineage lacks is left empty
(`c__`). The domain is NCBI's `domain` (2025 taxonomy) or `superkingdom`
(earlier); viruses, which have neither, get `d__Viruses`.
"""
import argparse
import os
import sys

PREFIXES = [("domain", "d__"), ("phylum", "p__"), ("class", "c__"), ("order", "o__"),
            ("family", "f__"), ("genus", "g__"), ("species", "s__")]
VIRUSES = "10239"


def read_dmp(path, columns):
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = [field.strip() for field in line.rstrip("\t|\n").split("\t|\t")]
            yield [fields[i] for i in columns]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxonomy", required=True, help="directory holding names.dmp and nodes.dmp (and merged.dmp)")
    parser.add_argument("--genomes", required=True, help="TSV: genome file name, taxid")
    parser.add_argument("--output", required=True, help="sylph-tax metadata to write")
    args = parser.parse_args()

    nodes = {taxid: (parent, rank) for taxid, parent, rank in read_dmp(os.path.join(args.taxonomy, "nodes.dmp"), (0, 1, 2))}
    names = {}
    for taxid, name, name_class in read_dmp(os.path.join(args.taxonomy, "names.dmp"), (0, 1, 3)):
        if name_class == "scientific name":
            names[taxid] = name
    merged = {}
    merged_path = os.path.join(args.taxonomy, "merged.dmp")
    if os.path.exists(merged_path):
        merged = {old: new for old, new in read_dmp(merged_path, (0, 1))}

    unknown = []
    rows = []
    with open(args.genomes, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            genome, taxid = line.rstrip("\n").split("\t")[:2]
            taxid = merged.get(taxid, taxid)
            if taxid not in nodes:
                unknown.append(f"{genome} ({taxid})")
                continue
            ranks = {}
            lineage = []
            node = taxid
            while True:
                lineage.append(node)
                parent, rank = nodes[node]
                if rank == "superkingdom":
                    rank = "domain"
                ranks.setdefault(rank, names.get(node, node))
                if parent == node or parent not in nodes:
                    break
                node = parent
            if "domain" not in ranks and VIRUSES in lineage:
                ranks["domain"] = "Viruses"
            rows.append((genome, ";".join(prefix + ranks.get(rank, "") for rank, prefix in PREFIXES)))

    if unknown:
        sys.exit(f"[db_lineage] {len(unknown)} genome(s) have a taxid the taxonomy does not know: "
                 + ", ".join(unknown[:20]) + ". Check the sheet's taxid column against the taxonomy used "
                 "(--db_taxonomy, or the NCBI taxdump downloaded for this run).")

    with open(args.output, "w", encoding="utf-8") as handle:
        for genome, lineage in rows:
            handle.write(f"{genome}\t{lineage}\n")
    print(f"[db_lineage] {len(rows)} genome(s) placed on the taxonomy.", file=sys.stderr)


if __name__ == "__main__":
    main()

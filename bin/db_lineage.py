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

3. With --mpa-output, the same lineages in MetaPhlAn's form, one row per taxid
   (genomes sharing a taxid are one species there): `k__` for the domain,
   names with spaces and punctuation as underscores, `t__TAX<taxid>` as the
   species-level genome bin, and the taxid of each rank beside it. A rank the
   lineage lacks is named after the nearest rank above that it has,
   `g__Thiotrichaceae_unclassified`, with an empty taxid - the convention of
   MetaPhlAn's own database. Viruses are left out: MetaPhlAn profiles them from
   whole viral genomes, not from markers. A eukaryote without a coding-sequence
   file stops the run here, before any gene is called: a prokaryotic gene finder
   would split every intron-containing gene.
"""
import argparse
import os
import re
import sys

MPA_RANKS = [("phylum", "p__"), ("class", "c__"), ("order", "o__"),
             ("family", "f__"), ("genus", "g__"), ("species", "s__")]
PREFIXES = [("domain", "d__"), ("phylum", "p__"), ("class", "c__"), ("order", "o__"),
            ("family", "f__"), ("genus", "g__"), ("species", "s__")]
VIRUSES = "10239"


def read_dmp(path, columns):
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = [field.strip() for field in line.rstrip("\t|\n").split("\t|\t")]
            yield [fields[i] for i in columns]


def mpa_name(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


def mpa_lineage(taxid, lineage_nodes, nodes, names):
    """(lineage, taxid lineage, domain) in MetaPhlAn's form, from root-ward nodes."""
    by_rank = {}
    for node in lineage_nodes:
        rank = nodes[node][1]
        rank = "domain" if rank == "superkingdom" else rank
        by_rank.setdefault(rank, node)
    domain = names.get(by_rank.get("domain", ""), "")
    levels, taxids = [f"k__{mpa_name(domain)}"], [by_rank.get("domain", "")]
    last = mpa_name(domain)
    for rank, prefix in MPA_RANKS:
        node = by_rank.get(rank)
        if node is not None:
            last = mpa_name(names.get(node, node))
            levels.append(prefix + last)
            taxids.append(node)
        else:
            levels.append(f"{prefix}{last}_unclassified")
            taxids.append("")
    levels.append(f"t__TAX{taxid}")
    return "|".join(levels), "|".join(taxids) + "|", domain


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--taxonomy", required=True, help="directory holding names.dmp and nodes.dmp (and merged.dmp)")
    parser.add_argument("--genomes", required=True, help="TSV: genome file name, taxid")
    parser.add_argument("--output", required=True, help="sylph-tax metadata to write")
    parser.add_argument("--mpa-output", default=None,
                        help="also write MetaPhlAn lineages, one row per taxid (needs the third, has-CDS column)")
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
    mpa_rows, no_cds, viruses = {}, [], []
    with open(args.genomes, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            genome, taxid = fields[:2]
            has_cds = len(fields) > 2 and fields[2] == "1"
            original = taxid
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
            if args.mpa_output:
                if VIRUSES in lineage:
                    viruses.append(genome)
                else:
                    mpa = mpa_lineage(original, lineage, nodes, names)
                    if mpa[2] == "Eukaryota" and not has_cds:
                        no_cds.append(genome)
                    mpa_rows[original] = mpa
            rows.append((genome, ";".join(prefix + ranks.get(rank, "") for rank, prefix in PREFIXES)))

    if unknown:
        sys.exit(f"[db_lineage] {len(unknown)} genome(s) have a taxid the taxonomy does not know: "
                 + ", ".join(unknown[:20]) + ". Check the sheet's taxid column against the taxonomy used "
                 "(--db_taxonomy, or the NCBI taxdump downloaded for this run).")

    if no_cds:
        sys.exit(f"[db_lineage] {len(no_cds)} eukaryotic genome(s) have no cds file: " + ", ".join(no_cds[:20])
                 + ". The MetaPhlAn build calls genes with pyrodigal, which is for prokaryotes and splits every "
                 "intron-containing gene; give each eukaryote its coding sequences in the sheet's cds column "
                 "(VEuPathDB, Ensembl and NCBI distribute them), or leave metaphlan out of --build_databases.")

    with open(args.output, "w", encoding="utf-8") as handle:
        for genome, lineage in rows:
            handle.write(f"{genome}\t{lineage}\n")
    if args.mpa_output:
        if viruses:
            print(f"[db_lineage] {len(viruses)} viral genome(s) left out of the MetaPhlAn database: "
                  + ", ".join(viruses[:10]), file=sys.stderr)
        if not mpa_rows:
            sys.exit("[db_lineage] no genome MetaPhlAn can profile from markers: every one is a virus.")
        with open(args.mpa_output, "w", encoding="utf-8") as handle:
            handle.write("taxid\tlineage\ttaxid_lineage\tdomain\n")
            for taxid, (lineage, taxids, domain) in sorted(mpa_rows.items()):
                handle.write(f"{taxid}\t{lineage}\t{taxids}\t{domain}\n")
    print(f"[db_lineage] {len(rows)} genome(s) placed on the taxonomy.", file=sys.stderr)


if __name__ == "__main__":
    main()

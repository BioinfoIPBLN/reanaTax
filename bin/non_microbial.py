#!/usr/bin/env python3
"""
non_microbial.py -- what is not a microbe, listed for removal from the final table.

Host reads that escape host depletion do not all land on the host. Kraken2 puts
them on whatever the database holds that they match: a relative (which
--drop_host_clade covers), or an animal or plant whose assembly carries a
conserved gene or a stretch of contaminating host sequence. With the human or
mouse host and --drop_host_clade order, the final tables still held, per run or
library:

    CLID (human cell lines, 70 runs)         4.5 non-microbial calls, 12.5 microbial
    stroke GSE225948 (mouse single-cell)     590 non-microbial, 65 microbial

mouse, zebrafish and turkey in human cells; human, rabbit, moths, sea urchins and
grasses in mouse brain. --drop_host_clade phylum removes the chordates, but the
invertebrates and plants are outside any rank of the host lineage.

This lists every taxon in the run's own Kraken2 reports that sits under one of
the anchors - by default Metazoa (33208), Viridiplantae (33090), other sequences
(28384: synthetic constructs, vectors) and unclassified sequences (12908:
metagenome placeholders) - so that ABUNDANCE_FILTER removes them with every
other verdict. Bacteria, archaea, viruses, fungi and the protists stay.

WHAT IT COSTS. Helminths are Metazoa: a parasitic worm goes with the moths. So
does a second animal that is the point of the experiment - a blood meal, a
xenograft's mouse cells, a co-culture. That is why it is opt-in.

The lineage comes from the reports - two spaces per level in the name field -
as in host_clade.py, so the anchors are the ones this database uses.
"""

import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from contig_evidence import parse_kreport
from host_clade import read_counts

DEFAULT_ANCHORS = ["33208", "33090", "28384", "12908"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", required=True, help="per-sample Kraken2 reports")
    parser.add_argument("--anchors", nargs="+", default=DEFAULT_ANCHORS,
                        help="drop everything under these taxids (default: Metazoa, Viridiplantae, "
                             "other sequences, unclassified sequences)")
    parser.add_argument("--drop-list", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--mqc")
    args = parser.parse_args()

    nodes = {}
    for path in args.reports:
        nodes.update(parse_kreport(path))
    counts = read_counts(args.reports)

    anchors = [a for a in args.anchors if a in nodes]
    members = {}
    for taxid, (_rank, _name, ancestors) in nodes.items():
        hit = next((a for a in anchors if taxid == a or a in ancestors), None)
        if hit is not None:
            members[taxid] = hit
    ordered = sorted(members, key=lambda t: -counts.get(t, 0.0))

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in ordered:
            handle.write(f"{taxid}\n")

    with open(args.evidence, "w", encoding="utf-8") as handle:
        handle.write("taxid\tname\trank\treads\tanchor\n")
        for taxid in ordered:
            rank, name, _a = nodes[taxid]
            handle.write(f"{taxid}\t{name}\t{rank}\t{counts.get(taxid, 0.0):.0f}\t{nodes[members[taxid]][1]}\n")

    per_anchor = collections.Counter()
    for taxid, anchor in members.items():
        per_anchor[nodes[anchor][1]] += counts.get(taxid, 0.0)
    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_non_microbial'\n"
                "# section_name: 'Non-microbial taxa'\n"
                "# description: 'Reads on animals, plants and synthetic or unclassified\n"
                "#     sequences, removed from the final tables by --drop_non_microbial.\n"
                "#     Mostly host reads that escaped depletion and matched another\n"
                "#     eukaryote; helminths are animals too and go with them.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_non_microbial_plot'\n"
                "#     title: 'reanaTax: non-microbial taxa removed'\n"
                "#     ylab: 'Reads'\n"
            )
            handle.write("Sample\tdropped_reads\n")
            for name, total in sorted(per_anchor.items()):
                handle.write(f"{name}\t{total:.0f}\n")

    print(
        f"[non_microbial] {len(members)} taxid(s) under "
        + ", ".join(f"{nodes[a][1]} ({a})" for a in anchors)
        + f" listed for removal, {sum(per_anchor.values()):.0f} read(s).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

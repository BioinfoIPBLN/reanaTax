#!/usr/bin/env python3
"""
host_clade.py -- the host's relatives are host reads wearing another name.

--drop_host_taxon removes one taxid. That is not where host leakage ends up: a
human read that misses GRCh38 does not vanish, it gets assigned to the nearest
thing in the database that it does match, and the nearest thing is a relative.
Measured on the CSI-Microbes 10x cohort, with Homo sapiens already dropped:

    clade                 uninfected   heat-killed   infected   infected(5')
    genus  Homo                 0.1%          0.1%       0.0%           0.0%
    family Hominidae           65.7%         60.8%      58.1%          52.4%
    order  Primates            78.8%         74.3%      72.1%          67.7%
    class  Mammalia            88.3%         88.5%      84.4%          71.0%
    phylum Chordata            93.4%         92.8%      88.6%          71.2%

...as a share of everything the pipeline was still calling microbial. Pan,
Pongo, Gorilla, Macaca, Tupaia. The genus rank catches almost nothing because
Homo has no other species in the database; the family rank is where it starts.

WHAT THIS IS NOT. It is a prediction from the taxonomy, not a measurement of the
reads, so it cannot see leakage that lands outside the clade - the same cohort
carries turbot, grouper, spruce and Naegleria, which no rank of the host lineage
contains. Those are conserved or low-complexity matches, and --host_kmer_filter
is the filter that judges them, by asking whether a taxon's reads are mostly
host k-mers. The two are complements: this one is free and deterministic, that
one is evidence-based and costs a pass over the read-level output.

THE RISK IS REAL AT HIGH RANKS. Dropping class Mammalia is harmless when the
sample is human tissue and the question is bacterial. It is destructive on a
xenograft, where mouse reads are biology rather than noise, and on any cohort
where a second vertebrate is the point. The rank is therefore a choice the user
makes, never a default.

The lineage comes from the run's own Kraken2 reports - two spaces per level in
the name field - so the ranks are the ones the database actually uses and no
external taxonomy dump has to be kept in step with it.
"""

import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from contig_evidence import open_maybe_gzip, parse_kreport

# Kraken2's major rank codes. The sub-ranks it interleaves (P1, C3, O4 ...) are
# deliberately not selectable: they are database bookkeeping, not ranks a user
# can reason about, and "order" should mean Primates rather than Simiiformes.
RANK_CODES = {"genus": "G", "family": "F", "order": "O", "class": "C", "phylum": "P"}


def read_counts(paths):
    """{taxid: direct reads summed over every report}."""
    counts = collections.defaultdict(float)
    for path in paths:
        with open_maybe_gzip(path) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                taxid = fields[-2].strip()
                if not taxid.isdigit():
                    continue
                try:
                    counts[taxid] += float(fields[2])
                except ValueError:
                    continue
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", required=True, help="per-sample Kraken2 reports")
    parser.add_argument("--taxid", required=True, help="the host taxon, e.g. 9606")
    parser.add_argument("--rank", required=True, choices=sorted(RANK_CODES),
                        help="drop everything under the host's ancestor at this rank")
    parser.add_argument("--drop-list", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--mqc")
    args = parser.parse_args()

    nodes = {}
    for path in args.reports:
        nodes.update(parse_kreport(path))

    host = str(args.taxid).strip()
    if host not in nodes:
        raise SystemExit(
            f"host_clade: taxid {host} appears in none of the {len(args.reports)} Kraken2 "
            "report(s), so its lineage cannot be read from them. Either the host taxid is "
            "wrong, or host depletion was clean enough that nothing was assigned to it - in "
            "which case there is no leakage to expand and --drop_host_clade should be unset."
        )

    wanted = RANK_CODES[args.rank]
    _rank, host_name, lineage = nodes[host]
    anchor = next((t for t in reversed(lineage) if nodes[t][0] == wanted), None)
    if anchor is None and nodes[host][0] == wanted:
        anchor = host
    if anchor is None:
        raise SystemExit(
            f"host_clade: {host_name} ({host}) has no ancestor at rank '{args.rank}' in this "
            "database's taxonomy. Its lineage is: "
            + " > ".join(f"{nodes[t][1]} [{nodes[t][0]}]" for t in lineage)
        )

    counts = read_counts(args.reports)
    members = [
        taxid
        for taxid, (_r, _n, ancestors) in nodes.items()
        if taxid == anchor or anchor in ancestors
    ]
    members.sort(key=lambda t: -counts.get(t, 0.0))

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in members:
            handle.write(f"{taxid}\n")

    anchor_name = nodes[anchor][1]
    with open(args.evidence, "w", encoding="utf-8") as handle:
        handle.write("taxid\tname\trank\treads\tclade_anchor\tanchor_rank\n")
        for taxid in members:
            rank, name, _a = nodes[taxid]
            handle.write(
                f"{taxid}\t{name}\t{rank}\t{counts.get(taxid, 0.0):.0f}\t{anchor_name}\t{args.rank}\n"
            )

    total = sum(counts.get(t, 0.0) for t in members)
    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_host_clade'\n"
                "# section_name: 'Host clade leakage'\n"
                "# description: 'Taxa removed as the host wearing another name. A host read\n"
                "#     that misses the host genome is assigned to the nearest relative in the\n"
                "#     database, so everything under the host lineage at the chosen rank is\n"
                "#     dropped. This is a prediction from the taxonomy, not a measurement of\n"
                "#     the reads: leakage landing outside the clade is what --host_kmer_filter\n"
                "#     is for.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_host_clade_plot'\n"
                "#     title: 'reanaTax: host clade leakage'\n"
                "#     ylab: 'Reads'\n"
            )
            handle.write("Sample\tdropped_reads\n")
            handle.write(f"{anchor_name} ({args.rank})\t{total:.0f}\n")

    print(
        f"[host_clade] {host_name} ({host}) -> {args.rank} {anchor_name} ({anchor}): "
        f"{len(members)} taxid(s) listed for removal, {total:.0f} read(s).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

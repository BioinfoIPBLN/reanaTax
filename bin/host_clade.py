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

SEVERAL HOSTS. --taxid takes a comma-separated list, as --host_kmer_taxid does:
a species complex (An. gambiae with An. coluzzii), or two hosts. Each listed
taxid contributes its clade at the rank, and the drop list is their union.

A HOST WITH NO READS OF ITS OWN. A host taxid can be missing from every report,
and that is not a reason to stop. On a host-free mock there is nothing to drop.
But a host whose genome is thin in the database - a tick, a sand fly - has its
reads assigned to relatives instead, which is the very leakage this removes. So
the host's lineage is then read from the database's nodes.dmp (--nodes), and its
clade is dropped wherever it does appear in the reports. Only a taxid neither
source knows is skipped, with a note; if every one is, the lists are empty.
"""

import argparse
import collections
import functools
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


@functools.lru_cache(maxsize=None)
def read_nodes(path):
    """{taxid: parent}, {taxid: rank} from an NCBI nodes.dmp."""
    parent, ranks = {}, {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            fields = [field.strip() for field in line.split("|", 3)]
            if len(fields) > 2:
                parent[fields[0]], ranks[fields[0]] = fields[1], fields[2]
    return parent, ranks


def nodes_ancestor(path, taxid, rank):
    """(known, ancestor): whether nodes.dmp knows taxid, and its ancestor (or itself) at rank."""
    parent, ranks = read_nodes(path)
    if taxid not in parent:
        return False, None
    node, seen = taxid, set()
    while node in parent and node not in seen:
        if ranks[node] == rank:
            return True, node
        seen.add(node)
        node = parent[node]
    return True, None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", nargs="+", required=True, help="per-sample Kraken2 reports")
    parser.add_argument("--taxid", required=True,
                        help="the host taxon, or several comma-separated, e.g. 9606 or 7165,44542")
    parser.add_argument("--rank", required=True, choices=sorted(RANK_CODES),
                        help="drop everything under the host's ancestor at this rank")
    parser.add_argument("--drop-list", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--mqc")
    parser.add_argument("--nodes", help="the database's nodes.dmp, for a host taxid no report contains")
    args = parser.parse_args()

    nodes = {}
    for path in args.reports:
        nodes.update(parse_kreport(path))

    wanted = RANK_CODES[args.rank]
    anchors, notes = [], []
    for host in [part.strip() for part in str(args.taxid).split(",") if part.strip()]:
        if host in nodes:
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
        elif args.nodes:
            known, anchor = nodes_ancestor(args.nodes, host, args.rank)
            if not known:
                notes.append(f"taxid {host} is in no report and not in {args.nodes}; skipped")
                continue
            if anchor is None:
                raise SystemExit(f"host_clade: taxid {host} has no ancestor at rank '{args.rank}' in {args.nodes}.")
            notes.append(f"taxid {host} is in no report; its {args.rank} ({anchor}) was read from {args.nodes}")
        else:
            notes.append(f"taxid {host} is in no report and no --nodes was given; skipped")
            continue
        if anchor not in anchors:
            anchors.append(anchor)

    for note in notes:
        print(f"[host_clade] {note}", file=sys.stderr)

    counts = read_counts(args.reports)
    clade_of = {}
    for taxid, (_r, _n, ancestors) in nodes.items():
        anchor = next((a for a in anchors if taxid == a or a in ancestors), None)
        if anchor is not None:
            clade_of[taxid] = anchor
    members = sorted(clade_of, key=lambda t: -counts.get(t, 0.0))
    label = {a: nodes[a][1] if a in nodes else f"taxid {a}" for a in anchors}

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in members:
            handle.write(f"{taxid}\n")

    with open(args.evidence, "w", encoding="utf-8") as handle:
        handle.write("taxid\tname\trank\treads\tclade_anchor\tanchor_rank\n")
        for taxid in members:
            rank, name, _a = nodes[taxid]
            handle.write(
                f"{taxid}\t{name}\t{rank}\t{counts.get(taxid, 0.0):.0f}\t{label[clade_of[taxid]]}\t{args.rank}\n"
            )

    total = sum(counts.get(t, 0.0) for t in members)
    if not members:
        print(f"[host_clade] nothing to drop: no report holds a taxon under the host's {args.rank}.", file=sys.stderr)
        return
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
            for anchor in anchors:
                if anchor in clade_of.values():
                    reads = sum(counts.get(t, 0.0) for t in members if clade_of[t] == anchor)
                    handle.write(f"{label[anchor]} ({args.rank})\t{reads:.0f}\n")

    print(
        f"[host_clade] {args.rank} {', '.join(f'{label[a]} ({a})' for a in anchors if a in clade_of.values())}: "
        f"{len(members)} taxid(s) listed for removal, {total:.0f} read(s).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

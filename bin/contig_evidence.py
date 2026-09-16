#!/usr/bin/env python3
"""
contig_evidence.py -- what the assembly says about each taxon the reads claimed.

Every other filter in this pipeline scores a taxon on read counts and k-mer
statistics. None of them ever produces a longer sequence, and that is the gap
this fills: a 1 kb contig that classifies to a taxon is categorically stronger
evidence than fifty 100 bp reads that do, because a false positive built from
index hopping or a conserved locus has nothing to assemble.

Three things make the accounting less obvious than it looks.

CONTIGS MUST ROLL UP. A contig assigned to Fusobacterium nucleatum is support
for genus Fusobacterium, and on the CSI-Microbes plate the genus row is where
the false positives lived. The ancestry comes from the contig kreport's own
indentation - two spaces per level - so no taxonomy dump is needed and the
lineage is exactly the one Kraken2 used.

REPORTS COME IN TWO SHAPES. With --kraken2_report_minimizer_data a kreport has
eight columns rather than six, and reading rank, taxid and name by fixed index
would put the minimizer counts where the rank belongs. Both are accepted.

ABSENCE IS WEAK EVIDENCE. A genuinely rare organism will not assemble either,
which is why nothing is dropped unless it had enough reads to have assembled
something: --min-reads is a floor below which a taxon is recorded and not
judged, the same contract bin/minimizer_filter.py uses.
"""

import argparse
import gzip
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from filter_abundance import read_table, taxid_column


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def parse_kreport(path):
    """{taxid: (rank, name, [ancestor taxids, outermost first])} from one kreport.

    Kraken2 writes six columns, or eight with --report-minimizer-data; the rank,
    taxid and name are the last three either way. Depth is two spaces per level
    in the name field, which is what makes the hierarchy recoverable.
    """
    nodes = {}
    stack = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            rank, taxid, name = fields[-3], fields[-2].strip(), fields[-1]
            if not taxid.isdigit():
                continue
            depth = (len(name) - len(name.lstrip(" "))) // 2
            del stack[depth:]
            nodes[taxid] = (rank.strip(), name.strip(), list(stack))
            stack.append(taxid)
    return nodes


def parse_assignments(path):
    """[(taxid, length)] for every classified contig in one Kraken2 output."""
    contigs = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] != "C":
                continue
            taxid = fields[2].strip()
            # Kraken2 writes `name (taxid 1234)` when --use-names is on.
            if not taxid.isdigit() and "taxid" in taxid:
                taxid = taxid.rsplit("taxid", 1)[1].strip(" )")
            if not taxid.isdigit():
                continue
            try:
                length = int(fields[3].split("|")[0])
            except ValueError:
                continue
            contigs.append((taxid, length))
    return contigs


def n50(lengths):
    if not lengths:
        return 0
    ordered = sorted(lengths, reverse=True)
    half = sum(ordered) / 2.0
    running = 0
    for length in ordered:
        running += length
        if running >= half:
            return length
    return ordered[-1]


def read_totals(path):
    """{taxid: clade reads} from a combined kreport, for the read side of the check."""
    _comments, _header_line, header, rows = read_table(path)
    column = taxid_column(header)
    if column is None:
        raise SystemExit(f"contig_evidence: {path} has no taxid column.")
    totals_column = None
    for candidate in ("tot_all", "tot_frac", "new_est_reads", "total_reads"):
        for i, name in enumerate(header):
            if name.strip().lower() == candidate:
                totals_column = i
                break
        if totals_column is not None:
            break
    if totals_column is None:
        raise SystemExit(
            f"contig_evidence: {path} has no total-reads column "
            "(looked for tot_all, new_est_reads, total_reads)."
        )
    totals = {}
    for row in rows:
        if column >= len(row) or totals_column >= len(row):
            continue
        taxid = row[column].strip()
        try:
            totals[taxid] = totals.get(taxid, 0) + int(float(row[totals_column]))
        except ValueError:
            continue
    return totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", nargs="+", required=True, help="Kraken2 per-contig output, one per pool")
    parser.add_argument("--reports", nargs="+", required=True, help="Kraken2 report of the contigs, one per pool")
    parser.add_argument("--table", required=True, help="filtered combined kreport, for the read counts")
    parser.add_argument("--output", required=True, help="the evidence table")
    parser.add_argument("--drop-list", help="taxa with reads but no contig support")
    parser.add_argument("--mqc", help="a MultiQC summary table")
    parser.add_argument("--min-reads", type=int, default=50,
                        help="below this a taxon is recorded and not judged (default: 50)")
    parser.add_argument("--min-length", type=int, default=500,
                        help="a contig must reach this to count as support (default: 500)")
    parser.add_argument("--min-contigs", type=int, default=1,
                        help="this many qualifying contigs are needed (default: 1)")
    args = parser.parse_args()

    nodes = {}
    for path in args.reports:
        nodes.update(parse_kreport(path))

    # Credit every contig to its own taxon and to each of its ancestors, so a
    # species-level contig counts as support for the genus row above it.
    lengths = {}
    pools = {}
    for path in args.assignments:
        pool = os.path.basename(path)
        for taxid, length in parse_assignments(path):
            for owner in list(nodes.get(taxid, ("", "", []))[2]) + [taxid]:
                lengths.setdefault(owner, []).append(length)
                pools.setdefault(owner, set()).add(pool)

    totals = read_totals(args.table)

    judged = supported = unsupported = 0
    dropped = []
    rows = []
    for taxid, reads in sorted(totals.items(), key=lambda item: -item[1]):
        rank, name, _ancestors = nodes.get(taxid, ("", "", []))
        own = lengths.get(taxid, [])
        qualifying = [length for length in own if length >= args.min_length]
        if reads < args.min_reads:
            verdict = "not_judged"
        elif len(qualifying) >= args.min_contigs:
            verdict = "supported"
            supported += 1
            judged += 1
        else:
            verdict = "unsupported"
            unsupported += 1
            judged += 1
            dropped.append(taxid)
        rows.append(
            [
                taxid,
                name,
                rank,
                str(reads),
                str(len(pools.get(taxid, ()))),
                str(len(own)),
                str(len(qualifying)),
                str(sum(own)),
                str(max(own) if own else 0),
                str(n50(own)),
                f"{statistics.fmean(own):.1f}" if own else "0.0",
                verdict,
            ]
        )

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(
            "taxid\tname\trank\treads\tpools\tcontigs\tcontigs_over_min\t"
            "contig_bp\tlongest_bp\tn50\tmean_bp\tverdict\n"
        )
        for row in rows:
            handle.write("\t".join(row) + "\n")

    if args.drop_list:
        with open(args.drop_list, "w", encoding="utf-8") as handle:
            for taxid in dropped:
                handle.write(f"{taxid}\n")

    # Every classified contig rolls up to root, so root's tally is the total.
    classified = len(lengths.get("1", []))

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_contig_evidence'\n"
                "# section_name: 'Contig evidence'\n"
                "# description: 'Taxa checked against a de novo assembly of the non-host\n"
                "#     fraction. A taxon is judged once it has enough reads to have assembled\n"
                "#     something, and is supported when a long enough contig classifies to it or\n"
                "#     to anything below it. Absence of a contig is weak evidence of absence - a\n"
                "#     genuinely rare organism does not assemble either - which is what the read\n"
                "#     gate is for.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_contig_evidence_plot'\n"
                "#     title: 'reanaTax: contig evidence'\n"
                "#     ylab: 'Taxa'\n"
            )
            handle.write("Sample\tsupported\tunsupported\tnot_judged\n")
            handle.write(f"all taxa\t{supported}\t{unsupported}\t{len(rows) - judged}\n")

    print(
        f"[contig_evidence] {len(rows)} taxa: {supported} supported, {unsupported} unsupported, "
        f"{len(rows) - judged} not judged (gate {args.min_reads} reads, "
        f"{args.min_contigs} contig(s) >= {args.min_length} bp); "
        f"{classified} classified contig(s) across {len(args.assignments)} pool(s); "
        f"{len(dropped)} taxid(s) listed for removal.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

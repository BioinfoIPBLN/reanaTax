#!/usr/bin/env python3
"""
rrna_share.py -- what share of each taxon's reads is ribosomal RNA.

The minimizer filter judges a taxon on the breadth of its evidence: how many of
the genome's distinct minimizers its reads reach. Reads from rRNA pile onto an
operon of about 5 kb, so a taxon seen mostly through its rRNA looks narrow -
the same shape as a contaminated reference contig - and can be dropped as low
coverage. On GSE152924 (mosquito RNA-seq, mostly mRNA) two low-abundance
Firmicutes were dropped with 46% and 62% of their reads rRNA, while the kept
bacteria were at 0-3%; on total-RNA libraries rRNA can be most of what a
bacterium contributes. Nothing here changes a verdict: the table puts each
taxon's rRNA share beside the minimizer verdict, so a drop can be read for what
it is.

  sample   one library: the rRNA read IDs SortMeRNA found among the reads
           Kraken2 classified, joined to Kraken2's per-read calls and rolled up
           the library's own report hierarchy, so a genus counts the reads of
           its species exactly as the clade counts the filter judges do.
  merge    the cohort: per-sample tables summed, with the minimizer filter's
           evidence joined on taxid when there is one.
"""
import argparse
import csv
import gzip
import re
import sys

NAMED = re.compile(r"\(taxid (\d+)\)\s*$")


def opener(path):
    with open(path, "rb") as probe:
        magic = probe.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path)


def read_id(token):
    token = token.split()[0]
    return token[:-2] if token.endswith(("/1", "/2")) else token


def report_tree(path):
    """Parent, rank and name of every taxon in a Kraken2 report, from its indentation."""
    parent, rank, name = {}, {}, {}
    stack = []
    with opener(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 8:
                taxid, taxon_rank, label = fields[6], fields[5], fields[7]
            elif len(fields) == 6:
                taxid, taxon_rank, label = fields[4], fields[3], fields[5]
            else:
                continue
            depth = (len(label) - len(label.lstrip(" "))) // 2
            taxid = taxid.strip()
            while len(stack) > depth:
                stack.pop()
            parent[taxid] = stack[-1] if stack else None
            rank[taxid] = taxon_rank.strip()
            name[taxid] = label.strip()
            stack.append(taxid)
    return parent, rank, name


def sample(args):
    rrna = set()
    with opener(args.ids) as handle:
        for line in handle:
            if line.strip():
                rrna.add(read_id(line))

    direct_total, direct_rrna = {}, {}
    with opener(args.assignments) as handle:
        for line in handle:
            fields = line.split("\t", 3)
            if len(fields) < 3 or fields[0] != "C":
                continue
            match = NAMED.search(fields[2])
            taxid = match.group(1) if match else fields[2].strip()
            direct_total[taxid] = direct_total.get(taxid, 0) + 1
            if read_id(fields[1]) in rrna:
                direct_rrna[taxid] = direct_rrna.get(taxid, 0) + 1

    parent, rank, name = report_tree(args.report)
    clade_total, clade_rrna = {}, {}
    for taxid, count in direct_total.items():
        hits = direct_rrna.get(taxid, 0)
        node, seen = taxid, set()
        while node is not None and node not in seen:
            seen.add(node)
            clade_total[node] = clade_total.get(node, 0) + count
            clade_rrna[node] = clade_rrna.get(node, 0) + hits
            node = parent.get(node)

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "taxid", "rank", "name", "reads", "rrna_reads"])
        for taxid in sorted(clade_total, key=lambda t: -clade_total[t]):
            writer.writerow([args.sample, taxid, rank.get(taxid, ""), name.get(taxid, ""),
                             clade_total[taxid], clade_rrna[taxid]])
    print(f"[rrna_share] {args.sample}: {len(rrna)} rRNA read(s) among "
          f"{sum(direct_total.values())} classified; {len(clade_total)} taxa.", file=sys.stderr)


def merge(args):
    totals, labels, per_sample = {}, {}, {}
    for path in args.tables:
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                taxid = row["taxid"]
                reads, hits = int(row["reads"]), int(row["rrna_reads"])
                cohort = totals.setdefault(taxid, [0, 0, 0])
                cohort[0] += reads
                cohort[1] += hits
                cohort[2] += 1
                labels.setdefault(taxid, (row["rank"], row["name"]))
                per_sample.setdefault(taxid, []).append(hits / reads if reads else 0.0)

    evidence = {}
    if args.evidence:
        with open(args.evidence) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                evidence[row["taxid"]] = row

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "rank", "name", "samples", "reads", "rrna_reads", "rrna_fraction",
                         "max_sample_rrna_fraction", "minimizer_verdict", "minimizer_coverage",
                         "minimizer_duplication"])
        for taxid in sorted(totals, key=lambda t: -totals[t][0]):
            reads, hits, samples = totals[taxid]
            if reads < args.min_reads:
                continue
            row = evidence.get(taxid, {})
            writer.writerow([
                taxid, labels[taxid][0], labels[taxid][1], samples, reads, hits,
                f"{hits / reads:.4f}" if reads else "NA",
                f"{max(per_sample[taxid]):.4f}",
                row.get("verdict", "NA"), row.get("coverage", "NA"), row.get("duplication", "NA"),
            ])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    one = commands.add_parser("sample")
    one.add_argument("--sample", required=True)
    one.add_argument("--ids", required=True, help="rRNA read IDs, one per line, plain or gzipped")
    one.add_argument("--assignments", required=True, help="Kraken2 per-read output")
    one.add_argument("--report", required=True, help="the same library's Kraken2 report")
    one.add_argument("--output", required=True)
    cohort = commands.add_parser("merge")
    cohort.add_argument("--tables", nargs="+", required=True)
    cohort.add_argument("--evidence", help="the minimizer filter's evidence table")
    cohort.add_argument("--min-reads", type=int, default=10,
                        help="leave out taxa with fewer reads than this across the cohort")
    cohort.add_argument("--output", required=True)
    args = parser.parse_args()
    sample(args) if args.command == "sample" else merge(args)


if __name__ == "__main__":
    main()

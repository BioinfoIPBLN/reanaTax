#!/usr/bin/env python3
"""
sparcc_network.py -- FastSpar's two square matrices as one edge table.

FastSpar writes a correlation matrix and a p-value matrix, both taxa by taxa and
both symmetric. That is the wrong shape for reading, for filtering and for
anything downstream, so this flattens the upper triangle into one row per pair,
attaches the taxon names, and adjusts the p-values across the pairs actually
tested.

The adjustment is the point of doing it here rather than in a spreadsheet.
FastSpar's p-values are per pair and uncorrected; a 200-taxon network is 19,900
pairs, so at a nominal 0.05 roughly a thousand edges are expected by chance
alone. An unadjusted SparCC network is mostly noise arranged in a suggestive
shape, which is exactly how co-occurrence networks acquire their reputation.

Benjamini-Hochberg across every pair, once. Not per taxon, and not per sign: an
edge is one test, and splitting the family by sign would let the positive edges
be judged against a smaller family than the negative ones for no reason but
their sign.

The p-value floor is reported alongside the result, because a bootstrap p-value
cannot go below 1/replicates and with enough pairs that floor can sit above what
BH requires. bin/sparcc_prepare.py refuses to start such a run; this repeats the
number so it is visible in the output as well as in the log.
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import benjamini_hochberg


def read_matrix(path):
    """(ids, {(i, j): value}) from a FastSpar square matrix."""
    with open(path, encoding="utf-8") as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not rows:
        sys.exit(f"[sparcc_network] '{path}' is empty.")
    ids = [field.strip() for field in rows[0][1:]]
    values = {}
    for row in rows[1:]:
        name = row[0].strip()
        if name not in ids:
            continue
        i = ids.index(name)
        for j, field in enumerate(row[1:]):
            try:
                values[(i, j)] = float(field)
            except ValueError:
                values[(i, j)] = float("nan")
    return ids, values


def read_taxa(path):
    mapping = {}
    if not path or not os.path.exists(path):
        return mapping
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["id"]] = (row.get("taxid", ""), row.get("name", ""))
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--correlation", required=True)
    parser.add_argument("--pvalues", required=True)
    parser.add_argument("--taxa", default=None, help="the id/taxid/name map from sparcc_prepare.py")
    parser.add_argument("--permutations", type=int, default=1000)
    parser.add_argument("--min-correlation", type=float, default=0.3,
                        help="|rho| an edge needs before it is reported at all")
    parser.add_argument("--p-threshold", type=float, default=0.05, help="on the adjusted p")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    ids, correlation = read_matrix(args.correlation)
    p_ids, pvalues = read_matrix(args.pvalues)
    if p_ids != ids:
        sys.exit(
            "[sparcc_network] the correlation and p-value matrices are over different taxa, so "
            "their cells do not correspond. They must come from the same FastSpar run."
        )
    taxa = read_taxa(args.taxa)

    edges = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            rho = correlation.get((i, j))
            p = pvalues.get((i, j))
            if rho is None or rho != rho:
                continue
            edges.append({
                "taxon_a": taxa.get(ids[i], ("", ids[i]))[1] or ids[i],
                "taxon_b": taxa.get(ids[j], ("", ids[j]))[1] or ids[j],
                "taxid_a": taxa.get(ids[i], ("", ""))[0],
                "taxid_b": taxa.get(ids[j], ("", ""))[0],
                "rho": round(rho, 4),
                "p": p if p is not None and p == p else None,
            })

    for edge, q in zip(edges, benjamini_hochberg([edge["p"] for edge in edges])):
        edge["q"] = q

    floor = 1.0 / args.permutations if args.permutations else 0.0
    significant = [
        edge for edge in edges
        if edge["q"] is not None and edge["q"] < args.p_threshold
        and abs(edge["rho"]) >= args.min_correlation
    ]

    columns = ["taxon_a", "taxon_b", "taxid_a", "taxid_b", "rho", "p", "q"]
    edges.sort(key=lambda edge: (edge["q"] if edge["q"] is not None else 1.0, -abs(edge["rho"])))
    for name, subset in (("edges_all", edges), ("edges_significant", significant)):
        with open(f"{args.prefix}.sparcc_{name}.tsv", "w", encoding="utf-8") as handle:
            handle.write("\t".join(columns) + "\n")
            for edge in subset:
                handle.write("\t".join(
                    "NA" if edge[column] is None else
                    (f"{edge[column]:.6g}" if column in ("p", "q") else str(edge[column]))
                    for column in columns
                ) + "\n")

    positive = sum(1 for edge in significant if edge["rho"] > 0)
    with open(f"{args.prefix}_sparcc_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_sparcc'",
            "# section_name: 'Taxon co-occurrence (SparCC)'",
            "# description: 'Correlations between taxa inferred by FastSpar. Taxonomic profiles are",
            "#     compositional - the counts sum to a library size that has nothing to do with the",
            "#     biology - so an ordinary correlation matrix manufactures negative associations",
            "#     everywhere; SparCC works on log-ratio variances, which are invariant to the",
            "#     total. p-values come from bootstrap resampling and are Benjamini-Hochberg",
            f"#     adjusted across all {len(edges)} pairs at once, because at a nominal 0.05 a",
            "#     network this size would otherwise carry hundreds of edges by chance alone.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_sparcc_plot'",
            "#     title: 'reanaTax: SparCC edges'",
            "#     ylab: 'Edges'",
            "Sample\tPositive\tNegative\tNot significant",
            f"all pairs\t{positive}\t{len(significant) - positive}\t{len(edges) - len(significant)}",
            "",
        ]))

    print(
        f"[sparcc_network] {len(ids)} taxa, {len(edges)} pairs; {len(significant)} edge(s) at "
        f"q < {args.p_threshold:g} and |rho| >= {args.min_correlation:g} "
        f"({positive} positive). Smallest attainable p from {args.permutations} replicates: "
        f"{floor:.3g}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

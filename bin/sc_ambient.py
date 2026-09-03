#!/usr/bin/env python3
"""
sc_ambient.py -- is this taxon associated with cells, or is it in the medium?

THIS DOES NOT DECIDE WHAT IS CONTAMINATION, and the distinction is the whole
reason it exists as an annotation rather than a filter.

An empty droplet contains at least four things a droplet cannot tell apart:
reagent and kit contaminants; genuinely extracellular microbes that were in the
tissue suspension; microbes released from cells lysed during dissociation; and
ambient nucleic-acid soup. A luminal, mucosal or biofilm organism - which for a
gut, skin, oral or vulvar sample is very often the organism the study is about -
lands in that pool by construction, not by artefact. Calling everything in the
empty droplets a contaminant would delete it.

What the comparison CAN answer is narrower and still useful: is this taxon
spatially associated with cells? Three outcomes, and only one of them is a
statement about contamination:

    ratio >> 1   cell-associated - intracellular, or tightly adherent. The
                 strongest evidence droplet data can offer.
    ratio ~= 1   ambient. Reagent contaminant OR genuine extracellular
                 organism. UNDECIDABLE from this data. Not a verdict.
    ratio << 1   depleted in cells - reagent, index hopping, or something that
                 does not survive the cells it came with.

The middle row does not resolve here and cannot be made to. The only thing that
separates a kit contaminant from a real extracellular organism is an external
measurement of the kit, which is what --decontam's blanks provide and what
nothing else in this pipeline does. Read the two together.

Two ratios are computed, because they can disagree and the disagreement is
informative:

  prevalence_ratio  share of cells carrying the taxon, over the share of empty
                    droplets carrying it. Robust, ignores how much is there.
  rate_ratio        the taxon's UMIs per host UMI in cells, over the same in
                    empty droplets. Sensitive to abundance, and the one that
                    moves when a taxon is present everywhere but concentrated
                    in cells.

The enrichment test is Fisher's exact on the 2x2 presence table, one-sided in
whichever direction the ratio points, adjusted across taxa with Benjamini-
Hochberg. Fisher rather than chi-square for the reason CSI-Microbes gives: the
expected counts are small, and the chi-square approximation needs them above 5.

`--drop-ambient` writes a drop list. It is off by default and it is only
correct when the question is about INTRACELLULAR microbes. If the study is
about the community in the tissue, ambient taxa are part of the answer.
"""
import argparse
import csv
import gzip
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import benjamini_hochberg, fisher_exact_greater


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(
        path, encoding="utf-8", errors="replace"
    )


def find_file(directory, *stems):
    for stem in stems:
        for candidate in (stem, stem + ".gz"):
            path = os.path.join(directory, candidate)
            if os.path.exists(path):
                return path
    return None


def locate(root, feature_type, filtering):
    """The MatrixMarket directory for one filtering level of a Solo.out tree."""
    preferred = os.path.join(root, feature_type, filtering)
    if find_file(preferred, "barcodes.tsv"):
        return preferred
    for base, _dirs, files in os.walk(root):
        if os.path.basename(base) == filtering and any(f.startswith("barcodes.tsv") for f in files):
            return base
    sys.exit(
        f"[sc_ambient] no {filtering}/ barcodes under '{root}'. This needs BOTH halves of a "
        f"STARsolo Solo.out tree: {feature_type}/filtered/ names the cells, {feature_type}/raw/ "
        "names every barcode including the empty droplets. Without raw/ there is nothing to "
        "compare cells against."
    )


def read_barcodes(directory):
    path = find_file(directory, "barcodes.tsv")
    if path is None:
        sys.exit(f"[sc_ambient] '{directory}' has no barcodes.tsv.")
    out = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            value = line.rstrip("\n").split("\t")[0]
            if value:
                out.append(value)
    return out


def host_umis_per_barcode(directory, barcodes):
    """Column sums of the matrix, streamed.

    The raw matrix of a 10x run is hundreds of millions of nonzeros and tens of
    gigabytes unpacked. Nothing here needs the matrix itself - only the total
    per barcode - so it is summed in one pass and never held.
    """
    path = find_file(directory, "matrix.mtx")
    if path is None:
        sys.exit(f"[sc_ambient] '{directory}' has no matrix.mtx.")
    totals = [0.0] * len(barcodes)
    with open_maybe_gzip(path) as handle:
        header = None
        for line in handle:
            if line.startswith("%"):
                continue
            if header is None:
                header = line.split()
                if int(header[1]) != len(barcodes):
                    sys.exit(
                        f"[sc_ambient] '{path}' declares {header[1]} columns but "
                        f"{len(barcodes)} barcodes are listed beside it."
                    )
                continue
            fields = line.split()
            totals[int(fields[1]) - 1] += float(fields[2])
    return totals


def read_cell_taxa(paths, sample, min_umis):
    """({taxid: {barcode: umis}}, names, ranks) for one sample."""
    counts, names, ranks = {}, {}, {}
    for path in paths:
        with open_maybe_gzip(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if sample and row.get("sample") and row["sample"] != sample:
                    continue
                try:
                    value = float(row["count"])
                except (KeyError, ValueError):
                    continue
                if value < min_umis:
                    continue
                taxid = row["taxid"]
                counts.setdefault(taxid, {})[row["barcode"]] = value
                names[taxid] = row.get("name", "")
                ranks[taxid] = row.get("rank", "")
    return counts, names, ranks


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix", required=True, help="STARsolo Solo.out tree (needs raw/ AND filtered/)")
    parser.add_argument("--feature-type", default="Gene")
    parser.add_argument("--cell-taxa", nargs="+", required=True)
    parser.add_argument("--sample", default=None)
    parser.add_argument("--min-umis", type=float, default=2,
                        help="UMIs of a taxon before a droplet counts as carrying it")
    parser.add_argument("--min-empty-umis", type=float, default=1,
                        help="host UMIs an empty droplet needs to be used at all")
    parser.add_argument("--max-empty-umis", type=float, default=0,
                        help="...and at most this many, 0 for no upper bound")
    parser.add_argument("--min-droplets", type=int, default=5,
                        help="a taxon must appear in this many droplets overall to be tested")
    parser.add_argument("--enriched-ratio", type=float, default=2.0,
                        help="prevalence ratio at or above which a taxon is called cell-associated")
    parser.add_argument("--depleted-ratio", type=float, default=0.5,
                        help="...and at or below which it is called cell-depleted")
    parser.add_argument("--p-threshold", type=float, default=0.05, help="on the adjusted p")
    parser.add_argument("--drop-ambient", action="store_true",
                        help="write a drop list of the ambient and cell-depleted taxa. ONLY correct "
                             "when the question is about intracellular microbes")
    parser.add_argument("--prefix", default="sc_ambient")
    args = parser.parse_args()

    filtered_dir = locate(args.matrix, args.feature_type, "filtered")
    raw_dir = locate(args.matrix, args.feature_type, "raw")
    cells = set(read_barcodes(filtered_dir))
    raw_barcodes = read_barcodes(raw_dir)
    if not cells:
        sys.exit(f"[sc_ambient] '{filtered_dir}' lists no cell. Nothing to compare against.")

    totals = host_umis_per_barcode(raw_dir, raw_barcodes)
    empty = []
    for barcode, total in zip(raw_barcodes, totals):
        if barcode in cells:
            continue
        if total < args.min_empty_umis:
            continue
        if args.max_empty_umis and total > args.max_empty_umis:
            continue
        empty.append(barcode)

    if not empty:
        sys.exit(
            f"[sc_ambient] no barcode qualifies as an empty droplet: none outside the "
            f"{len(cells)} cells has between {args.min_empty_umis:g} and "
            f"{args.max_empty_umis or 'unlimited'} host UMIs. Lower --sc_ambient_min_empty_umis, "
            "or check that STARsolo wrote a raw/ matrix at all - with --soloCellFilter None it "
            "writes only one matrix and there is no empty-droplet pool to compare against."
        )

    cell_totals = {b: t for b, t in zip(raw_barcodes, totals) if b in cells}
    empty_set = set(empty)
    empty_totals = {b: t for b, t in zip(raw_barcodes, totals) if b in empty_set}
    host_in_cells = sum(cell_totals.values())
    host_in_empty = sum(empty_totals.values())

    print(
        f"[sc_ambient] {len(cells)} cell(s) carrying {host_in_cells:.0f} host UMIs against "
        f"{len(empty)} empty droplet(s) carrying {host_in_empty:.0f}.",
        file=sys.stderr,
    )
    if host_in_empty <= 0:
        sys.exit("[sc_ambient] the empty droplets hold no host UMIs, so no rate can be formed.")

    counts, names, ranks = read_cell_taxa(args.cell_taxa, args.sample, args.min_umis)
    rows = []
    for taxid, per_barcode in counts.items():
        in_cells = [b for b in per_barcode if b in cells]
        in_empty = [b for b in per_barcode if b in empty_set]
        if len(in_cells) + len(in_empty) < args.min_droplets:
            continue
        umis_cells = sum(per_barcode[b] for b in in_cells)
        umis_empty = sum(per_barcode[b] for b in in_empty)

        prevalence_cells = len(in_cells) / len(cells)
        prevalence_empty = len(in_empty) / len(empty)
        rate_cells = umis_cells / host_in_cells
        rate_empty = umis_empty / host_in_empty

        prevalence_ratio = prevalence_cells / prevalence_empty if prevalence_empty > 0 else None
        rate_ratio = rate_cells / rate_empty if rate_empty > 0 else None

        # One-sided in whichever direction the prevalence points, so the test
        # matches the claim being made about that taxon rather than testing
        # enrichment for taxa that are plainly depleted.
        if prevalence_ratio is not None and prevalence_ratio < 1:
            p = fisher_exact_greater(len(in_empty), len(empty) - len(in_empty),
                                     len(in_cells), len(cells) - len(in_cells))
        else:
            p = fisher_exact_greater(len(in_cells), len(cells) - len(in_cells),
                                     len(in_empty), len(empty) - len(in_empty))

        rows.append({
            "taxid": taxid,
            "name": names.get(taxid, ""),
            "rank": ranks.get(taxid, ""),
            "cells_positive": len(in_cells),
            "cells_total": len(cells),
            "empty_positive": len(in_empty),
            "empty_total": len(empty),
            "umis_in_cells": round(umis_cells, 2),
            "umis_in_empty": round(umis_empty, 2),
            "prevalence_ratio": "inf" if prevalence_ratio is None else round(prevalence_ratio, 4),
            "rate_ratio": "inf" if rate_ratio is None else round(rate_ratio, 4),
            "p": p,
        })

    if not rows:
        sys.exit(
            f"[sc_ambient] no taxon appears in at least {args.min_droplets} droplets at "
            f">= {args.min_umis:g} UMI(s). Nothing to annotate."
        )

    for row, q in zip(rows, benjamini_hochberg([row["p"] for row in rows])):
        row["q"] = q
        ratio = row["prevalence_ratio"]
        significant = q is not None and q < args.p_threshold
        if ratio == "inf" or (ratio >= args.enriched_ratio and significant):
            row["verdict"] = "cell_associated"
        elif ratio <= args.depleted_ratio and significant:
            row["verdict"] = "cell_depleted"
        else:
            # Deliberately not "contaminant". See the module docstring: a taxon
            # as common in the medium as in the cells may be a reagent
            # contaminant or a genuine extracellular organism, and this
            # comparison cannot tell them apart.
            row["verdict"] = "ambient_undecided"

    columns = ["taxid", "name", "rank", "cells_positive", "cells_total", "empty_positive",
               "empty_total", "umis_in_cells", "umis_in_empty", "prevalence_ratio",
               "rate_ratio", "p", "q", "verdict"]
    rows.sort(key=lambda row: -(row["prevalence_ratio"] if row["prevalence_ratio"] != "inf" else 1e9))
    with open(f"{args.prefix}.sc_ambient.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(
                "NA" if row[column] is None else
                (f"{row[column]:.6g}" if column in ("p", "q") else str(row[column]))
                for column in columns
            ) + "\n")

    associated = sum(1 for row in rows if row["verdict"] == "cell_associated")
    depleted = sum(1 for row in rows if row["verdict"] == "cell_depleted")
    ambient = len(rows) - associated - depleted

    if args.drop_ambient:
        with open(f"{args.prefix}.sc_ambient_drop.txt", "w", encoding="utf-8") as handle:
            for row in rows:
                if row["verdict"] in ("ambient_undecided", "cell_depleted"):
                    handle.write(f"{row['taxid']}\n")
        print(
            f"[sc_ambient] --drop-ambient: {ambient + depleted} taxa listed for removal. This is "
            "only the right call if the question is about INTRACELLULAR microbes - an "
            "extracellular organism genuinely present in the tissue is ambient by construction "
            "and is on that list.",
            file=sys.stderr,
        )

    with open(f"{args.prefix}_sc_ambient_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_sc_ambient'",
            "# section_name: 'Cell association vs ambient pool'",
            "# description: 'Whether a taxon is found in cell-containing droplets more often than",
            "#     in empty ones. This is NOT a contamination call. An empty droplet holds reagent",
            "#     contaminants, ambient soup, microbes released by lysis AND genuine extracellular",
            "#     organisms from the tissue, and no droplet-level comparison separates them - so",
            "#     the middle category is reported as undecided rather than as contamination.",
            "#     Only an external measurement of the kit settles that, which is what --decontam",
            "#     blanks provide. Read the two together.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_sc_ambient_plot'",
            "#     title: 'reanaTax: cell association'",
            "#     ylab: 'Taxa'",
            "Sample\tCell-associated\tAmbient (undecided)\tCell-depleted",
            f"{args.sample or 'all cells'}\t{associated}\t{ambient}\t{depleted}",
            "",
        ]))

    print(
        f"[sc_ambient] {len(rows)} taxa: {associated} cell-associated, {ambient} ambient "
        f"(undecided - contaminant or genuinely extracellular), {depleted} cell-depleted.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

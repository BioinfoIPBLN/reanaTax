#!/usr/bin/env python3
"""
sc_pseudobulk.py -- the infected-vs-bystander host response, with the donor as
the unit of replication.

bin/sc_host_de.py tests infected against bystander cells inside one library,
with cells as the units. That is the comparison SAHMI makes, and it is a
statement about that library: thousands of cells from one donor are not
thousands of replicates, and a cohort-level claim built on cell-level p-values
inherits that pseudoreplication (Squair et al., Nat Commun 2021). The standard
answer is pseudobulk: sum each donor's cells per state, and test the sums with
the donor as a blocking factor, so the response has to recur across donors to
count.

This script does the summing-up across libraries. sc_host_de.py --pseudobulk
writes, for every taxon x cell type it tested, the summed raw counts of the
infected cells and of the bystander cells. Here those tables are set side by
side per taxon x cell type, one column per donor and state, and a group is kept
only when enough donors have both states. The test itself is bin/gene_de.R
(DESeq2/edgeR) with `--grouping state --covariates donor`, i.e. a paired design.

The donor is the library (the `sample` the tables carry). A donor sequenced in
several libraries would count several times; merge such libraries upstream.

A group is also kept only if it would survive gene_de.R's own expression filter
(at least two genes with --min-count reads in --min-samples columns): that
script stops with an error on anything less, and one thin group must not take
the run down with it.
"""
import argparse
import csv
import re
import sys


def safe(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tables", nargs="+", required=True, help="*.sc_host_de_pseudobulk.tsv, one per library")
    parser.add_argument("--min-donors", type=int, default=3,
                        help="donors that must have both infected and bystander cells in a group")
    parser.add_argument("--min-count", type=float, default=10, help="as gene_de.R --min-count")
    parser.add_argument("--min-samples", type=int, default=2, help="as gene_de.R --min-samples")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    if args.min_donors < 2:
        sys.exit("[sc_pseudobulk] --min-donors must be at least 2: a paired test needs two donors.")

    # (taxid, cell_type) -> {"taxon": name, "counts": {(donor, state): {gene: count}}, "cells": {...}}
    groups = {}
    for path in args.tables:
        with open(path, encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                key = (row["taxid"], row["cell_type"])
                group = groups.setdefault(key, {"taxon": row["taxon"], "counts": {}, "cells": {}})
                column = (row["sample"], row["state"])
                genes = group["counts"].setdefault(column, {})
                # Duplicated gene names (one symbol, several ids) are summed: at
                # the pseudobulk level they are one more slice of the same reads.
                genes[row["gene"]] = genes.get(row["gene"], 0) + int(row["count"])
                group["cells"][column] = int(row["n_cells"])

    metadata = set()
    summary = []
    for (taxid, cell_type), group in sorted(groups.items()):
        donors = sorted({donor for donor, _state in group["counts"]})
        paired = [donor for donor in donors
                  if (donor, "infected") in group["counts"] and (donor, "bystander") in group["counts"]]
        status = "tested"
        if len(paired) < args.min_donors:
            status = f"skipped: {len(paired)} donor(s) with both states, --min-donors {args.min_donors}"
        columns = [(donor, state) for donor in paired for state in ("infected", "bystander")]
        genes = sorted({gene for column in columns for gene in group["counts"][column]})
        expressed = sum(
            1 for gene in genes
            if sum(1 for column in columns if group["counts"][column].get(gene, 0) >= args.min_count) >= args.min_samples
        )
        if status == "tested" and expressed < 2:
            status = f"skipped: {expressed} gene(s) with --min-count {args.min_count:g} in --min-samples {args.min_samples}"
        summary.append((taxid, group["taxon"], cell_type, len(donors), len(paired), status))
        if status != "tested":
            continue

        names = [f"{donor}.{state}" for donor, state in columns]
        with open(f"{args.prefix}.{safe(taxid)}.{safe(cell_type)}.pseudobulk_counts.tsv", "w", encoding="utf-8") as handle:
            handle.write("gene\t" + "\t".join(names) + "\n")
            for gene in genes:
                handle.write(gene + "\t" + "\t".join(str(group["counts"][column].get(gene, 0)) for column in columns) + "\n")
        for donor, state in columns:
            metadata.add((f"{donor}.{state}", donor, state))

    # One metadata sheet for every group: a column name means the same donor
    # and state in each, and gene_de.R keeps only the columns its matrix has.
    with open(f"{args.prefix}.pseudobulk_metadata.tsv", "w", encoding="utf-8") as handle:
        handle.write("column\tdonor\tstate\n")
        for name, donor, state in sorted(metadata):
            handle.write(f"{name}\t{donor}\t{state}\n")

    with open(f"{args.prefix}.pseudobulk_groups.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxid\ttaxon\tcell_type\tdonors\tdonors_with_both\tstatus\n")
        for entry in summary:
            handle.write("\t".join(str(value) for value in entry) + "\n")

    tested = sum(1 for entry in summary if entry[-1] == "tested")
    print(f"[sc_pseudobulk] {len(summary)} taxon x cell-type group(s); {tested} with at least "
          f"{args.min_donors} donors carrying both states, written for the paired test.", file=sys.stderr)


if __name__ == "__main__":
    main()

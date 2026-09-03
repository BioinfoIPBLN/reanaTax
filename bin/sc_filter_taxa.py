#!/usr/bin/env python3
"""
sc_filter_taxa.py -- apply the evidence filters' drop lists to a cell-by-taxon
matrix.

The drop lists produced by the minimizer filter, the shuffled-read control,
decontam and the host k-mer scan were only ever applied to the COHORT tables -
the combined Kraken2 report and the combined Bracken table. The cell-by-taxon
matrix went to the enrichment test untouched, so on a single-cell run
`--minimizer_filter` changed the abundance tables and nothing else.

That is worst exactly where it matters most. The enrichment test corrects for
multiple testing across every taxon x cell type pair, and on a real cohort the
overwhelming majority of those taxa are background. On the CSI-Microbes plate
the true signal - Fusobacterium in the infected epithelial cells - was the
single strongest result in the table at raw p = 1.6e-05 and still did not
survive BH across 10,901 tests, almost all of them contaminants. Filtering the
matrix does not merely tidy the output; it decides whether a real signal can be
detected at all.

A drop list removes whole taxa and nothing less, because that is the only kind
of claim it makes: a statement about a taxon's evidence across the cohort.
Applying one per cell would invent a verdict nothing computed.

--drop-cells is the exception that proves the rule. The negative-control filter
DOES compute a per-library verdict - it measures how much of each taxon arrives
without a sample, so the same taxon can be signal in one library and carryover
in the next - and hands that verdict over explicitly rather than having it
inferred here. Matched on the barcode first and on the sample second, which is
the same distinction the two single-cell modes draw: in plate mode a library is
a well and the two are the same string, while in droplet mode a library holds
many cells and the control measurement applies to all of them.
"""
import argparse
import csv
import sys


def read_drop_lists(paths):
    """Union of every drop list. One taxid per line, blanks and comments ignored."""
    drop = set()
    for path in paths:
        try:
            with open(path, encoding="utf-8") as handle:
                for line in handle:
                    taxid = line.strip()
                    if taxid and not taxid.startswith("#"):
                        drop.add(taxid.split("\t")[0].strip())
        except OSError as exc:
            print(f"[sc_filter_taxa] cannot read {path}: {exc}", file=sys.stderr)
    return drop


def read_drop_cells(paths):
    """{(taxid, sample)} from the negative-control filter's per-cell verdicts."""
    cells = set()
    for path in paths:
        try:
            with open(path, encoding="utf-8") as handle:
                for line in handle:
                    fields = [f.strip() for f in line.rstrip("\n").split("\t")]
                    if len(fields) < 2 or fields[0] in ("", "taxid") or fields[0].startswith("#"):
                        continue
                    cells.add((fields[0], fields[1]))
        except OSError as exc:
            print(f"[sc_filter_taxa] cannot read {path}: {exc}", file=sys.stderr)
    return cells


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", nargs="+", required=True, help="cell x taxon TSV(s)")
    parser.add_argument("--drop", nargs="*", default=[], help="drop list(s), one taxid per line")
    parser.add_argument("--drop-cells", nargs="*", default=[],
                        help="taxid<TAB>sample verdicts from the negative-control filter")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    drop = read_drop_lists(args.drop)
    drop_cells = read_drop_cells(args.drop_cells)

    header, kept_rows = None, []
    removed_reads, kept_reads = 0.0, 0.0
    removed_taxa, kept_taxa = set(), set()
    removed_names = {}
    removed_cells = 0

    for path in sorted(args.counts):
        with open(path, encoding="utf-8") as handle:
            reader = csv.reader(handle, delimiter="\t")
            for row in reader:
                if not row:
                    continue
                if row[0] == "sample" and header is None:
                    header = row
                    continue
                if row[0] == "sample":
                    continue
                if header is None:
                    sys.exit(f"[sc_filter_taxa] '{path}' has no header row.")
                record = dict(zip(header, row))
                taxid = (record.get("taxid") or "").strip()
                try:
                    count = float(record.get("count") or 0)
                except ValueError:
                    count = 0.0
                cell = (
                    taxid,
                    (record.get("barcode") or "").strip(),
                    (record.get("sample") or "").strip(),
                )
                if taxid in drop or (cell[0], cell[1]) in drop_cells or (cell[0], cell[2]) in drop_cells:
                    if taxid not in drop:
                        removed_cells += 1
                    removed_taxa.add(taxid)
                    removed_names[taxid] = (record.get("name") or "").strip()
                    removed_reads += count
                else:
                    kept_taxa.add(taxid)
                    kept_reads += count
                    kept_rows.append(row)

    if header is None:
        sys.exit("[sc_filter_taxa] no input rows.")

    with open(f"{args.prefix}.cell_taxa.tsv", "w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(kept_rows)

    with open(f"{args.prefix}.cell_taxa_removed.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxid\tname\n")
        for taxid in sorted(removed_taxa):
            handle.write(f"{taxid}\t{removed_names.get(taxid, '')}\n")

    total_reads = kept_reads + removed_reads
    with open(f"{args.prefix}_sc_filter_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_sc_filter'",
            "# section_name: 'Evidence filters on the cell matrix'",
            "# description: 'The drop lists from the minimizer filter, the shuffled-read control,",
            "#     decontam and the host k-mer scan, applied to the cell-by-taxon matrix as well as",
            "#     to the cohort tables. This matters most for the enrichment test, whose",
            "#     multiple-testing correction is levied across every taxon x cell type pair: a",
            "#     background taxon that survives here costs power everywhere.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_sc_filter_plot'",
            "#     title: 'reanaTax: cell matrix filtering'",
            "#     ylab: 'Taxa'",
            "Sample\tKept\tRemoved",
            f"cell matrix\t{len(kept_taxa)}\t{len(removed_taxa)}",
            "",
        ]))

    print(
        f"[sc_filter_taxa] {len(kept_taxa)} taxa kept, {len(removed_taxa)} removed "
        f"({100.0 * removed_reads / total_reads if total_reads else 0:.3f}% of assigned reads); "
        f"drop lists held {len(drop)} taxid(s)"
        + (f", and {len(drop_cells)} per-library verdict(s) removed {removed_cells} further cell(s)"
           if drop_cells else "")
        + ".",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

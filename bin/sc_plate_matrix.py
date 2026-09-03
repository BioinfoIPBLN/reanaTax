#!/usr/bin/env python3
"""
sc_plate_matrix.py -- a cell-by-taxon matrix when every cell is its own library.

CSI-Microbes (Robinson et al., Sci Adv 2024) runs on two kinds of single-cell
data and treats them differently, because they are different experiments
wearing the same name. In droplet data one library holds thousands of cells and
the barcode is what separates them. In plate-based data - Smart-seq2 and its
relatives - each well is sequenced as its own library, there is no barcode, and
the cell is the sample.

That second layout needs none of the barcode machinery in bin/sc_taxa_counts.py
and none of STARsolo. Every cell has already been classified by the ordinary
bulk route, so the cell-by-taxon matrix is just those per-cell reports stacked,
which is what this does. Everything downstream - the presence threshold, the
cell-type enrichment, the co-occurrence test - is unchanged, and reads the same
long-format table bin/sc_taxa_counts.py produces.

Two things are genuinely different and are enforced rather than documented.

  READS, NOT UMIs.  Smart-seq2 has no UMIs, so the presence threshold counts
  reads. A read threshold and a UMI threshold are not interchangeable: a single
  cDNA molecule amplified across a well can contribute hundreds of reads and
  one UMI, so CSI-Microbes' ">=2 UMIs" is a much stronger rule than ">=2 reads".
  The default here is therefore higher, and the column is named for what it
  holds.

  THE SAMPLE IS THE PATIENT, NOT THE CELL.  bin/sc_enrichment.py refuses to run
  on a pooled matrix because pooling cell types across patients inverts
  enrichment (Simpson's paradox, which CSI-Microbes devotes a methods paragraph
  and a counterexample to). In plate-based data the pipeline's `sample` is one
  cell, so a metadata column naming the patient, plate or donor each cell came
  from is REQUIRED - without it every cell would be its own stratum and the
  per-sample-then-combine design would silently collapse to a pooled test.

Reads Kraken2 reports (clade counts, so a genus row includes its species) at
the requested ranks.
"""
import argparse
import csv
import gzip
import os
import sys


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(
        path, encoding="utf-8", errors="replace"
    )


def cell_name(path):
    base = os.path.basename(path)
    for suffix in (".kraken2.report.txt", ".kreport.txt", ".report.txt", ".txt", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def parse_report(path, ranks):
    """[(taxid, rank, name, clade_reads)] for the requested rank codes."""
    out = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("%"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            rank = fields[-3].strip()
            if rank not in ranks:
                continue
            try:
                reads = int(fields[1])
            except ValueError:
                continue
            if reads <= 0:
                continue
            out.append((fields[-2].strip(), rank, fields[-1].strip(), reads))
    return out


def read_metadata(path):
    """{cell: (sample, cell_type)}. The first column is the cell id."""
    mapping = {}
    with open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            sys.exit(f"[sc_plate_matrix] '{path}' is empty.")
        lookup = {name.strip().lower(): name for name in reader.fieldnames}
        cell_key = lookup.get("barcode") or lookup.get("cell") or lookup.get("cell_id") or reader.fieldnames[0]
        sample_key = lookup.get("sample") or lookup.get("patient") or lookup.get("donor") or lookup.get("plate")
        type_key = lookup.get("cell_type") or lookup.get("celltype")
        if not sample_key:
            sys.exit(
                f"[sc_plate_matrix] '{path}' has no sample/patient/donor/plate column. In "
                "plate-based data the pipeline's sample IS one cell, so without a column naming "
                "the patient each cell came from every cell becomes its own stratum and the "
                "per-sample-then-combine enrichment test collapses to a pooled one - the exact "
                "failure CSI-Microbes' methods warn about. Columns present: "
                + ", ".join(reader.fieldnames)
            )
        for row in reader:
            cell = (row.get(cell_key) or "").strip()
            if not cell:
                continue
            mapping[cell] = (
                (row.get(sample_key) or "").strip(),
                (row.get(type_key) or "").strip() if type_key else "",
            )
    return mapping


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reports", nargs="+", required=True, help="one Kraken2 report per cell")
    parser.add_argument("--metadata", required=True, help="TSV: cell id, sample/patient[, cell_type]")
    parser.add_argument("--ranks", default="S,G", help="comma-separated Kraken rank codes to tabulate")
    parser.add_argument("--min-reads", type=int, default=10,
                        help="reads of a taxon before a cell counts as carrying it")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    ranks = {rank.strip() for rank in args.ranks.split(",") if rank.strip()}
    metadata = read_metadata(args.metadata)

    cells = [cell_name(path) for path in args.reports]
    unknown = [cell for cell in cells if cell not in metadata]
    if len(unknown) == len(cells):
        sys.exit(
            "[sc_plate_matrix] not one of the classified cells appears in the metadata. The "
            "first metadata column must hold the sample ids the pipeline used, one per cell. "
            f"Cells: {', '.join(sorted(cells)[:5])}{'...' if len(cells) > 5 else ''}"
        )
    if unknown:
        print(f"[sc_plate_matrix] {len(unknown)} cell(s) not in the metadata, excluded: "
              f"{', '.join(sorted(unknown)[:10])}{'...' if len(unknown) > 10 else ''}", file=sys.stderr)

    rows = []
    kept_cells = set()
    dropped_pairs = 0
    for path, cell in zip(args.reports, cells):
        if cell not in metadata:
            continue
        sample, _cell_type = metadata[cell]
        if not sample:
            sys.exit(f"[sc_plate_matrix] cell '{cell}' has an empty sample/patient value.")
        for taxid, rank, name, reads in parse_report(path, ranks):
            if reads < args.min_reads:
                dropped_pairs += 1
                continue
            # Same long format as bin/sc_taxa_counts.py, so bin/sc_enrichment.py
            # and everything else downstream cannot tell the two routes apart.
            # `barcode` holds the cell id; `sample` holds the patient.
            rows.append([sample, cell, taxid, rank, name, reads])
            kept_cells.add(cell)

    with open(f"{args.prefix}.cell_taxa.tsv", "w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "barcode", "taxid", "rank", "name", "count"])
        writer.writerows(rows)

    samples = {row[0] for row in rows}
    taxa = {row[2] for row in rows}
    with open(f"{args.prefix}.cell_taxa_summary.tsv", "w", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["cells_classified", len(cells)])
        writer.writerow(["cells_in_metadata", len(cells) - len(unknown)])
        writer.writerow(["cells_carrying_a_taxon", len(kept_cells)])
        writer.writerow(["patients", len(samples)])
        writer.writerow(["taxa", len(taxa)])
        writer.writerow(["cell_taxon_pairs", len(rows)])
        writer.writerow(["pairs_below_min_reads", dropped_pairs])
        writer.writerow(["min_reads", args.min_reads])

    if len(samples) < 2:
        print(
            f"[sc_plate_matrix] WARNING: every cell maps to the single patient "
            f"'{next(iter(samples), '')}'. The enrichment step stratifies by patient and combines "
            "across them, so with one stratum it reduces to a single Fisher test - valid, but the "
            "protection against Simpson's paradox is not doing anything.",
            file=sys.stderr,
        )

    print(
        f"[sc_plate_matrix] {len(kept_cells)}/{len(cells)} cell(s) carry at least one taxon at "
        f">= {args.min_reads} read(s); {len(taxa)} taxa across {len(samples)} patient(s).",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

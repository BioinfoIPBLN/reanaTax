#!/usr/bin/env python3
"""
sc_host_de.py -- what the host cell does when it is carrying something.

This is the half of SAHMI that reanaTax has been missing. SAHMI's published
pipeline does not stop at a denoised cell-by-taxon matrix: the matrix exists so
that the barcodes carrying a taxon can be set against the barcodes that do not,
and the HOST transcriptome compared between them. That comparison is the result
in the paper - bacteria-positive tumour cells with an altered inflammatory
programme - and everything upstream is machinery for making it trustworthy.

The comparison here is the one SAHMI makes: infected cells against bystander
cells, Wilcoxon rank-sum on log-normalised counts, which is what Seurat's
FindMarkers does by default and what SAHMI calls.

WITHIN A CELL TYPE, always. This is the whole design and it is not negotiable.
Infection is not distributed at random over cell types - measuring that is the
entire point of the enrichment step in bin/sc_enrichment.py - so a pooled
infected-vs-uninfected test recovers the difference between the cell types that
happen to be infected and the cell types that do not, and reports it as a
response to infection. The two are indistinguishable in a pooled test. So
--cells is required, and running without it needs --force-pooled and stamps
every row with `ALL_POOLED` so nobody can read the output as a within-type
result by accident.

What "infected" means here is a threshold on the cell-by-taxon matrix
(--min-umis, the same rule CSI-Microbes uses), applied AFTER whatever denoising
ran upstream. What "bystander" means is: a cell of the same type, in the same
sample, carrying none of that taxon. Not "carrying no microbe at all" - a cell
positive for a different taxon is still a valid bystander for this one, and
excluding it would silently restrict the comparison to the cleanest cells in
the library.

Normalisation is Seurat's LogNormalize: counts per cell scaled to 10,000, then
log1p. Fold change is reported on the linear scale, as log2 of the ratio of
mean expm1 values - again Seurat's definition, so the numbers are comparable to
what anyone would get by loading the same matrix into Seurat.

Reads STARsolo's MatrixMarket output directly. No Seurat, no scanpy, no
AnnData: the single-cell modules here run in a plain Python container and
adding a whole single-cell framework to run one rank-sum test would be a
strange trade - the same reasoning that produced bin/sahmi_stats.py.
"""
import argparse
import csv
import gzip
import math
import os
import sys
from array import array

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import benjamini_hochberg, ranksums

SCALE_FACTOR = 10000.0


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


def locate_matrix(root, feature_type, filtering):
    """The MatrixMarket directory inside a STARsolo Solo.out tree, or root itself."""
    if find_file(root, "matrix.mtx"):
        return root
    preferred = os.path.join(root, feature_type, filtering)
    if find_file(preferred, "matrix.mtx"):
        return preferred
    for base, _dirs, files in os.walk(root):
        if any(name.startswith("matrix.mtx") for name in files):
            if os.path.basename(base) == filtering:
                return base
    sys.exit(
        f"[sc_host_de] no matrix.mtx under '{root}'. Expected a STARsolo Solo.out tree holding "
        f"{feature_type}/{filtering}/matrix.mtx, or a directory that is one."
    )


def read_lines(path, column=0):
    values = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            values.append(fields[column] if column < len(fields) else fields[0])
    return values


def read_matrix(directory):
    """(genes, barcodes, gene_index, cell_index, values, cell_totals).

    The three triplet arrays are sorted by gene, so every gene's nonzero entries
    are one contiguous slice - which is what makes a per-gene test over tens of
    thousands of genes affordable without a sparse-matrix library.
    """
    matrix_path = find_file(directory, "matrix.mtx")
    features_path = find_file(directory, "features.tsv", "genes.tsv")
    barcodes_path = find_file(directory, "barcodes.tsv")
    for label, path in (("matrix.mtx", matrix_path), ("features.tsv", features_path),
                        ("barcodes.tsv", barcodes_path)):
        if path is None:
            sys.exit(f"[sc_host_de] '{directory}' has no {label}.")

    # STARsolo's features.tsv is id, name, type; the gene NAME is column 2 when
    # it is there, because a table of ENSG ids helps nobody read the result.
    genes = []
    with open_maybe_gzip(features_path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if not fields or not fields[0]:
                continue
            genes.append(fields[1] if len(fields) > 1 and fields[1] else fields[0])
    barcodes = read_lines(barcodes_path)

    gene_index, cell_index, values = array("i"), array("i"), array("f")
    with open_maybe_gzip(matrix_path) as handle:
        header = None
        for line in handle:
            if line.startswith("%"):
                continue
            if header is None:
                header = line.split()
                rows, columns = int(header[0]), int(header[1])
                if rows != len(genes) or columns != len(barcodes):
                    sys.exit(
                        f"[sc_host_de] matrix.mtx declares {rows}x{columns} but there are "
                        f"{len(genes)} feature(s) and {len(barcodes)} barcode(s)."
                    )
                continue
            row, column, value = line.split()
            gene_index.append(int(row) - 1)
            cell_index.append(int(column) - 1)
            values.append(float(value))

    cell_totals = [0.0] * len(barcodes)
    for column, value in zip(cell_index, values):
        cell_totals[column] += value

    # Counting sort into gene order. The MatrixMarket spec does not promise an
    # ordering and STARsolo's is not documented, so it is imposed rather than
    # assumed - reading a gene's row from the wrong slice would produce a table
    # of confident nonsense.
    counts = [0] * (len(genes) + 1)
    for gene in gene_index:
        counts[gene + 1] += 1
    for position in range(1, len(counts)):
        counts[position] += counts[position - 1]
    offsets = list(counts)
    sorted_cells = array("i", [0]) * len(gene_index)
    sorted_values = array("f", [0.0]) * len(values)
    cursor = list(counts[:-1])
    for gene, column, value in zip(gene_index, cell_index, values):
        position = cursor[gene]
        sorted_cells[position] = column
        sorted_values[position] = value
        cursor[gene] = position + 1

    return genes, barcodes, offsets, sorted_cells, sorted_values, cell_totals


def read_cell_taxa(paths, sample, min_umis):
    """{taxid: {barcode}} for one sample, plus the taxon's name and rank."""
    positive, names, ranks = {}, {}, {}
    for path in paths:
        with open_maybe_gzip(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if sample and row.get("sample") and row["sample"] != sample:
                    continue
                try:
                    count = float(row["count"])
                except (KeyError, ValueError):
                    continue
                if count < min_umis:
                    continue
                taxid = row["taxid"]
                positive.setdefault(taxid, set()).add(row["barcode"])
                names[taxid] = row.get("name", "")
                ranks[taxid] = row.get("rank", "")
    return positive, names, ranks


def read_cell_types(path, sample):
    """{barcode: cell_type} for one sample."""
    types = {}
    with open_maybe_gzip(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            keys = {key.strip().lower(): key for key in row}
            barcode = row.get(keys.get("barcode", ""), "").strip()
            if not barcode:
                continue
            if sample and "sample" in keys:
                if row[keys["sample"]].strip() != sample:
                    continue
            cell_type = row.get(keys.get("cell_type", keys.get("celltype", "")), "").strip()
            if cell_type:
                types[barcode] = cell_type
    return types


def normalised(gene, offsets, cells, values, factors, position, size):
    """Log-normalised expression of one gene over a group, zeros included.

    `position` maps a matrix column to its index in the group, and is built once
    per group rather than once per gene - with tens of thousands of genes that
    is the difference between one pass over the members and tens of thousands.
    The returned list is in group order, so two groups are directly comparable.
    """
    out = [0.0] * size
    for offset in range(offsets[gene], offsets[gene + 1]):
        index = position.get(cells[offset])
        if index is not None:
            out[index] = math.log1p(values[offset] * factors[cells[offset]])
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--matrix", required=True, help="STARsolo Solo.out tree, or a MatrixMarket directory")
    parser.add_argument("--feature-type", default="Gene", help="Gene, GeneFull, ... (inside Solo.out)")
    parser.add_argument("--filtering", default="filtered", choices=("filtered", "raw"))
    parser.add_argument("--cell-taxa", nargs="+", required=True, help="cell_taxa.tsv from SCTAXA_COUNTS")
    parser.add_argument("--sample", default=None, help="restrict every input to this sample id")
    parser.add_argument("--cells", default=None, help="TSV with barcode, cell_type[, sample]")
    parser.add_argument("--force-pooled", action="store_true",
                        help="run without --cells, pooling every cell type (see the module docstring)")
    parser.add_argument("--min-umis", type=float, default=2, help="UMIs of a taxon before a cell counts as infected")
    parser.add_argument("--min-cells", type=int, default=10, help="infected cells a group needs to be tested")
    parser.add_argument("--min-pct", type=float, default=0.1, help="a gene must be detected in this fraction of one side")
    parser.add_argument("--logfc-threshold", type=float, default=0.25, help="minimum |log2 fold change| to test")
    parser.add_argument("--top-taxa", type=int, default=20, help="most-infecting taxa to test (0 = all)")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--prefix", default="sc_host_de")
    args = parser.parse_args()

    if not args.cells and not args.force_pooled:
        sys.exit(
            "[sc_host_de] REFUSING: --cells was not given.\n"
            "  Infected cells are not a random sample of the library - which cell types carry a\n"
            "  taxon is itself a result, produced by the enrichment step. Pooling every cell type\n"
            "  into one infected-vs-uninfected test therefore recovers the difference BETWEEN\n"
            "  cell types and reports it as a response TO infection, and nothing in the output\n"
            "  distinguishes the two. Supply cell-type annotations, or pass --force-pooled to\n"
            "  accept a pooled comparison; every row is then stamped ALL_POOLED."
        )

    matrix_dir = locate_matrix(args.matrix, args.feature_type, args.filtering)
    genes, barcodes, offsets, cells, values, cell_totals = read_matrix(matrix_dir)
    print(f"[sc_host_de] {len(genes)} feature(s) x {len(barcodes)} barcode(s) from '{matrix_dir}'.",
          file=sys.stderr)

    factors = [SCALE_FACTOR / total if total > 0 else 0.0 for total in cell_totals]
    column_of = {barcode: index for index, barcode in enumerate(barcodes)}

    positive, names, ranks = read_cell_taxa(args.cell_taxa, args.sample, args.min_umis)
    if not positive:
        sys.exit(
            f"[sc_host_de] no cell carries any taxon at >= {args.min_umis:g} UMI(s)"
            + (f" in sample '{args.sample}'" if args.sample else "")
            + ". Nothing to compare."
        )

    cell_types = read_cell_types(args.cells, args.sample) if args.cells else {}
    if args.cells and not cell_types:
        sys.exit(f"[sc_host_de] '{args.cells}' names no cell of "
                 f"{'sample ' + args.sample if args.sample else 'this run'}.")

    # Only barcodes present in BOTH the matrix and the annotation are usable.
    # The filtered matrix has already had empty droplets removed, so a barcode
    # carrying a taxon but absent from it is ambient signal on a non-cell and is
    # correctly excluded here.
    universe = [barcode for barcode in barcodes if barcode in column_of]
    if cell_types:
        universe = [barcode for barcode in universe if barcode in cell_types]
    by_type = {}
    for barcode in universe:
        by_type.setdefault(cell_types.get(barcode, "ALL_POOLED"), []).append(barcode)

    ordering = sorted(positive, key=lambda taxid: -len(positive[taxid] & set(universe)))
    if args.top_taxa:
        ordering = ordering[: args.top_taxa]

    rows = []
    tested_groups = []
    for taxid in ordering:
        carriers = positive[taxid]
        for cell_type, members in sorted(by_type.items()):
            infected = [column_of[b] for b in members if b in carriers]
            bystander = [column_of[b] for b in members if b not in carriers]
            if len(infected) < args.min_cells or len(bystander) < args.min_cells:
                continue
            tested_groups.append((taxid, cell_type, len(infected), len(bystander)))
            group = infected + bystander
            group_set = set(group)
            group_position = {column: index for index, column in enumerate(group)}

            # Detection rates first: they need only the nonzero entries, and
            # they throw out most of the transcriptome before the expensive
            # per-gene vectors are built for anything.
            infected_columns = set(infected)
            keep = []
            for gene in range(len(genes)):
                hits_infected = 0
                hits_bystander = 0
                for offset in range(offsets[gene], offsets[gene + 1]):
                    column = cells[offset]
                    if column in infected_columns:
                        hits_infected += 1
                    elif column in group_set:
                        hits_bystander += 1
                pct_infected = hits_infected / len(infected)
                pct_bystander = hits_bystander / len(bystander)
                if max(pct_infected, pct_bystander) >= args.min_pct:
                    keep.append((gene, pct_infected, pct_bystander))

            group_rows = []
            for gene, pct_infected, pct_bystander in keep:
                expression = normalised(gene, offsets, cells, values, factors,
                                        group_position, len(group))
                left = expression[: len(infected)]
                right = expression[len(infected):]
                mean_infected = sum(math.expm1(v) for v in left) / len(left)
                mean_bystander = sum(math.expm1(v) for v in right) / len(right)
                log2fc = math.log2((mean_infected + 1e-9) / (mean_bystander + 1e-9))
                if abs(log2fc) < args.logfc_threshold:
                    continue
                _z, p = ranksums(left, right)
                group_rows.append({
                    "taxid": taxid,
                    "taxon": names.get(taxid, ""),
                    "rank": ranks.get(taxid, ""),
                    "cell_type": cell_type,
                    "gene": genes[gene],
                    "n_infected": len(infected),
                    "n_bystander": len(bystander),
                    "pct_infected": round(pct_infected, 4),
                    "pct_bystander": round(pct_bystander, 4),
                    "mean_infected": round(mean_infected, 4),
                    "mean_bystander": round(mean_bystander, 4),
                    "log2fc": round(log2fc, 4),
                    "p": p,
                })
            # Adjusted WITHIN the group. Each (taxon, cell type) is a separate
            # question asked of a separate set of cells, and pooling the genes
            # of every group into one family would let a group with thousands of
            # tested genes set the threshold for a group with a hundred.
            for row, q in zip(group_rows, benjamini_hochberg([row["p"] for row in group_rows])):
                row["q"] = q
            rows.extend(group_rows)

    columns = ["taxid", "taxon", "rank", "cell_type", "gene", "n_infected", "n_bystander",
               "pct_infected", "pct_bystander", "mean_infected", "mean_bystander", "log2fc", "p", "q"]
    rows.sort(key=lambda row: (row["q"] if row["q"] is not None else 1.0, -abs(row["log2fc"])))
    with open(f"{args.prefix}.sc_host_de.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(
                "NA" if row[column] is None else
                (f"{row[column]:.6g}" if column in ("p", "q") else str(row[column]))
                for column in columns
            ) + "\n")

    with open(f"{args.prefix}.sc_host_de_groups.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxid\ttaxon\tcell_type\tn_infected\tn_bystander\tgenes_significant\n")
        for taxid, cell_type, n_infected, n_bystander in tested_groups:
            hits = sum(
                1 for row in rows
                if row["taxid"] == taxid and row["cell_type"] == cell_type
                and row["q"] is not None and row["q"] < args.p_threshold
            )
            handle.write(f"{taxid}\t{names.get(taxid, '')}\t{cell_type}\t{n_infected}\t"
                         f"{n_bystander}\t{hits}\n")

    significant = sum(1 for row in rows if row["q"] is not None and row["q"] < args.p_threshold)
    with open(f"{args.prefix}_sc_host_de_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_sc_host_de'",
            "# section_name: 'Host response in infected cells'",
            "# description: 'Host genes differing between cells carrying a taxon and bystander",
            "#     cells OF THE SAME TYPE in the same library - the comparison SAHMI makes once",
            "#     its cell-by-taxon matrix is denoised. Wilcoxon rank-sum on log-normalised",
            "#     counts, adjusted within each taxon x cell type. Cells are matched on type",
            "#     because which cell types carry a taxon is itself a result, so a pooled test",
            "#     would report the difference between cell types as a response to infection.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_sc_host_de_plot'",
            "#     title: 'reanaTax: infected vs bystander'",
            "#     ylab: 'Genes'",
            "Sample\tSignificant\tTested",
            f"{args.sample or 'all cells'}\t{significant}\t{len(rows) - significant}",
            "",
        ]))

    print(
        f"[sc_host_de] {len(tested_groups)} taxon x cell-type group(s) tested; "
        f"{significant}/{len(rows)} gene call(s) at q < {args.p_threshold:g}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

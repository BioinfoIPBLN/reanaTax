#!/usr/bin/env python3
"""
sc_enrichment.py -- which host cell types carry which microbes?

This is the question the whole single-cell branch exists to answer, and the
statistics are CSI-Microbes' (Robinson et al., Sci Adv 2024) rather than
SAHMI's, for three reasons the paper argues explicitly:

  Fisher, not chi-square   "for the chi-square approximation to be valid, the
      expected frequency should be at least 5", and that fails for most
      sample x cell type x taxon combinations. Llorens-Rico et al. use a pooled
      chi-square; on sparse data Fisher is the correct exact alternative.

  Presence, not abundance  the matrix is >90% zeros, so the useful question is
      what FRACTION OF CELLS carry a taxon, not how much of it they carry.
      That is also why ALDEx2/ANCOM-BC2, which model compositional abundance
      across samples, are the wrong tools here and are not used.

  Per-sample, then combined  never pooled. The paper devotes a methods
      paragraph and a counterexample to this: two samples where cell type 2 is
      clearly enriched within each can pool to make cell type 1 look enriched,
      purely from differing cell-type composition. This is Simpson's paradox,
      and it is exactly the failure mode a pipeline that concatenates public
      datasets will hit. So --cells MUST carry a sample column, and this script
      refuses to run on a pooled matrix rather than silently producing it.

Effect size is log2(observed / expected), where expected infected cells for a
type in a sample is (cells of that type / cells in sample) x infected cells in
sample. Across the cohort the observed counts are summed AND the expected counts
are summed, and the ratio taken only at the end - never an average of ratios,
which would weight a two-cell sample like a two-thousand-cell one.
"""
import argparse
import csv
import itertools
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import (
    benjamini_hochberg, fisher_exact_greater, hypergeom_sf, ranksums, stouffer,
)


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            yield row


def column(row, *names):
    for name in names:
        for key in row:
            if key.strip().lower() == name:
                return row[key].strip()
    return None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--counts", nargs="+", required=True,
                        help="cell x taxon tables from SCTAXA_COUNTS")
    parser.add_argument("--cells", required=True,
                        help="TSV of barcode, sample, cell_type [, host_umis]")
    parser.add_argument("--enrichment", required=True)
    parser.add_argument("--cooccurrence")
    parser.add_argument("--doublet-test")
    parser.add_argument("--mqc")
    parser.add_argument("--min-cells", type=int, default=10,
                        help="a taxon must infect this many cells to enter the co-occurrence test")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    args = parser.parse_args()

    # ---- host cells -------------------------------------------------------
    #
    # The aliases must stay in step with sc_plate_matrix.py: the same
    # --sc_cell_metadata file feeds both, and a name accepted by one and not the
    # other produced a file that was valid for the plate matrix and invisible
    # here. The columns are also resolved BEFORE the rows are walked, so a
    # naming mismatch reports which column is missing instead of an empty result.
    rows = list(read_tsv(args.cells))
    present = sorted({key.strip() for row in rows for key in row})
    def resolve(*names):
        lowered = {key.strip().lower() for row in rows[:1] for key in row}
        return next((n for n in names if n in lowered), None)

    barcode_key = resolve("barcode", "cb", "cell", "cell_barcode", "cell_id")
    sample_key = resolve("sample", "sample_id", "sampleid", "patient", "donor", "plate")
    type_key = resolve("cell_type", "celltype", "cluster", "annotation")
    if rows and not barcode_key:
        raise SystemExit(
            "sc_enrichment: --cells has no cell-id column (barcode/cb/cell/cell_barcode/"
            f"cell_id). Columns present: {', '.join(present)}"
        )
    if rows and not sample_key:
        raise SystemExit(
            "sc_enrichment: --cells has no `sample` column (sample/sample_id/patient/donor/"
            "plate). Every test here is computed within a sample and combined afterwards, "
            "because pooling cells across samples with different cell-type compositions "
            "inverts enrichment results (Simpson's paradox; see CSI-Microbes, Robinson et al. "
            f"2024). Columns present: {', '.join(present)}"
        )
    if rows and not type_key:
        raise SystemExit(
            "sc_enrichment: --cells has no `cell_type` column (cell_type/celltype/cluster/"
            f"annotation). Columns present: {', '.join(present)}"
        )

    cells = {}
    for row in rows:
        barcode = column(row, barcode_key)
        sample = column(row, sample_key)
        cell_type = column(row, type_key)
        if not barcode or not cell_type or not sample:
            continue
        host_umis = column(row, "host_umis", "numi", "n_umi", "ncount_rna")
        cells[(sample, barcode)] = {
            "cell_type": cell_type,
            "host_umis": float(host_umis) if host_umis and host_umis.replace(".", "", 1).isdigit() else None,
        }
    if not cells:
        raise SystemExit(f"sc_enrichment: no usable rows in {args.cells}.")

    # ---- infected cells ---------------------------------------------------
    # The counts table has already had the presence rule (--sc_min_umis)
    # applied, so every row here IS a call of "this taxon is in this cell".
    infected = {}
    names = {}
    for path in args.counts:
        for row in read_tsv(path):
            sample = column(row, "sample")
            barcode = column(row, "barcode")
            taxid = column(row, "taxid")
            if not (sample and barcode and taxid):
                continue
            if (sample, barcode) not in cells:
                continue
            infected.setdefault(taxid, set()).add((sample, barcode))
            names.setdefault(taxid, column(row, "name") or "")

    if not infected:
        raise SystemExit(
            "sc_enrichment: no cell in the counts table matches a barcode in --cells. "
            "STARsolo barcodes usually carry a '-1' suffix; check the two files agree."
        )

    by_sample = {}
    for (sample, barcode), meta in cells.items():
        by_sample.setdefault(sample, []).append((barcode, meta["cell_type"]))

    # ---- cell-type enrichment --------------------------------------------
    cell_types = sorted({meta["cell_type"] for meta in cells.values()})
    rows, pvalues = [], []

    for taxid in sorted(infected, key=lambda t: -len(infected[t])):
        positives = infected[taxid]
        for cell_type in cell_types:
            per_sample_p, weights = [], []
            observed_total = expected_total = 0.0
            samples_used = 0
            for sample, members in by_sample.items():
                n_sample = len(members)
                n_type = sum(1 for _bc, ct in members if ct == cell_type)
                n_infected = sum(1 for bc, _ct in members if (sample, bc) in positives)
                if not n_type or not n_infected or n_type == n_sample:
                    continue
                a = sum(1 for bc, ct in members if ct == cell_type and (sample, bc) in positives)
                b = n_type - a
                c = n_infected - a
                d = n_sample - n_type - c
                pvalue = fisher_exact_greater(a, b, c, d)
                expected = n_type * (n_infected / n_sample)
                per_sample_p.append(pvalue)
                weights.append(expected)
                observed_total += a
                expected_total += expected
                samples_used += 1
            if samples_used == 0 or expected_total <= 0:
                continue
            _z, combined = stouffer(per_sample_p, weights)
            if combined is None:
                continue
            # Sum numerators, sum denominators, divide last.
            log2fc = math.log2(observed_total / expected_total) if observed_total > 0 else float("-inf")
            rows.append({
                "taxid": taxid, "name": names.get(taxid, ""), "cell_type": cell_type,
                "samples": samples_used, "observed": int(observed_total),
                "expected": expected_total, "log2fc": log2fc, "p": combined,
            })
            pvalues.append(combined)

    for row, qvalue in zip(rows, benjamini_hochberg(pvalues)):
        row["p_adj"] = qvalue

    rows.sort(key=lambda r: (r["p_adj"], -r["observed"]))
    with open(args.enrichment, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "name", "cell_type", "samples", "observed_cells",
                         "expected_cells", "log2fc", "p_stouffer", "p_adj", "significant"])
        for row in rows:
            writer.writerow([
                row["taxid"], row["name"], row["cell_type"], row["samples"],
                row["observed"], f"{row['expected']:.3f}",
                "-Inf" if row["log2fc"] == float("-inf") else f"{row['log2fc']:.4f}",
                f"{row['p']:.6g}", f"{row['p_adj']:.6g}",
                "yes" if row["p_adj"] <= args.p_threshold else "no",
            ])

    # ---- co-occurrence ----------------------------------------------------
    pairs = []
    if args.cooccurrence:
        eligible = [t for t in infected if len(infected[t]) >= args.min_cells]
        pair_p = []
        for first, second in itertools.combinations(sorted(eligible), 2):
            per_sample_p, weights = [], []
            observed_total = expected_total = 0.0
            for sample, members in by_sample.items():
                n_sample = len(members)
                in_first = {bc for bc, _ct in members if (sample, bc) in infected[first]}
                in_second = {bc for bc, _ct in members if (sample, bc) in infected[second]}
                if not in_first or not in_second:
                    continue
                both = len(in_first & in_second)
                expected = len(in_first) * len(in_second) / n_sample
                per_sample_p.append(hypergeom_sf(both, n_sample, len(in_first), len(in_second)))
                weights.append(expected)
                observed_total += both
                expected_total += expected
            if not per_sample_p or expected_total <= 0:
                continue
            _z, combined = stouffer(per_sample_p, weights)
            if combined is None:
                continue
            pairs.append({
                "a": first, "b": second, "observed": int(observed_total),
                "expected": expected_total, "p": combined,
            })
            pair_p.append(combined)
        for pair, qvalue in zip(pairs, benjamini_hochberg(pair_p)):
            pair["p_adj"] = qvalue
        pairs.sort(key=lambda p: p["p_adj"])
        with open(args.cooccurrence, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["taxid_a", "name_a", "taxid_b", "name_b", "observed_cells",
                             "expected_cells", "p_stouffer", "p_adj", "significant"])
            for pair in pairs:
                writer.writerow([
                    pair["a"], names.get(pair["a"], ""), pair["b"], names.get(pair["b"], ""),
                    pair["observed"], f"{pair['expected']:.3f}",
                    f"{pair['p']:.6g}", f"{pair['p_adj']:.6g}",
                    "yes" if pair["p_adj"] <= args.p_threshold else "no",
                ])

    # ---- doublet control --------------------------------------------------
    # Mandatory companion to any co-occurrence claim: cells carrying two taxa
    # may simply be undetected doublets, which would explain the co-occurrence
    # entirely. The check is whether they carry more HOST UMIs than
    # singly-infected cells.
    if args.doublet_test:
        per_cell = {}
        for taxid, positives in infected.items():
            for key in positives:
                per_cell[key] = per_cell.get(key, 0) + 1
        multi = [cells[k]["host_umis"] for k, n in per_cell.items() if n > 1 and cells[k]["host_umis"] is not None]
        single = [cells[k]["host_umis"] for k, n in per_cell.items() if n == 1 and cells[k]["host_umis"] is not None]
        with open(args.doublet_test, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["metric", "value"])
            writer.writerow(["cells_with_multiple_taxa", len([1 for n in per_cell.values() if n > 1])])
            writer.writerow(["cells_with_one_taxon", len([1 for n in per_cell.values() if n == 1])])
            if len(multi) >= 3 and len(single) >= 3:
                z, pvalue = ranksums(multi, single)
                median_multi = sorted(multi)[len(multi) // 2]
                median_single = sorted(single)[len(single) // 2]
                writer.writerow(["median_host_umis_multi", f"{median_multi:.1f}"])
                writer.writerow(["median_host_umis_single", f"{median_single:.1f}"])
                writer.writerow(["ranksum_z", "NA" if z is None else f"{z:.4f}"])
                writer.writerow(["ranksum_p", "NA" if pvalue is None else f"{pvalue:.6g}"])
                verdict = "co-occurrence may be driven by doublets" \
                    if pvalue is not None and pvalue < 0.05 and median_multi > median_single \
                    else "no doublet signal"
                writer.writerow(["verdict", verdict])
            else:
                writer.writerow(["verdict", "not enough cells, or no host_umis column, to test"])

    if args.mqc:
        significant = sum(1 for row in rows if row["p_adj"] <= args.p_threshold)
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_sc_enrichment'\n"
                "# section_name: 'Cell-type enrichment'\n"
                "# description: 'Host cell types carrying each taxon more often than their share of\n"
                "#     cells predicts. One-sided Fisher exact per sample, combined across samples by\n"
                "#     Stouffer Z weighted by expected infected cells. Computed WITHIN samples and\n"
                "#     combined afterwards, never pooled - pooling inverts enrichment when samples\n"
                "#     differ in cell-type composition.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_sc_enrichment_plot'\n"
                "#     title: 'reanaTax: cell-type enrichment'\n"
                "#     ylab: 'taxon x cell-type tests'\n"
                "Sample\tSignificant\tNot significant\n"
                f"all\t{significant}\t{len(rows) - significant}\n"
            )

    print(
        f"[sc_enrichment] {len(infected)} taxa x {len(cell_types)} cell types over "
        f"{len(by_sample)} sample(s): {len(rows)} tests, "
        f"{sum(1 for r in rows if r['p_adj'] <= args.p_threshold)} significant"
        + (f"; {len(pairs)} co-occurrence pairs tested." if args.cooccurrence else "."),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

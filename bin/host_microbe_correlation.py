#!/usr/bin/env python3
"""Correlate host gene expression against microbial abundance across samples.

This is the third layer Monteleone et al. describe and the one thing that makes
--quantify_host more than a QC by-product: host expression and microbial
composition measured in the SAME library, so a correlation between them is not
confounded by sample handling the way two separate assays would be.

Both sides are made scale-free before correlating. Host counts go to CPM, then
log; microbial abundances are already relative. Spearman is used rather than
Pearson because neither side is close to normal and abundance tables are
dominated by a few taxa.

THE FEASIBILITY GUARD IS THE POINT OF THIS SCRIPT.

With n samples, the smallest two-sided Spearman p-value attainable is 2/n!, no
matter how clean the data. Testing g genes against t taxa is g*t tests, and
under Benjamini-Hochberg the most significant of them must clear alpha/(g*t) to
be called. For a five-sample cohort the floor is 2/120 = 0.0167, so at alpha=0.05
no more than 3 pairs could EVER be significant - against the millions a
gene-by-taxon grid contains. Run blind, such an analysis returns an empty table
that reads like a negative result and is nothing of the kind.

So the number of testable pairs is computed up front and compared against the
floor. If the grid cannot produce a single significant pair the run is refused,
with the arithmetic, rather than reported as "no associations found".
"""

import argparse
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sahmi_stats import spearman, benjamini_hochberg  # noqa: E402


ANNOTATION_COLUMNS = {"name", "taxonomy_id", "taxonomy_lvl", "taxid", "rank", "lvl_type", "gene_id"}


def read_host(path):
    """gene x sample counts as written by bin/merge_featurecounts.py: the first
    column is the gene id and every other column is one sample."""
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:]
        rows = {}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            values = []
            for cell in fields[1:]:
                try:
                    values.append(float(cell))
                except ValueError:
                    values.append(0.0)
            rows[fields[0]] = values
    return samples, rows


def read_microbial(path):
    """The combined Bracken/Kraken table, which does NOT have one column per
    sample: Bracken writes `<sample>_num` and `<sample>_frac` side by side and
    combine_kreports writes `<sample>_all` and `<sample>_lvl`. The same
    convention bin/filter_abundance.py follows is used here - prefer the
    fractions, because a correlation wants abundance free of library size, and
    fall back to counts only when no fractions are present."""
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")

        def candidates(suffixes):
            return [
                i
                for i, h in enumerate(header)
                if h.strip().lower() not in ANNOTATION_COLUMNS
                and any(h.strip().endswith(s) for s in suffixes)
            ]

        cols = candidates(("_frac",))
        suffix = "_frac"
        if not cols:
            cols = candidates(("_num", "_all"))
            suffix = None
        if not cols:
            sys.exit(f"[host_microbe] {path} has no per-sample abundance columns "
                     f"(looked for *_frac, *_num, *_all in: {', '.join(header[:8])}...)")

        samples = []
        for i in cols:
            name = header[i].strip()
            for s in ("_frac", "_num", "_all"):
                if name.endswith(s):
                    name = name[: -len(s)]
                    break
            samples.append(name)

        rows = {}
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(cols):
                continue
            values = []
            for i in cols:
                try:
                    values.append(float(fields[i]))
                except ValueError:
                    values.append(0.0)
            rows[fields[0]] = values
    return samples, rows


def cpm_log(counts):
    """Counts to log2(CPM+1), per sample."""
    totals = [0.0] * len(next(iter(counts.values()), []))
    for values in counts.values():
        for i, v in enumerate(values):
            totals[i] += v
    out = {}
    for gene, values in counts.items():
        out[gene] = [
            math.log2((v * 1e6 / totals[i] if totals[i] else 0.0) + 1.0)
            for i, v in enumerate(values)
        ]
    return out


def variable(values, min_nonzero):
    return sum(1 for v in values if v > 0) >= min_nonzero and len(set(values)) > 1


def spearman_floor(n):
    """Smallest attainable two-sided Spearman p-value with n observations."""
    return 2.0 / math.factorial(n) if n > 1 else 1.0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--host", required=True, help="gene x sample count matrix")
    parser.add_argument("--microbial", required=True, help="combined taxon x sample abundance table")
    parser.add_argument("--output", required=True)
    parser.add_argument("--mqc")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--min-nonzero", type=int, default=3,
                        help="a gene or taxon must be non-zero in this many samples")
    parser.add_argument("--top-taxa", type=int, default=50,
                        help="keep only the N most abundant taxa (0 = all)")
    parser.add_argument("--top-genes", type=int, default=2000,
                        help="keep only the N most variable genes (0 = all)")
    parser.add_argument("--force", action="store_true",
                        help="run even when no pair could reach significance")
    args = parser.parse_args()

    host_samples, host_rows = read_host(args.host)
    mic_samples, mic_rows = read_microbial(args.microbial)

    shared = [s for s in host_samples if s in set(mic_samples)]
    if len(shared) < 3:
        sys.exit(f"[host_microbe] only {len(shared)} sample(s) appear in both tables "
                 f"({', '.join(shared) or 'none'}). Host ids: {', '.join(host_samples[:4])}...; "
                 f"microbial ids: {', '.join(mic_samples[:4])}... A correlation needs the "
                 "same libraries on both sides.")

    hi = [host_samples.index(s) for s in shared]
    mi = [mic_samples.index(s) for s in shared]
    n = len(shared)

    host = {g: [v[i] for i in hi] for g, v in host_rows.items()}
    host = cpm_log(host)
    mic = {t: [v[i] for i in mi] for t, v in mic_rows.items()}

    host = {g: v for g, v in host.items() if variable(v, args.min_nonzero)}
    mic = {t: v for t, v in mic.items() if variable(v, args.min_nonzero)}

    if args.top_taxa and len(mic) > args.top_taxa:
        mic = dict(sorted(mic.items(), key=lambda kv: -sum(kv[1]))[: args.top_taxa])
    if args.top_genes and len(host) > args.top_genes:
        def spread(values):
            mean = sum(values) / len(values)
            return sum((v - mean) ** 2 for v in values)
        host = dict(sorted(host.items(), key=lambda kv: -spread(kv[1]))[: args.top_genes])

    pairs = len(host) * len(mic)
    if pairs == 0:
        sys.exit("[host_microbe] nothing left to test after filtering: "
                 f"{len(host)} genes x {len(mic)} taxa")

    floor = spearman_floor(n)
    max_testable = int(args.p_threshold / floor) if floor > 0 else pairs
    message = (f"{n} samples -> the smallest attainable Spearman p is 2/{n}! = {floor:.3g}; "
               f"under BH at alpha={args.p_threshold} at most {max_testable} of the "
               f"{pairs:,} gene x taxon pairs could ever be called significant")
    if max_testable < 1 and not args.force:
        sys.exit(f"[host_microbe] REFUSING: {message}. An empty result here would mean "
                 "'this cohort is too small to ask the question', not 'there are no "
                 "associations' - and the two are indistinguishable in the output file. "
                 "Add samples, cut the grid with --top-genes/--top-taxa, or pass --force "
                 "if you want the table anyway.")
    print(f"[host_microbe] {message}", file=sys.stderr)

    results = []
    for taxon, tv in mic.items():
        for gene, gv in host.items():
            rho, p = spearman(gv, tv)
            if rho is None:
                continue
            results.append([gene, taxon, rho, p])

    if not results:
        sys.exit("[host_microbe] every pair was undefined (a constant vector on one side)")

    qs = benjamini_hochberg([r[3] for r in results])
    for row, q in zip(results, qs):
        row.append(q)
    results.sort(key=lambda r: (r[4], r[3]))

    n_sig = sum(1 for r in results if r[4] < args.p_threshold)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(f"# {message}\n")
        handle.write(f"# {len(host)} genes x {len(mic)} taxa over {n} samples: "
                     f"{', '.join(shared)}\n")
        handle.write("gene_id\ttaxon\trho\tp_value\tq_value\tsignificant\n")
        for gene, taxon, rho, p, q in results:
            handle.write(f"{gene}\t{taxon}\t{rho:.4f}\t{p:.6g}\t{q:.6g}\t"
                         f"{'yes' if q < args.p_threshold else 'no'}\n")

    print(f"[host_microbe] {n_sig} of {len(results):,} pairs at q<{args.p_threshold}",
          file=sys.stderr)

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write("# id: 'host_microbe_correlation'\n")
            handle.write("# section_name: 'Host-microbe correlation'\n")
            handle.write("# format: 'tsv'\n")
            handle.write("# plot_type: 'bargraph'\n")
            handle.write(f"# description: 'Spearman over {n} samples, {len(host)} genes x "
                         f"{len(mic)} taxa. {message}.'\n")
            handle.write("Sample\tSignificant\tTested\n")
            handle.write(f"gene-taxon pairs\t{n_sig}\t{len(results) - n_sig}\n")


if __name__ == "__main__":
    main()

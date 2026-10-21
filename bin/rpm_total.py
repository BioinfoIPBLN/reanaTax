#!/usr/bin/env python3
"""
rpm_total.py -- taxon counts as reads per million of the whole library.

Every abundance the pipeline reports is a share of something smaller than the
library. Bracken's `_frac` is a share of the reads Bracken kept, and the
filtered tables renormalise over the taxa that survived. In host RNA-seq, where
the host is 99% of the library or more, those denominators hide the scale
entirely: one read in a library left with a hundred classified reads is 1%.

This writes the same table with every count divided by the sample's trimmed
read total (fastp's after-filtering reads, from read_accounting.tsv) and
multiplied by a million: RPM of the library, the unit clinical metagenomics
reports (e.g. the MetaMIC pipeline; Brincus et al., Neuromuscul Disord 2026).
It is a unit, not a filter: nothing is dropped, and the matrix-specific RPM
cut-offs that clinical studies validate are calibrated for their own assays.

Both layouts the pipeline writes are read: Bracken's combined table (`<s>_num`
counts) and combine_kreports' (`<s>_all` clade counts, so an internal node's
RPM includes its descendants', as its count does).
"""
import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from filter_abundance import count_columns, read_table, to_float


def read_library_sizes(path, column):
    """{sample: total reads} from read_accounting.tsv; raw reads when trimming is absent."""
    sizes = {}
    with open(path, encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            total = to_float(row.get(column)) or to_float(row.get("raw_reads"))
            if total > 0:
                sizes[sample] = total
    return sizes


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--table", required=True, help="combined Bracken or Kraken2 table (filtered or not)")
    parser.add_argument("--accounting", required=True, help="read_accounting.tsv")
    parser.add_argument("--column", default="trimmed_reads",
                        help="read_accounting column holding the library size (default: trimmed_reads)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    _comments, _header_line, header, rows = read_table(args.table)
    sizes = read_library_sizes(args.accounting, args.column)
    cols = count_columns(header)
    if not cols:
        sys.exit(f"[rpm_total] {args.table} has no per-sample count columns (*_num or *_all).")

    # The name, taxid and rank columns are carried as they are, and each
    # sample's count column becomes <sample>_rpm; fractions and cohort totals
    # are left out. A sample missing from the accounting table gets NA rather
    # than a made-up denominator.
    annotation = [i for i, h in enumerate(header)
                  if h.strip().lower() in ("name", "taxid", "taxonomy_id", "taxonomy_lvl", "lvl_type")]
    samples = {}
    for col in cols:
        name = header[col].strip()
        for suffix in ("_num", "_all"):
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        samples[col] = name
    missing = sorted(name for name in samples.values() if name not in sizes)
    if missing:
        print(f"[rpm_total] no library size for {len(missing)} sample(s), written as NA: {', '.join(missing[:5])}",
              file=sys.stderr)

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("\t".join([header[i].strip() for i in annotation] + [f"{samples[c]}_rpm" for c in cols]) + "\n")
        for row in rows:
            values = []
            for col in cols:
                size = sizes.get(samples[col])
                count = to_float(row[col]) if col < len(row) else 0.0
                values.append(f"{count * 1e6 / size:.6g}" if size else "NA")
            handle.write("\t".join([row[i] if i < len(row) else "" for i in annotation] + values) + "\n")

    print(f"[rpm_total] {len(rows)} taxa x {len(cols)} sample(s) written as reads per million of the library "
          f"({args.column}).", file=sys.stderr)


if __name__ == "__main__":
    main()

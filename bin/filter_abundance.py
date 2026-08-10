#!/usr/bin/env python3
"""
filter_abundance.py -- drop taxa that never reach a meaningful abundance.

The single filter that every low-biomass microbiome study applies and that costs
nothing to run: keep a taxon only if its relative abundance exceeds a threshold
in at least N samples (Monteleone et al. 2026 use >0.1% in >=1 sample). Without
it, the long tail of one-read hits - reagent contamination, database artefacts,
mis-assigned host reads - dominates any distance metric, and a taxon seen once at
0.001% carries the same weight in a Bray-Curtis matrix as a real signal.

This is deliberately NOT a contamination caller. It has no negative controls and
no concentration data, so it cannot tell a reagent contaminant from a rare real
organism; it only removes taxa that are too sparse for any downstream statistic
to say anything about. Taxa it drops are written out, not silently discarded.

Handles both table layouts this pipeline produces:
  * Bracken combined  - `name  taxonomy_id  taxonomy_lvl  <s>_num  <s>_frac ...`
  * Kraken2 combined  - `#lvl_type  name  taxid  tot_all  tot_frac  <s>_all ...`
"""
import argparse
import re
import sys


def read_table(path):
    """Return (comment_lines, header_line, header_fields, rows). Leading '#'
    lines are kept verbatim: combine_kreports puts its per-sample legend there,
    and its header line is itself '#'-prefixed."""
    comments, header_line, header, rows = [], None, None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if header is None and line.startswith("#") and "\t" not in line:
                comments.append(line)
                continue
            if header is None:
                header_line = line
                header = line.lstrip("#").split("\t")
                continue
            rows.append(line.split("\t"))
    if header is None:
        raise SystemExit(f"filter_abundance: {path} has no header row")
    return comments, header_line, header, rows


# Columns that look like per-sample values by their suffix but are not:
# Bracken's rank column ends in `_lvl`, and combine_kreports' cross-sample
# totals end in `_all`/`_frac`. Counting either would let a taxon pass on a
# number that is not a sample.
ANNOTATION_COLUMNS = {
    "name",
    "taxid",
    "taxonomy_id",
    "taxonomy_lvl",
    "lvl_type",
    "tot_all",
    "tot_frac",
    "tot_lvl",
}


def value_columns(header):
    """Column indices holding a per-sample abundance, and whether they are
    already fractions. Prefer fraction columns when the table has both."""

    def candidates(suffixes):
        return [
            i
            for i, h in enumerate(header)
            if h.strip().lower() not in ANNOTATION_COLUMNS and any(h.endswith(s) for s in suffixes)
        ]

    frac = candidates(("_frac", "_lvl"))
    if frac:
        return frac, True
    num = candidates(("_num", "_all"))
    if num:
        return num, False
    # Fall back to every column after the annotation block.
    return [i for i, h in enumerate(header) if i >= 3 and h.strip().lower() not in ANNOTATION_COLUMNS], False


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def main():
    parser = argparse.ArgumentParser(description="Abundance/prevalence filter for combined taxonomic tables.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--removed", help="write the dropped rows here, for the record")
    parser.add_argument("--min-rel-abundance", type=float, default=0.001,
                        help="relative abundance a taxon must exceed (0.001 = 0.1%%)")
    parser.add_argument("--min-samples", type=int, default=1,
                        help="in at least this many samples")
    args = parser.parse_args()

    comments, header_line, header, rows = read_table(args.input)
    cols, already_frac = value_columns(header)
    if not cols:
        raise SystemExit(f"filter_abundance: no per-sample columns found in {args.input}")

    # Counts have to be turned into per-sample relative abundances first, or the
    # threshold would mean something different in every column.
    totals = {}
    if not already_frac:
        for col in cols:
            totals[col] = sum(to_float(row[col]) for row in rows if col < len(row)) or 1.0

    kept, dropped = [], []
    for row in rows:
        hits = 0
        for col in cols:
            if col >= len(row):
                continue
            value = to_float(row[col])
            rel = value if already_frac else value / totals[col]
            if rel > args.min_rel_abundance:
                hits += 1
        (kept if hits >= args.min_samples else dropped).append(row)

    with open(args.output, "w", encoding="utf-8") as handle:
        for line in comments:
            handle.write(line + "\n")
        # Re-emit the header exactly as it came in, '#' and all: downstream
        # readers of a combine_kreports table expect it.
        handle.write(header_line + "\n")
        for row in kept:
            handle.write("\t".join(row) + "\n")

    if args.removed:
        with open(args.removed, "w", encoding="utf-8") as handle:
            handle.write("\t".join(header) + "\n")
            for row in dropped:
                handle.write("\t".join(row) + "\n")

    unit = "fraction" if already_frac else "count"
    print(
        f"[filter_abundance] {args.input}: kept {len(kept)}/{len(rows)} taxa "
        f"(> {args.min_rel_abundance:.4g} relative abundance in >= {args.min_samples} sample(s); "
        f"{len(cols)} {unit} column(s)); dropped {len(dropped)}.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

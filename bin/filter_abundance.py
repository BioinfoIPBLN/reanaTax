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
  * Kraken2 combined  - `#perc  tot_all  tot_lvl  <s>_all  <s>_lvl  lvl_type  taxid  name`
"""
import argparse
import os
import re
import sys


def read_table(path):
    """Return (comment_lines, header_line, header_fields, rows).

    combine_kreports writes a two-line preamble, then one `#<sample>\t<file>`
    legend line per sample, then a header that is itself '#'-prefixed. The
    legend lines carry tabs too, so the header is identified as the LAST leading
    '#' line rather than the first one holding a tab. Bracken's combined table
    has no '#' lines at all and its header is simply the first line.
    """
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    if not lines:
        raise SystemExit(f"filter_abundance: {path} is empty")

    leading = 0
    while leading < len(lines) and lines[leading].startswith("#"):
        leading += 1

    if leading:
        comments = lines[: leading - 1]
        header_line = lines[leading - 1]
        header = header_line.lstrip("#").split("\t")
    else:
        comments = []
        header_line = lines[0]
        header = header_line.split("\t")

    rows = [line.split("\t") for line in lines[max(leading, 1) :]]
    return comments, header_line, header, rows


# Columns that look like per-sample values by their suffix but are not:
# Bracken's rank column is `taxonomy_lvl`, and combine_kreports' cross-sample
# columns are `perc`/`tot_*`. Counting any of them would let a taxon pass on a
# number that is not a sample.
ANNOTATION_COLUMNS = {
    "name",
    "taxid",
    "taxonomy_id",
    "taxonomy_lvl",
    "lvl_type",
    "perc",
    "tot_all",
    "tot_frac",
    "tot_lvl",
}


def value_columns(header):
    """Column indices holding a per-sample abundance, and whether they are
    already fractions.

    Bracken publishes `<s>_num` and `<s>_frac` side by side, so the fractions
    are taken directly. combine_kreports publishes `<s>_all` (reads anywhere in
    the clade) and `<s>_lvl` (reads at exactly this rank); both are counts, and
    `_all` is the one that means abundance - `_lvl` is zero for every internal
    node, so filtering on it would drop every taxon above species.
    """

    def candidates(suffixes):
        return [
            i
            for i, h in enumerate(header)
            if h.strip().lower() not in ANNOTATION_COLUMNS and any(h.strip().endswith(s) for s in suffixes)
        ]

    frac = candidates(("_frac",))
    if frac:
        return frac, True
    return candidates(("_num", "_all")), False


def count_columns(header):
    """Per-sample COUNT columns, whichever layout the table is in.

    Separate from value_columns() on purpose. Bracken publishes `_num` and
    `_frac` side by side and value_columns() prefers the fractions, so the
    read floor cannot piggyback on it - the floor needs reads, and asking a
    fraction how many reads it represents is meaningless.
    """
    return [
        i
        for i, h in enumerate(header)
        if h.strip().lower() not in ANNOTATION_COLUMNS
        and any(h.strip().endswith(suffix) for suffix in ("_num", "_all"))
    ]


def column_totals(header, rows, cols):
    """Per-sample denominators for count columns.

    A Bracken table is flat - one rank, every read counted once - so the column
    sums to the sample total. A combine_kreports table is a hierarchy, where
    every read is counted again at each of its ancestors, so the column sum is
    several times the sample total. There the root row's clade count is used
    instead, which is the classified-read total and so matches the denominator
    behind Bracken's `_frac`.
    """
    lvl_type = next((i for i, h in enumerate(header) if h.strip().lower() == "lvl_type"), None)
    if lvl_type is None:
        return {col: sum(to_float(row[col]) for row in rows if col < len(row)) or 1.0 for col in cols}

    root = next(
        (row for row in rows if lvl_type < len(row) and row[lvl_type].strip().upper() == "R"),
        None,
    )
    totals = {}
    for col in cols:
        if root is not None and col < len(root):
            totals[col] = to_float(root[col]) or 1.0
        else:
            totals[col] = max((to_float(row[col]) for row in rows if col < len(row)), default=0.0) or 1.0
    return totals


def taxid_column(header):
    """Bracken calls it `taxonomy_id`, combine_kreports calls it `taxid`."""
    for name in ("taxonomy_id", "taxid"):
        for i, h in enumerate(header):
            if h.strip().lower() == name:
                return i
    return None


def drop_taxa(header, rows, taxids):
    """Split out the rows belonging to one taxon.

    Host carry-over is the reason this exists: reads the depletion step missed
    are still classified, and while they sit in the table every microbial
    abundance is a fraction of host + microbes rather than of microbes. Removing
    the row and letting the fractions be recomputed against the remaining total
    is exact for a Bracken table, which is flat and holds each read once.

    A combine_kreports table is a hierarchy, so only the taxon's own rows go;
    its ancestors (root, Eukaryota, ...) keep their original clade counts and
    would still include it. That table is a report, not the analysis input -
    docs/output.md says to take Bracken downstream - so the discrepancy is
    documented rather than papered over by rewriting ancestor rows.
    """
    column = taxid_column(header)
    if column is None:
        print("[filter_abundance] no taxid column; --drop-taxid/--drop-taxids ignored.", file=sys.stderr)
        return rows, []

    wanted = {str(taxid).strip() for taxid in taxids}
    kept, dropped = [], []
    for row in rows:
        target = dropped if column < len(row) and row[column].strip() in wanted else kept
        target.append(row)
    return kept, dropped


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
    parser.add_argument("--min-reads", type=float, default=0,
                        help="AND at least this many reads summed across samples. A relative "
                             "threshold alone is cleared by a handful of reads once host "
                             "depletion has removed 99%% of the library, which is exactly the "
                             "regime where 0.1%% of the microbial fraction is three reads on one "
                             "conserved locus. 0 disables.")
    parser.add_argument("--drop-taxid", type=int,
                        help="remove this taxon (typically the host) and renormalise the "
                             "fraction columns, so abundances are fractions of what is left")
    parser.add_argument("--drop-taxids", nargs="+", default=[],
                        help="one or more files of taxids to remove, one per line - the lists "
                             "bin/minimizer_filter.py and bin/host_kmer_filter.py write. Several "
                             "are taken as a union: each condemns a taxon for its own reason, and "
                             "surviving one test is no argument against the other.")
    args = parser.parse_args()

    comments, header_line, header, rows = read_table(args.input)
    dropped_host = []
    # The host taxon and anything the minimizer evidence condemned are removed
    # together: both are "this is not part of the microbial profile", and both
    # want the fraction columns renormalised against what survives.
    drop = set()
    if args.drop_taxid is not None:
        drop.add(str(args.drop_taxid))
    for path in args.drop_taxids:
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8") as handle:
            drop.update(line.strip() for line in handle if line.strip() and not line.startswith("#"))
    if drop:
        rows, dropped_host = drop_taxa(header, rows, drop)
    cols, already_frac = value_columns(header)
    if not cols:
        raise SystemExit(f"filter_abundance: no per-sample columns found in {args.input}")

    # Counts have to be turned into per-sample relative abundances first, or the
    # threshold would mean something different in every column.
    totals = {} if already_frac else column_totals(header, rows, cols)

    # With the host gone the fraction columns no longer sum to 1, so rescale
    # them against the surviving total. Count columns need no rescaling: they
    # are turned into fractions above, against a total that already excludes
    # the dropped rows.
    if dropped_host and already_frac:
        for col in cols:
            total = sum(to_float(row[col]) for row in rows if col < len(row))
            if total <= 0:
                continue
            for row in rows:
                if col < len(row):
                    row[col] = f"{to_float(row[col]) / total:.10g}"

    # An absolute floor alongside the relative one, applied as AND. Every
    # careful reanalysis in this area imposes one: Gihawi et al. (mBio 2023)
    # filter to genus counts >= 10 "on the assumption that smaller values
    # likely represent noise or contamination", and PRISM will not even attempt
    # confirmation below 10 reads. The reason is that the relative threshold is
    # computed against the MICROBIAL total, and after host depletion has removed
    # 99% of the library that total is small - Salzberg's TCGA reanalysis left a
    # median of 2.6M non-host reads, 0.48% of the library - so 0.1% of it can be
    # three reads on one conserved locus.
    #
    # Only meaningful on a count table; a MetaPhlAn-style fraction table has no
    # reads to floor, and the threshold is skipped rather than misapplied.
    count_cols = count_columns(header)
    read_floor = args.min_reads if (args.min_reads > 0 and count_cols) else 0
    if args.min_reads > 0 and not count_cols:
        print("[filter_abundance] --min-reads ignored: no count columns in this table "
              "(a MetaPhlAn-style profile carries fractions only).", file=sys.stderr)

    kept, dropped, below_floor = [], [], 0
    for row in rows:
        hits = 0
        for col in cols:
            if col >= len(row):
                continue
            value = to_float(row[col])
            rel = value if already_frac else value / totals[col]
            if rel > args.min_rel_abundance:
                hits += 1
        total_reads = sum(to_float(row[col]) for col in count_cols if col < len(row))
        passes = hits >= args.min_samples
        if passes and read_floor and total_reads < read_floor:
            passes = False
            below_floor += 1
        (kept if passes else dropped).append(row)

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
            for row in dropped + dropped_host:
                handle.write("\t".join(row) + "\n")

    if dropped_host:
        print(
            f"[filter_abundance] removed {len(drop)} taxon/taxa "
            f"({len(dropped_host)} row(s)) and renormalised.",
            file=sys.stderr,
        )
    unit = "fraction" if already_frac else "count"
    print(
        f"[filter_abundance] {args.input}: kept {len(kept)}/{len(rows)} taxa "
        f"(> {args.min_rel_abundance:.4g} relative abundance in >= {args.min_samples} sample(s)"
        + (f", and >= {read_floor:g} reads total" if read_floor else "")
        + f"; {len(cols)} {unit} column(s)); dropped {len(dropped)}"
        + (f", of which {below_floor} for the read floor alone" if read_floor else "") + ".",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
shuffle_compare.py -- what the classifier found in sequence that was not there.

bin/shuffle_reads.py made a copy of every library with the sequence destroyed
and the composition intact; the pipeline classified that copy against the same
database with the same settings. This compares the two classifications.

The comparison is per taxon, across the cohort, and the statistic is not a
p-value. It is a RATIO: reads the taxon collected from shuffled reads, divided
by reads it collected from real ones. A ratio near zero means the taxon needs
genuine sequence to be found, which is what a real organism looks like. A ratio
near one means the database would have reported it from noise of the same base
composition, and its real count is not evidence.

Three things this deliberately does not do.

  It does not test significance.  The shuffled library is one draw, and a
  per-taxon Poisson test against it would dress a single observation up as an
  inference. The ratio is a measurement, reported as such.

  It does not compare per sample.  Chance matches are rare per taxon per
  sample, so a per-sample ratio is mostly 0/0. Counts are pooled over the
  cohort, and the per-sample columns are written for inspection only.

  It does not know about the host.  A taxon whose reads are host carry-over
  survives this control perfectly - carry-over reads are real sequence and
  shuffling them removes them. That is the host-k-mer filter's job, and the
  two are not substitutes.

The detection floor is computed and enforced. If the shuffled library was
subsampled, one chance read in it stands for `scale` reads in the real one, so
no taxon with fewer than scale/max_ratio real reads can ever be flagged. That
number is printed, and taxa below it are marked `untested` rather than `clean`
- an untested taxon has not passed anything.
"""
import argparse
import os
import sys

# Kraken2 rows for the two pseudo-taxa that are not organisms.
UNCLASSIFIED_TAXID = "0"
ROOT_TAXID = "1"


def parse_report(path):
    """(clade_reads_by_taxid, name_by_taxid, rank_by_taxid, total_fragments).

    Handles both report layouts this pipeline can produce: the plain six-column
    Kraken2 report and the eight-column one from `--report-minimizer-data`. The
    trailing three fields are always rank, taxid, name, so they are read from
    the end and the extra columns in the middle are simply not looked at.
    """
    clade, names, ranks = {}, {}, {}
    total = 0
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("%"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            taxid = fields[-2].strip()
            names[taxid] = fields[-1].strip()
            ranks[taxid] = fields[-3].strip()
            try:
                reads = int(fields[1])
            except ValueError:
                continue
            clade[taxid] = clade.get(taxid, 0) + reads
            # Every fragment is either unclassified or under the root, and the
            # two rows never overlap - so this is the library size, taken from
            # the report itself rather than assumed from elsewhere.
            if taxid in (UNCLASSIFIED_TAXID, ROOT_TAXID):
                total += reads
    return clade, names, ranks, total


def sample_name(path):
    base = os.path.basename(path)
    for suffix in (".kraken2.report.txt", ".kreport.txt", ".report.txt", ".txt", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--real", nargs="+", required=True, help="Kraken2 reports of the real reads")
    parser.add_argument("--shuffled", nargs="+", required=True, help="...and of the shuffled copies")
    parser.add_argument("--max-ratio", type=float, default=0.1,
                        help="flag a taxon once shuffled/real reaches this")
    parser.add_argument("--min-reads", type=int, default=10,
                        help="below this many real reads a taxon is not judged at all")
    parser.add_argument("--prefix", default="shuffle")
    args = parser.parse_args()

    real_by_sample = {sample_name(path): path for path in args.real}
    shuffled_by_sample = {sample_name(path): path for path in args.shuffled}
    shared = sorted(set(real_by_sample) & set(shuffled_by_sample))
    if not shared:
        sys.exit(
            "[shuffle_compare] no sample has both a real and a shuffled report. The two sets are "
            f"real={sorted(real_by_sample)} shuffled={sorted(shuffled_by_sample)}. They are paired "
            "by file basename, so the shuffled reports must be named after the same sample ids."
        )
    for orphan in sorted(set(real_by_sample) ^ set(shuffled_by_sample)):
        print(f"[shuffle_compare] '{orphan}' has only one of the two classifications; skipped.",
              file=sys.stderr)

    real_total, shuffled_total = {}, {}
    real_fragments, shuffled_fragments = 0, 0
    names, ranks = {}, {}
    per_sample = {}

    for sample in shared:
        real_clade, real_names, real_ranks, real_n = parse_report(real_by_sample[sample])
        shuffled_clade, shuffled_names, shuffled_ranks, shuffled_n = parse_report(shuffled_by_sample[sample])
        names.update(real_names)
        names.update(shuffled_names)
        ranks.update(real_ranks)
        ranks.update(shuffled_ranks)
        real_fragments += real_n
        shuffled_fragments += shuffled_n
        for taxid, reads in real_clade.items():
            real_total[taxid] = real_total.get(taxid, 0) + reads
        for taxid, reads in shuffled_clade.items():
            shuffled_total[taxid] = shuffled_total.get(taxid, 0) + reads
        per_sample[sample] = (real_clade, shuffled_clade, real_n, shuffled_n)

    if shuffled_fragments == 0:
        sys.exit("[shuffle_compare] the shuffled libraries hold no fragments at all.")

    # One shuffled read stands for this many real ones. Exactly 1.0 when the
    # whole library was shuffled, which is the default.
    scale = real_fragments / shuffled_fragments

    # The smallest real count at which a taxon could possibly be flagged: it
    # takes one chance read, worth `scale` after scaling, to reach max_ratio.
    floor = scale / args.max_ratio if args.max_ratio > 0 else float("inf")

    print(
        f"[shuffle_compare] {len(shared)} sample(s); {real_fragments} real fragment(s) against "
        f"{shuffled_fragments} shuffled (scale {scale:.4f}). A taxon needs at least "
        f"{floor:.0f} real reads before a single chance match could flag it; below that it is "
        "reported as untested.",
        file=sys.stderr,
    )

    rows = []
    for taxid in sorted(set(real_total) | set(shuffled_total), key=lambda t: -real_total.get(t, 0)):
        if taxid in (UNCLASSIFIED_TAXID, ROOT_TAXID):
            continue
        real = real_total.get(taxid, 0)
        shuffled = shuffled_total.get(taxid, 0)
        scaled = shuffled * scale
        ratio = scaled / real if real else float("inf") if shuffled else 0.0
        if real < args.min_reads:
            verdict = "untested_low_depth"
        elif real < floor:
            verdict = "untested_below_floor"
        elif ratio >= args.max_ratio:
            verdict = "composition_only"
        else:
            verdict = "clean"
        rows.append({
            "taxid": taxid,
            "name": names.get(taxid, ""),
            "rank": ranks.get(taxid, ""),
            "real_reads": real,
            "shuffled_reads": shuffled,
            "shuffled_scaled": round(scaled, 2),
            "ratio": round(ratio, 4) if ratio != float("inf") else "inf",
            "samples_real": sum(1 for s in shared if per_sample[s][0].get(taxid, 0) > 0),
            "samples_shuffled": sum(1 for s in shared if per_sample[s][1].get(taxid, 0) > 0),
            "verdict": verdict,
        })

    columns = ["taxid", "name", "rank", "real_reads", "shuffled_reads", "shuffled_scaled",
               "ratio", "samples_real", "samples_shuffled", "verdict"]
    with open(f"{args.prefix}.shuffle_evidence.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    dropped = [row for row in rows if row["verdict"] == "composition_only"]
    with open(f"{args.prefix}.shuffle_drop.txt", "w", encoding="utf-8") as handle:
        for row in dropped:
            handle.write(f"{row['taxid']}\n")

    untested = sum(1 for row in rows if row["verdict"].startswith("untested"))
    clean = len(rows) - len(dropped) - untested
    with open(f"{args.prefix}_shuffle_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_shuffle'",
            "# section_name: 'Shuffled-read negative control'",
            "# description: 'The same libraries were classified twice: once as sequenced, and once",
            "#     with every read shuffled so that its length, GC content and dinucleotide",
            "#     frequencies survive and its k-mers do not. A taxon that still collects reads",
            "#     from the shuffled copy is being called on base composition rather than on",
            "#     homology. Taxa are counted as composition-only once their shuffled count reaches",
            f"#     {args.max_ratio:g} of their real one. Untested taxa are those with too few real",
            "#     reads for a single chance match to reach that ratio - they have not passed",
            "#     anything, which is why they are counted apart.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_shuffle_plot'",
            "#     title: 'reanaTax: shuffled-read control'",
            "#     ylab: 'Taxa'",
            "Sample\tComposition only\tClean\tUntested",
            f"all taxa\t{len(dropped)}\t{clean}\t{untested}",
            "",
        ]))

    print(
        f"[shuffle_compare] {len(dropped)}/{len(rows)} taxa reported by the classifier from "
        f"shuffled sequence at ratio >= {args.max_ratio:g}; {untested} untested; "
        f"{len(dropped)} taxid(s) listed for removal.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

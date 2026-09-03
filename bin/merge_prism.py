#!/usr/bin/env python3
"""
merge_prism.py -- PRISM's per-sample verdicts into one table per cohort.

PRISM writes `<sample>-counts.csv` with four columns: `tax_name`, `staxids`,
`n` (reads confirmed for that taxon) and `pred` (the XGBoost model's
probability that the taxon is genuinely present rather than a contaminant).

Two matrices come out, taxa by sample: the read counts, and the scores. They
are kept apart because they answer different questions and because a score has
no meaning where there were no reads to score - a taxon absent from a sample
gets an empty cell in the score matrix, not a zero. A zero would read as
"PRISM was confident this is a contaminant", which is the opposite of "PRISM
never saw it".

The drop list uses the score, and the threshold is the caller's. PRISM's
authors do not publish one cut-off for every context; what they do say, and
what this pipeline already encodes in `--min_reads`, is that below about ten
reads there is nothing to confirm either way. So a taxon is only condemned
where it was actually scored, and taxa carried by too few reads are reported as
untested rather than as either verdict.

Aggregating across samples: a taxon is dropped when its MEDIAN score across the
samples that saw it falls below the threshold. Not the mean, which one
confidently-scored library can drag either way, and not "any sample", which on
a large cohort condemns everything.
"""
import argparse
import csv
import os
import statistics
import sys


def sample_name(path):
    base = os.path.basename(path)
    for suffix in ("-counts.csv", ".counts.csv", ".csv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def read_counts(path):
    """{taxid: (name, reads, score_or_None)}."""
    out = {}
    with open(path, encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return out
        lookup = {name.strip().lower(): name for name in reader.fieldnames}
        for key in ("staxids", "n"):
            if key not in lookup:
                sys.exit(
                    f"[merge_prism] '{path}' has no `{key}` column; its header is "
                    f"{reader.fieldnames}. Expected PRISM's <sample>-counts.csv."
                )
        for row in reader:
            taxid = (row.get(lookup["staxids"]) or "").strip()
            if not taxid:
                continue
            try:
                reads = float(row[lookup["n"]])
            except (KeyError, ValueError):
                continue
            score = None
            if "pred" in lookup:
                raw = (row.get(lookup["pred"]) or "").strip()
                try:
                    score = float(raw)
                except ValueError:
                    score = None
            name = (row.get(lookup.get("tax_name", ""), "") or "").strip()
            out[taxid] = (name, reads, score)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", nargs="+", required=True, help="PRISM <sample>-counts.csv files")
    parser.add_argument("--score-threshold", type=float, default=0.5,
                        help="median score below which a taxon is called a contaminant")
    parser.add_argument("--min-reads", type=float, default=10,
                        help="below this many reads in every sample, leave the taxon untested")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    samples, per_sample, names = [], {}, {}
    for path in sorted(args.counts):
        sample = sample_name(path)
        if sample in per_sample:
            sys.exit(f"[merge_prism] two files resolve to the sample id '{sample}'.")
        table = read_counts(path)
        samples.append(sample)
        per_sample[sample] = table
        for taxid, (name, _reads, _score) in table.items():
            names.setdefault(taxid, name)

    if not names:
        sys.exit(
            "[merge_prism] PRISM confirmed no taxon in any sample. That is a result, not an "
            "error, but it is also what an unprepared reference directory produces - check the "
            "per-sample logs before reading it as biology."
        )

    taxids = sorted(names, key=lambda t: -sum(
        per_sample[s].get(t, ("", 0.0, None))[1] for s in samples
    ))

    for label, index in (("reads", 1), ("score", 2)):
        with open(f"{args.prefix}.prism_{label}.tsv", "w", encoding="utf-8") as handle:
            handle.write("\t".join(["taxonomy_id", "name"] + samples) + "\n")
            for taxid in taxids:
                row = [taxid, names[taxid]]
                for sample in samples:
                    entry = per_sample[sample].get(taxid)
                    if entry is None:
                        # Absent, not zero. See the module docstring.
                        row.append("0" if label == "reads" else "NA")
                    else:
                        value = entry[index]
                        row.append("NA" if value is None else f"{value:.6g}")
                handle.write("\t".join(row) + "\n")

    rows, dropped = [], []
    for taxid in taxids:
        seen = [per_sample[s][taxid] for s in samples if taxid in per_sample[s]]
        reads = [entry[1] for entry in seen]
        scores = [entry[2] for entry in seen if entry[2] is not None]
        median = statistics.median(scores) if scores else None
        if not reads or max(reads) < args.min_reads:
            verdict = "untested_low_depth"
        elif median is None:
            verdict = "unscored"
        elif median < args.score_threshold:
            verdict = "contaminant"
            dropped.append(taxid)
        else:
            verdict = "present"
        rows.append({
            "taxid": taxid,
            "name": names[taxid],
            "samples_seen": len(seen),
            "total_reads": round(sum(reads), 2),
            "max_reads": round(max(reads), 2) if reads else 0,
            "median_score": "NA" if median is None else round(median, 4),
            "min_score": "NA" if not scores else round(min(scores), 4),
            "max_score": "NA" if not scores else round(max(scores), 4),
            "verdict": verdict,
        })

    columns = ["taxid", "name", "samples_seen", "total_reads", "max_reads",
               "median_score", "min_score", "max_score", "verdict"]
    with open(f"{args.prefix}.prism_evidence.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    with open(f"{args.prefix}.prism_drop.txt", "w", encoding="utf-8") as handle:
        for taxid in dropped:
            handle.write(f"{taxid}\n")

    present = sum(1 for row in rows if row["verdict"] == "present")
    untested = sum(1 for row in rows if row["verdict"] in ("untested_low_depth", "unscored"))
    with open(f"{args.prefix}_prism_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_prism'",
            "# section_name: 'PRISM confirmation'",
            "# description: 'PRISM scores each candidate taxon with an XGBoost model over forty",
            "#     features - BLAST alignment against nt, the GenBank annotation those alignments",
            "#     land on, host mapping, and k-mer composition at every rank - and returns the",
            "#     probability that it is genuinely present rather than a contaminant. Taxa are",
            f"#     called contaminants below a median score of {args.score_threshold:g} across the",
            "#     samples that saw them; taxa with too few reads anywhere to confirm are counted",
            "#     apart, because they have not passed anything.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_prism_plot'",
            "#     title: 'reanaTax: PRISM'",
            "#     ylab: 'Taxa'",
            "Sample\tPresent\tContaminant\tUntested",
            f"all taxa\t{present}\t{len(dropped)}\t{untested}",
            "",
        ]))

    print(
        f"[merge_prism] {len(rows)} taxa across {len(samples)} sample(s): {present} present, "
        f"{len(dropped)} called contaminants at median score < {args.score_threshold:g}, "
        f"{untested} untested.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

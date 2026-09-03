#!/usr/bin/env python3
"""
read_accounting.py -- one row per sample tracking what happened to every read.

raw -> after trimming -> aligned to each host reference -> non-host -> classified

Each of those numbers already exists somewhere (fastp's JSON, HISAT2's summary
log, the Kraken2 report), but never side by side, which is exactly when you need
them: to decide whether a sample has enough non-host reads for its taxonomic
profile to mean anything, and whether the host depletion actually worked.

It also computes the HOST CARRY-OVER: the share of supposedly non-host reads
that Kraken2 still assigns to the host taxon. Kraken2's standard databases
contain the host genome on purpose, so this is a free check on the aligner --
a few percent is normal, tens of percent means host depletion is failing (wrong
reference, or a genome whose repeats and structural variants are missing from
the assembly you used).

Outputs a plain TSV plus two MultiQC custom-content files: a stacked bargraph of
read fate and a general-statistics column for the carry-over.
"""
import argparse
import json
import os
import re
import sys

# Kraken2 rank codes for a "real" taxon row; U is unclassified, R is the root.
HOST_DEFAULT_TAXID = 9606  # Homo sapiens


def sample_name(path, suffixes):
    """Strip the tool-specific suffix off a file name to recover the sample id."""
    name = os.path.basename(path)
    for suffix in sorted(suffixes, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return re.sub(r"\.(tsv|txt|log|json)$", "", name)


def read_fastp(path):
    """Reads before and after trimming. fastp counts READS, both mates included."""
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, ValueError) as exc:
        print(f"[read_accounting] cannot parse {path}: {exc}", file=sys.stderr)
        return {}
    after = data.get("summary", {}).get("after_filtering", {}) or {}
    return {
        "raw_reads": (data.get("summary", {}).get("before_filtering", {}) or {}).get("total_reads", 0),
        "trimmed_reads": after.get("total_reads", 0),
        "trimmed_mean_length": after.get("read1_mean_length", 0),
        "per_fragment": 2 if after.get("read2_mean_length") else 1,
    }


def read_hisat2(path):
    """HISAT2 writes a Bowtie2-style summary. It counts PAIRS for paired-end
    input and reads for single-end, so normalise to reads to stay comparable
    with fastp and Kraken2."""
    try:
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
    except OSError as exc:
        print(f"[read_accounting] cannot read {path}: {exc}", file=sys.stderr)
        return {}

    total = re.search(r"^\s*(\d+)\s+reads;\s+of these", text, re.M)
    paired = re.search(r"^\s*(\d+)\s+\([\d.]+%\)\s+were paired", text, re.M)
    unaligned_conc = re.search(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly 0 times", text, re.M)
    unaligned_se = re.search(r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned 0 times", text, re.M)
    rate = re.search(r"([\d.]+)%\s+overall alignment rate", text)

    if not total:
        return {}
    fragments = int(total.group(1))
    is_paired = bool(paired) and int(paired.group(1)) > 0
    per_fragment = 2 if is_paired else 1

    if is_paired:
        unaligned = int(unaligned_conc.group(1)) if unaligned_conc else 0
    else:
        unaligned = int(unaligned_se.group(1)) if unaligned_se else 0

    return {
        "input_reads": fragments * per_fragment,
        "unaligned_reads": unaligned * per_fragment,
        "aligned_reads": (fragments - unaligned) * per_fragment,
        "alignment_rate": float(rate.group(1)) if rate else 0.0,
        "per_fragment": per_fragment,
    }


def read_sortmerna(path):
    """SortMeRNA's log. `Total reads passing E-value threshold` are the reads
    that ALIGNED to the rRNA references, i.e. the rRNA; the ones `failing` are
    what survives. (Same reading as MultiQC's sortmerna module.) These are
    per-read counts of the E-value test, so with `--paired_in` they do not equal
    the number of reads actually written out - a pair is dropped whole when
    either mate is rRNA. They are reported for information only; no percentage
    downstream is derived from them, because SortMeRNA runs BEFORE host
    depletion and HISAT2's own input count is the exact figure."""
    try:
        with open(path, encoding="utf-8") as handle:
            text = handle.read()
    except OSError as exc:
        print(f"[read_accounting] cannot read {path}: {exc}", file=sys.stderr)
        return {}

    total = re.search(r"Total reads\s*=\s*(\d+)", text)
    rrna = re.search(r"Total reads passing[^=]*=\s*(\d+)", text)
    if not total or not rrna:
        return {}
    return {
        "rrna_input_reads": int(total.group(1)),
        "rrna_reads": int(rrna.group(1)),
    }


def read_kraken2(path, host_taxid):
    """A Kraken-style report, from either Kraken2 or KrakenUniq.

    Column 2 is the count covered by the clade rooted at that taxon in BOTH
    formats. What differs is where the taxid sits: Kraken2 ends `... rank taxid
    name`, KrakenUniq ends `... taxid rank name` (and spells its ranks out as
    'species' rather than 'S'). So the taxid is found by trying the second-to-
    last field and falling back to the third-to-last - the same dual-branch
    trick Bracken and KrakenTools use - and the three numbers wanted here are
    keyed on TAXID rather than on a rank code, which is identical across the
    two: 0 is unclassified, 1 is the root, i.e. everything classified.

    That also survives `--report-minimizer-data`, which inserts two columns in
    the middle and leaves both ends alone.

    Those counts are FRAGMENTS: in paired mode a pair is classified once and
    reported once. Every other source in this table counts mates, so derive()
    scales them by per_fragment; leaving them unscaled halves both
    classified_pct and host_carryover_pct on paired data.
    """
    classified = unclassified = host = 0
    try:
        with open(path, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("#") or line.startswith("%"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                try:
                    clade_reads = int(fields[1])
                except ValueError:
                    continue  # header or preamble
                try:
                    taxid = int(fields[-2])          # Kraken2: rank, taxid, name
                except ValueError:
                    try:
                        taxid = int(fields[-3])      # KrakenUniq: taxid, rank, name
                    except ValueError:
                        continue

                if taxid == 0:
                    unclassified = max(unclassified, clade_reads)
                elif taxid == 1:
                    classified = max(classified, clade_reads)
                if taxid == host_taxid:
                    host = max(host, clade_reads)
    except OSError as exc:
        print(f"[read_accounting] cannot read {path}: {exc}", file=sys.stderr)
        return {}
    return {
        "classified_reads": classified,
        "unclassified_reads": unclassified,
        "host_carryover_reads": host,
    }


def collect(args):
    rows = {}

    def slot(sample):
        return rows.setdefault(sample, {"sample": sample})

    for path in args.fastp:
        slot(sample_name(path, [".fastp.json"])).update(read_fastp(path))

    # HISAT2 logs are named <sample>.<label>.hisat2.summary.log, where the label
    # says which host reference the pass used. Keep the passes apart so a
    # two-reference run shows what each one removed.
    for path in args.hisat2:
        name = os.path.basename(path)
        match = re.match(r"^(.+?)\.([^.]+)\.hisat2\.summary\.log$", name)
        if match:
            sample, label = match.group(1), match.group(2)
        else:
            sample, label = sample_name(path, [".hisat2.summary.log"]), "host"
        stats = read_hisat2(path)
        entry = slot(sample)
        for key, value in stats.items():
            entry[f"{label}_{key}"] = value

    for path in args.sortmerna:
        slot(sample_name(path, [".sortmerna.log"])).update(read_sortmerna(path))

    for path in args.kraken2:
        slot(sample_name(path, [".kraken2.report.txt", ".kraken2.report", ".report.txt"])).update(
            read_kraken2(path, args.host_taxid)
        )

    return [rows[key] for key in sorted(rows)]


def pass_order(label):
    """Depletion passes in the order they ran. Intermediate passes are labelled
    `host1`, `host2`, ...; the last host pass is always plain `host`; and
    `univec`, when --univec is on, runs after every host pass. Sorting
    alphabetically would put `host` before `host1` and silently report the
    wrong pass as final, and would put `host` after `univec` only by luck.
    `nonhost_reads` is read from the LAST entry of this ordering, so getting it
    wrong misreports the size of the fraction everything downstream sees."""
    match = re.fullmatch(r"host(\d+)", label)
    if match:
        return (0, int(match.group(1)))
    return (2, 0) if label == "univec" else (1, 0)


def derive(row):
    """Fill in the totals that only make sense once every source is joined."""
    host_labels = sorted({k[: -len("_aligned_reads")] for k in row if k.endswith("_aligned_reads")}, key=pass_order)
    # Deliberately not `host_aligned_reads`: that is the FINAL pass's own key,
    # and overwriting it would throw away a per-pass number.
    row["host_removed_reads"] = sum(row.get(f"{label}_aligned_reads", 0) for label in host_labels)
    # The reads that survived every host pass: the last pass's unaligned count,
    # or the trimmed count when no host depletion ran at all.
    if host_labels:
        row["nonhost_reads"] = row.get(f"{host_labels[-1]}_unaligned_reads", 0)
    else:
        row["nonhost_reads"] = row.get("trimmed_reads", 0)

    # Kraken2 counted fragments; everything above counts mates. fastp knows the
    # library layout, and the HISAT2 summary is the fallback when trimming was
    # skipped (its key is prefixed with the pass label).
    per_fragment = row.get("per_fragment") or next(
        (row[key] for key in row if key.endswith("_per_fragment")), 1
    )
    for key in ("classified_reads", "unclassified_reads", "host_carryover_reads"):
        if key in row:
            row[key] *= per_fragment
    # Keep one column saying what the layout was; drop the per-pass duplicates.
    for key in [k for k in list(row) if k.endswith("_per_fragment")]:
        del row[key]
    row["per_fragment"] = per_fragment

    # How many reads rRNA depletion actually removed. Taken as trimmed minus
    # what the FIRST host pass received, not from the SortMeRNA log: with
    # `--paired_in` a pair is dropped whole when either mate is rRNA, so the
    # log's per-read E-value counts overstate what left the step. The
    # subtraction is exact whenever both numbers exist, and the log is only the
    # fallback for a run with no host depletion at all.
    if "rrna_input_reads" in row:
        trimmed = row.get("trimmed_reads", 0)
        first_pass = row.get(f"{host_labels[0]}_input_reads") if host_labels else None
        if trimmed and first_pass is not None:
            row["rrna_removed_reads"] = max(trimmed - first_pass, 0)
        else:
            row["rrna_removed_reads"] = row.get("rrna_reads", 0)
        denominator = trimmed or row.get("rrna_input_reads", 0)
        row["rrna_pct"] = round(100.0 * row["rrna_removed_reads"] / denominator, 2) if denominator else 0.0
    for key in ("rrna_input_reads", "rrna_reads"):
        row.pop(key, None)

    nonhost = row["nonhost_reads"] or 0
    classified = row.get("classified_reads", 0)
    host_carry = row.get("host_carryover_reads", 0)
    row["classified_pct"] = round(100.0 * classified / nonhost, 2) if nonhost else 0.0
    row["host_carryover_pct"] = round(100.0 * host_carry / nonhost, 2) if nonhost else 0.0
    raw = row.get("raw_reads", 0)
    row["nonhost_pct_of_raw"] = round(100.0 * nonhost / raw, 2) if raw else 0.0
    return row


COLUMNS = [
    "sample",
    "raw_reads",
    "trimmed_reads",
    "trimmed_mean_length",
    "rrna_removed_reads",
    "rrna_pct",
    "host_removed_reads",
    "nonhost_reads",
    "nonhost_pct_of_raw",
    "classified_reads",
    "unclassified_reads",
    "classified_pct",
    "host_carryover_reads",
    "host_carryover_pct",
]


def write_tsv(path, rows):
    extra = sorted({k for row in rows for k in row} - set(COLUMNS))
    header = COLUMNS + extra
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, 0)) for col in header) + "\n")


def write_mqc_bargraph(path, rows):
    """Stacked bargraph of where every raw read ended up. The categories are
    disjoint and sum to the raw count, so the bar length is the library size."""
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "# id: 'reanatax_read_fate'\n"
            "# section_name: 'Read fate'\n"
            "# description: 'What happened to every raw read: removed by fastp, removed as rRNA,\n"
            "#     aligned to a host reference, or carried through to classification. Categories are\n"
            "#     disjoint, so the bar length is the raw library size.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'reanatax_read_fate_plot'\n"
            "#     title: 'reanaTax: read fate'\n"
            "#     ylab: 'Reads'\n"
            "# section_href: 'https://github.com/BioinfoIPBLN/reanatax'\n"
        )
        handle.write(
            "Sample\tFiltered by fastp\tRemoved as rRNA\tAligned to host"
            "\tNon-host, classified\tNon-host, unclassified\n"
        )
        for row in rows:
            raw = row.get("raw_reads", 0)
            trimmed = row.get("trimmed_reads", 0)
            filtered = max(raw - trimmed, 0)
            rrna = row.get("rrna_removed_reads", 0)
            host = row.get("host_removed_reads", 0)
            classified = row.get("classified_reads", 0)
            unclassified = row.get("unclassified_reads", 0)
            # If Kraken2 did not run, everything non-host lands in one bucket.
            if not classified and not unclassified:
                unclassified = row.get("nonhost_reads", 0)
            handle.write(f"{row['sample']}\t{filtered}\t{rrna}\t{host}\t{classified}\t{unclassified}\n")


def write_mqc_generalstats(path, rows):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "# id: 'reanatax_host_carryover'\n"
            "# section_name: 'reanaTax host carry-over'\n"
            "# plot_type: 'generalstats'\n"
            "# pconfig:\n"
            "#     - host_carryover_pct:\n"
            "#         title: 'Host carry-over'\n"
            "#         description: 'Non-host reads that Kraken2 still assigns to the host taxon. High values mean host depletion is leaking.'\n"
            "#         suffix: '%'\n"
            "#         min: 0\n"
            "#         scale: 'OrRd'\n"
            "#         format: '{:,.2f}'\n"
            "#     - nonhost_pct_of_raw:\n"
            "#         title: 'Non-host'\n"
            "#         description: 'Share of the raw library that survived trimming and host depletion.'\n"
            "#         suffix: '%'\n"
            "#         min: 0\n"
            "#         max: 100\n"
            "#         scale: 'BuGn'\n"
            "#         format: '{:,.2f}'\n"
        )
        handle.write("Sample\thost_carryover_pct\tnonhost_pct_of_raw\n")
        for row in rows:
            handle.write(f"{row['sample']}\t{row['host_carryover_pct']}\t{row['nonhost_pct_of_raw']}\n")


def main():
    parser = argparse.ArgumentParser(description="Per-sample read accounting across the pipeline.")
    parser.add_argument("--fastp", nargs="*", default=[], help="fastp *.fastp.json files")
    parser.add_argument("--hisat2", nargs="*", default=[], help="HISAT2 *.hisat2.summary.log files")
    parser.add_argument("--sortmerna", nargs="*", default=[], help="SortMeRNA *.sortmerna.log files")
    parser.add_argument("--kraken2", nargs="*", default=[], help="Kraken2 *.report.txt files")
    parser.add_argument("--host-taxid", type=int, default=HOST_DEFAULT_TAXID,
                        help="taxid counted as host carry-over (default: 9606, Homo sapiens)")
    parser.add_argument("--prefix", default="reanatax", help="output file prefix")
    args = parser.parse_args()

    rows = [derive(row) for row in collect(args)]
    if not rows:
        print("[read_accounting] no inputs matched; nothing written.", file=sys.stderr)
        return 0

    write_tsv(f"{args.prefix}.read_accounting.tsv", rows)
    write_mqc_bargraph(f"{args.prefix}_read_fate_mqc.tsv", rows)
    write_mqc_generalstats(f"{args.prefix}_host_carryover_mqc.tsv", rows)
    print(f"[read_accounting] wrote {len(rows)} sample(s).", file=sys.stderr)

    worst = max(rows, key=lambda r: r["host_carryover_pct"])
    if worst["host_carryover_pct"] >= 10:
        print(
            f"[read_accounting] WARNING: {worst['sample']} has {worst['host_carryover_pct']}% host carry-over "
            f"after depletion. Check that the host reference matches the organism, and consider adding a "
            f"second assembly (see --fasta / --host_accession, which accept two references).",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())

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
    return {
        "raw_reads": (data.get("summary", {}).get("before_filtering", {}) or {}).get("total_reads", 0),
        "trimmed_reads": (data.get("summary", {}).get("after_filtering", {}) or {}).get("total_reads", 0),
        "trimmed_mean_length": (data.get("summary", {}).get("after_filtering", {}) or {}).get("read1_mean_length", 0),
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
    }


def read_kraken2(path, host_taxid):
    """Kraken2 report: column 2 is reads covered by the clade rooted at that
    taxon, column 3 reads assigned directly to it, column 5 the taxid.
    `--report-minimizer-data` inserts two extra columns, so index from the
    rank-code column rather than assuming a fixed layout."""
    classified = unclassified = host = 0
    try:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 5:
                    continue
                # Locate the rank-code column; taxid is the one after it.
                rank_idx = next(
                    (i for i, f in enumerate(fields) if re.fullmatch(r"[URDKPCOFGS][0-9]*", f.strip())),
                    None,
                )
                if rank_idx is None or rank_idx + 1 >= len(fields):
                    continue
                try:
                    clade_reads = int(fields[1])
                    taxid = int(fields[rank_idx + 1])
                except ValueError:
                    continue
                rank = fields[rank_idx].strip()
                if rank == "U":
                    unclassified = clade_reads
                elif rank == "R":
                    classified = clade_reads
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

    for path in args.kraken2:
        slot(sample_name(path, [".kraken2.report.txt", ".kraken2.report", ".report.txt"])).update(
            read_kraken2(path, args.host_taxid)
        )

    return [rows[key] for key in sorted(rows)]


def pass_order(label):
    """Depletion passes in the order they ran. Intermediate passes are labelled
    `host1`, `host2`, ...; the FINAL pass is always plain `host`, because that
    is the one whose leftovers get classified. Sorting alphabetically would put
    `host` before `host1` and silently report the wrong pass as final."""
    match = re.fullmatch(r"host(\d+)", label)
    return (0, int(match.group(1))) if match else (1, 0)


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
            "# description: 'What happened to every raw read: removed by fastp, aligned to a host\n"
            "#     reference, or carried through to classification. Categories are disjoint, so the\n"
            "#     bar length is the raw library size.'\n"
            "# plot_type: 'bargraph'\n"
            "# pconfig:\n"
            "#     id: 'reanatax_read_fate_plot'\n"
            "#     title: 'reanaTax: read fate'\n"
            "#     ylab: 'Reads'\n"
            "# section_href: 'https://github.com/BioinfoIPBLN/reanatax'\n"
        )
        handle.write("Sample\tFiltered by fastp\tAligned to host\tNon-host, classified\tNon-host, unclassified\n")
        for row in rows:
            raw = row.get("raw_reads", 0)
            trimmed = row.get("trimmed_reads", 0)
            filtered = max(raw - trimmed, 0)
            host = row.get("host_removed_reads", 0)
            classified = row.get("classified_reads", 0)
            unclassified = row.get("unclassified_reads", 0)
            # If Kraken2 did not run, everything non-host lands in one bucket.
            if not classified and not unclassified:
                unclassified = row.get("nonhost_reads", 0)
            handle.write(f"{row['sample']}\t{filtered}\t{host}\t{classified}\t{unclassified}\n")


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

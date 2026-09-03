#!/usr/bin/env python3
"""
host_kmer_filter.py -- how much of a taxon's evidence is really host sequence?

Host depletion is an alignment problem and alignment is not exhaustive. Whatever
HISAT2 fails to place stays in the non-host fraction, gets classified, and is
reported as a microbe. The abundance and minimizer filters cannot see this: a
carried-over host read is a genuine read, it has genuine k-mers, and if enough
of them accumulate on some taxon that taxon looks well supported.

SAHMI (Ghaddar, Blaser & De, Nat Comput Sci 2023) attacks it one level down. Its
sckmer step discards any read carrying EVEN ONE k-mer assigned to the host,
regardless of what the read as a whole was classified as. That is much stricter
than dropping reads whose final assignment is the host - a chimeric or repetitive
read can be assigned to a bacterium while most of its k-mers are host - and it is
the mechanism behind the carry-over rather than its symptom.

This applies the same test to bulk data, which is where the pipeline can already
supply what it needs: Kraken2's --output file (--kraken2_save_readclassifications)
records, per read, the taxon each run of k-mers was assigned to. Reads are not
rewritten; what comes out is the accounting. For every taxon: how many reads it
was given, and how many of those carried host k-mers. A taxon whose reads are
mostly host-tainted is host leakage wearing a species name, and its taxid goes on
the drop list the abundance filter already consumes.

The one prerequisite is that the host must be IN the Kraken2 database, or no
k-mer can ever be assigned to it and every count here is zero. core_nt contains
Anopheles gambiae (taxid 7165); PlusPF does not. The filter refuses to write a
drop list when it sees no host k-mers at all, rather than silently reporting a
clean run.

Two passes, because the --output files are one line per read and large:
  scan   one task per sample, streaming, writes a small per-taxon table
  merge  one task for the cohort, sums the tables and makes the drop list
"""
import argparse
import csv
import os
import re
import sys

# `Name (taxid 1234)` when Kraken2 ran with --use-names, a bare integer without.
NAMED_TAXID = re.compile(r"^(?P<name>.*)\(taxid\s+(?P<taxid>\d+)\)\s*$")


def parse_assignment(field):
    """(taxid, name) from the third column of a Kraken2 --output line."""
    text = field.strip()
    match = NAMED_TAXID.match(text)
    if match:
        return int(match.group("taxid")), match.group("name").strip()
    try:
        return int(text), ""
    except ValueError:
        return None, ""


def scan(path, host_taxids):
    """{taxid: [reads, host_tainted, name]} for one sample, plus totals.

    The host test is a substring match on the k-mer string, which is how SAHMI
    does it and is far cheaper than tokenising every read. It is exact because
    of the shape of the string: k-mer runs are written `taxid:count`, separated
    by single spaces, with ` |:| ` between mates. So a host taxid can only ever
    appear as `<host>:` at the very start of the string or as ` <host>:` after a
    separator - a count can never be mistaken for it, and 7165 cannot match
    inside 17165.
    """
    prefixes = tuple(f"{taxid}:" for taxid in host_taxids)
    needles = tuple(f" {taxid}:" for taxid in host_taxids)

    per_taxon = {}
    reads_total = 0
    reads_tainted = 0

    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 4)
            if len(fields) < 5:
                continue
            taxid, name = parse_assignment(fields[2])
            # Taxid 0 is the unclassified bucket: a bookkeeping row, not a
            # taxon, and counting it would put unclassified reads in the
            # denominator of a "share of classified reads" figure.
            if not taxid:
                continue
            kmers = fields[4]
            tainted = kmers.startswith(prefixes) or any(n in kmers for n in needles)

            reads_total += 1
            reads_tainted += tainted
            entry = per_taxon.get(taxid)
            if entry is None:
                per_taxon[taxid] = [1, int(tainted), name]
            else:
                entry[0] += 1
                entry[1] += tainted
                if name and not entry[2]:
                    entry[2] = name
    return per_taxon, reads_total, reads_tainted


def write_scan(path, sample, per_taxon, reads_total, reads_tainted):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(f"# sample\t{sample}\n")
        handle.write(f"# reads\t{reads_total}\n")
        handle.write(f"# reads_with_host_kmers\t{reads_tainted}\n")
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "name", "reads", "host_tainted"])
        for taxid in sorted(per_taxon, key=lambda t: -per_taxon[t][0]):
            reads, tainted, name = per_taxon[taxid]
            writer.writerow([taxid, name, reads, tainted])


def read_scan(path):
    """(sample, reads, tainted, {taxid: [reads, tainted, name]}) from one table."""
    sample = os.path.basename(path)
    totals = [0, 0]
    per_taxon = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                parts = line.lstrip("#").strip().split("\t")
                if len(parts) == 2 and parts[0] == "sample":
                    sample = parts[1]
                elif len(parts) == 2 and parts[0] == "reads":
                    totals[0] = int(parts[1])
                elif len(parts) == 2 and parts[0] == "reads_with_host_kmers":
                    totals[1] = int(parts[1])
                continue
            fields = line.split("\t")
            if len(fields) < 4 or fields[0] == "taxid":
                continue
            per_taxon[int(fields[0])] = [int(fields[2]), int(fields[3]), fields[1]]
    return sample, totals[0], totals[1], per_taxon


def merge(paths, args):
    combined = {}
    per_sample_totals = []
    for path in paths:
        sample, reads, tainted, per_taxon = read_scan(path)
        per_sample_totals.append((sample, reads, tainted))
        for taxid, (taxon_reads, taxon_tainted, name) in per_taxon.items():
            entry = combined.setdefault(taxid, [0, 0, name, 0])
            entry[0] += taxon_reads
            entry[1] += taxon_tainted
            entry[3] += 1
            if name and not entry[2]:
                entry[2] = name
    return combined, per_sample_totals


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reads", nargs="+", required=True,
                        help="scan: one Kraken2 --output file. merge: the per-sample tables.")
    parser.add_argument("--merge", action="store_true",
                        help="combine per-sample tables instead of scanning reads")
    parser.add_argument("--sample", default="sample", help="scan: the sample id to record")
    parser.add_argument("--host-taxid", default="9606",
                        help="comma-separated taxids counted as host. The host species alone by "
                             "default; add its genus or family to catch k-mers Kraken2 could only "
                             "place higher up, at the cost of condemning their real relatives too.")
    parser.add_argument("--output", required=True, help="table to write")
    parser.add_argument("--drop-list", help="merge: taxids to drop, one per line")
    parser.add_argument("--mqc", help="merge: MultiQC custom-content section to write")
    parser.add_argument("--max-host-fraction", type=float, default=0.5,
                        help="merge: drop a taxon once this fraction of its reads carry host k-mers")
    parser.add_argument("--min-reads", type=int, default=10,
                        help="merge: leave taxa below this many reads to the abundance filter")
    args = parser.parse_args()

    host_taxids = [int(part) for part in str(args.host_taxid).split(",") if part.strip()]
    if not host_taxids:
        raise SystemExit("host_kmer_filter: --host-taxid is required and must be numeric.")

    if not args.merge:
        per_taxon, reads_total, reads_tainted = scan(args.reads[0], host_taxids)
        write_scan(args.output, args.sample, per_taxon, reads_total, reads_tainted)
        share = (reads_tainted / reads_total * 100) if reads_total else 0.0
        print(
            f"[host_kmer_filter] {args.sample}: {reads_tainted}/{reads_total} classified reads "
            f"({share:.3f}%) carry k-mers of taxid {'/'.join(map(str, host_taxids))}.",
            file=sys.stderr,
        )
        return 0

    combined, per_sample_totals = merge(args.reads, args)

    tainted_anywhere = sum(entry[1] for entry in combined.values())
    if tainted_anywhere == 0:
        raise SystemExit(
            f"host_kmer_filter: not one read carries a k-mer of taxid "
            f"{'/'.join(map(str, host_taxids))}. That is a database question, not a clean result: "
            "the test is vacuous unless the host is present in the Kraken2 database (core_nt has "
            "Anopheles gambiae, PlusPF does not). Check --host_kmer_taxid against the database, or "
            "turn the filter off."
        )

    dropped = []
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxid", "name", "samples", "reads", "host_tainted",
                         "host_fraction", "verdict"])
        for taxid in sorted(combined, key=lambda t: -combined[t][0]):
            reads, tainted, name, samples = combined[taxid]
            fraction = tainted / reads if reads else 0.0
            if taxid in host_taxids:
                call = "host"
            elif reads < args.min_reads:
                call = "below_read_gate"
            elif fraction >= args.max_host_fraction:
                call = "host_kmer_leakage"
                dropped.append(taxid)
            else:
                call = "kept"
            writer.writerow([taxid, name, samples, reads, tainted, f"{fraction:.6f}", call])

    if args.drop_list:
        with open(args.drop_list, "w", encoding="utf-8") as handle:
            for taxid in dropped:
                handle.write(f"{taxid}\n")

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_host_kmer'\n"
                "# section_name: 'Host k-mer carry-over'\n"
                "# description: 'Share of each library's CLASSIFIED reads carrying at least one\n"
                "#     k-mer assigned to the host, whatever the read itself was classified as.\n"
                "#     Alignment-based depletion cannot remove these and an abundance filter cannot\n"
                "#     see them, because they are real reads with real k-mers. Requires the host to\n"
                "#     be present in the Kraken2 database.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_host_kmer_plot'\n"
                "#     title: 'reanaTax: reads carrying host k-mers'\n"
                "#     ylab: 'Classified reads'\n"
                "Sample\tWith host k-mers\tWithout\n"
            )
            for sample, reads, tainted in sorted(per_sample_totals):
                handle.write(f"{sample}\t{tainted}\t{max(reads - tainted, 0)}\n")

    print(
        f"[host_kmer_filter] {len(combined)} taxa over {len(per_sample_totals)} samples; "
        f"{tainted_anywhere} host-tainted read assignments; "
        f"{len(dropped)} taxid(s) listed for removal at >= {args.max_host_fraction:g} tainted.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

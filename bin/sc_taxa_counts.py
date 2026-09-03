#!/usr/bin/env python3
"""
sc_taxa_counts.py -- a cell-by-taxon matrix from Kraken2's read-level output.

The single-cell branch classifies the reads STARsolo could not place on the
host, with the corrected cell barcode folded into each read name by
STARSOLO_UNMAPPED. This turns that into counts per (barcode, taxon), applying
the filters SAHMI's sckmer.r/taxa_counts.r apply and two it does not.

What is taken from SAHMI:

  host k-mer exclusion  a read carrying EVEN ONE k-mer assigned to the host is
                        dropped, whatever the read as a whole was called. Far
                        stricter than dropping reads assigned to the host, and
                        the mechanism behind carry-over rather than its symptom.
  min_frac              a read counts towards a taxon only if at least this
                        fraction of its k-mers fall inside that taxon's lineage.
                        Guards against a read scattered over the tree being
                        credited to whichever tip won the LCA.
  homopolymer filter    reads with a run of one nucleotide longer than
                        --max-homopolymer are dropped. SAHMI's nFilter, applied
                        to the read rather than to a whole-read base count.

Where this departs from SAHMI:

  positional barcodes   sckmer.r reads the barcode as substr(R1, 1, cb_len),
                        with no whitelist and no error correction, so one
                        sequencing error in the barcode manufactures a new
                        "cell". The barcode here is STARsolo's CB, already
                        corrected against the whitelist.
  UMI handling          SAHMI's taxa_counts.r DOES deduplicate - it builds
                        (barcode, umi, taxid) and takes unique() - so its matrix
                        is already UMI counts. What it does not do is carry that
                        through to sckmer.r, whose k-mer statistics, and so the
                        barcode-level denoising built on them, are computed over
                        undeduplicated reads. Here one (barcode, UMI, taxon)
                        triple is counted once throughout, so the matrix and the
                        statistics agree. --no-umi-dedup restores read counting.

Lineage comes from the Kraken2 report, whose leading indentation encodes the
tree: two spaces per level. No taxonomy dump is needed.
"""
import argparse
import csv
import gzip
import os
import re
import sys

NAMED_TAXID = re.compile(r"^(?P<name>.*)\(taxid\s+(?P<taxid>\d+)\)\s*$")
BARCODE_IN_NAME = re.compile(r"\|CB:(?P<cb>[^|]*)\|UB:(?P<ub>[^|]*)$")


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def parse_report(path):
    """(ancestors, names, ranks) from a Kraken2 report.

    `ancestors[taxid]` is the set of taxids on the path from the root down to
    and including that taxon. The report's leading indentation is two spaces per
    rank level, which is the tree; reconstructing it from there avoids shipping
    a taxonomy dump alongside every run.
    """
    ancestors, names, ranks = {}, {}, {}
    stack = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("%"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            raw_name = fields[-1]
            depth = (len(raw_name) - len(raw_name.lstrip(" "))) // 2
            # Layout-agnostic, as everywhere else here: Kraken2 ends
            # `rank taxid name`, KrakenUniq `taxid rank name`.
            try:
                taxid = int(fields[-2])
                rank = fields[-3]
            except ValueError:
                try:
                    taxid = int(fields[-3])
                    rank = fields[-2]
                except ValueError:
                    continue
            del stack[depth:]
            stack.append(taxid)
            ancestors[taxid] = set(stack)
            names[taxid] = raw_name.strip()
            ranks[taxid] = rank.strip()
    return ancestors, names, ranks


def parse_assignment(field):
    text = field.strip()
    match = NAMED_TAXID.match(text)
    if match:
        return int(match.group("taxid")), match.group("name").strip()
    try:
        return int(text), ""
    except ValueError:
        return None, ""


def kmer_runs(kmers):
    """[(taxid_or_None, count)] from Kraken2's per-read k-mer string."""
    runs = []
    for token in kmers.replace("|:|", " ").split():
        head, _, tail = token.rpartition(":")
        if not head:
            continue
        try:
            count = int(tail)
        except ValueError:
            continue
        try:
            runs.append((int(head), count))
        except ValueError:
            # `A` (ambiguous nucleotides) and any other non-numeric label.
            runs.append((None, count))
    return runs


def fastq_stream(path):
    """Yield (read_id, sequence) in file order."""
    with open_maybe_gzip(path) as handle:
        read_id = None
        for index, line in enumerate(handle):
            if index % 4 == 0:
                read_id = line[1:].split()[0] if len(line) > 1 else None
            elif index % 4 == 1 and read_id:
                yield read_id, line.strip().upper()


class SequenceReader:
    """Sequences for Kraken2 output lines, without holding the FASTQ in memory.

    Kraken2 writes one --output line per input sequence, in input order, so the
    two files advance together and the sequences never need to be indexed. The
    reader still tolerates a Kraken2 line having no counterpart - it scans
    forward a bounded distance and gives up rather than silently pairing a read
    with the wrong sequence, which is the failure that would matter.
    """

    def __init__(self, path, lookahead=4096):
        self.stream = fastq_stream(path)
        self.lookahead = lookahead
        self.exhausted = False
        self.misses = 0

    def get(self, read_id):
        if self.exhausted:
            return None
        for _ in range(self.lookahead):
            try:
                current_id, sequence = next(self.stream)
            except StopIteration:
                self.exhausted = True
                return None
            if current_id == read_id:
                return sequence
        self.misses += 1
        return None


def longest_homopolymer(sequence):
    longest = run = 0
    previous = ""
    for base in sequence:
        run = run + 1 if base == previous else 1
        previous = base
        if run > longest:
            longest = run
    return longest


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reads", required=True, help="Kraken2 --output file for this sample")
    parser.add_argument("--report", required=True, help="the matching Kraken2 report, for the lineage")
    parser.add_argument("--fastq", help="the classified FASTQ, only needed for --max-homopolymer")
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--output", required=True, help="long-format barcode x taxon table")
    parser.add_argument("--summary", required=True, help="per-barcode summary")
    parser.add_argument("--host-taxid", default="", help="comma-separated taxids to treat as host")
    parser.add_argument("--min-frac", type=float, default=0.5,
                        help="fraction of a read's k-mers that must fall in the taxon's lineage")
    parser.add_argument("--max-homopolymer", type=int, default=0,
                        help="drop reads with a longer single-nucleotide run (0 disables; needs --fastq)")
    parser.add_argument("--ranks", default="S,G", help="ranks to tabulate, comma-separated codes")
    parser.add_argument("--no-umi-dedup", action="store_true",
                        help="count reads rather than distinct UMIs (SAHMI's behaviour)")
    parser.add_argument("--min-umis", type=int, default=2,
                        help="UMIs a taxon needs in a cell to be called present there "
                             "(CSI-Microbes' rule; 1 keeps every observation)")
    parser.add_argument("--sweep", help="write a sensitivity sweep of --min-umis 1..5 here")
    args = parser.parse_args()

    host_taxids = {int(part) for part in args.host_taxid.split(",") if part.strip()}
    wanted_ranks = {part.strip().upper() for part in args.ranks.split(",") if part.strip()}
    ancestors, names, ranks = parse_report(args.report)

    if args.max_homopolymer and not args.fastq:
        raise SystemExit("sc_taxa_counts: --max-homopolymer needs --fastq to read the sequences.")
    # Streamed alongside the Kraken2 output rather than indexed: the two files
    # are in the same order, and a deep library's sequences do not fit in memory
    # comfortably. Only opened when a homopolymer limit is actually set.
    sequences = SequenceReader(args.fastq) if (args.max_homopolymer and args.fastq) else None

    # (barcode, taxid) -> set of UMIs, or a count when dedup is off.
    observed = {}
    stats = {
        "reads": 0, "no_barcode": 0, "unclassified": 0, "host_kmer": 0,
        "low_complexity": 0, "below_min_frac": 0, "wrong_rank": 0, "counted": 0,
    }

    with open_maybe_gzip(args.reads) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 4)
            if len(fields) < 5:
                continue
            stats["reads"] += 1
            taxid, _name = parse_assignment(fields[2])
            if not taxid:
                stats["unclassified"] += 1
                continue

            read_id = fields[1]
            match = BARCODE_IN_NAME.search(read_id)
            if not match or match.group("cb") == "NA":
                stats["no_barcode"] += 1
                continue
            barcode, umi = match.group("cb"), match.group("ub")

            if sequences is not None:
                sequence = sequences.get(read_id)
                if sequence and longest_homopolymer(sequence) > args.max_homopolymer:
                    stats["low_complexity"] += 1
                    continue

            runs = kmer_runs(fields[4])
            if host_taxids and any(run_taxid in host_taxids for run_taxid, _count in runs):
                stats["host_kmer"] += 1
                continue

            lineage = ancestors.get(taxid)
            if lineage is None:
                continue
            total = sum(count for _taxid, count in runs)
            inside = sum(count for run_taxid, count in runs if run_taxid in lineage)
            if not total or inside / total < args.min_frac:
                stats["below_min_frac"] += 1
                continue

            # Tabulate at the requested ranks: a read assigned below species
            # counts towards its species, and towards its genus.
            targets = [
                ancestor for ancestor in lineage
                if ranks.get(ancestor, "").upper() in wanted_ranks
                or ranks.get(ancestor, "")[:1].upper() in wanted_ranks and ranks.get(ancestor, "")[1:].isdigit()
            ]
            if not targets:
                stats["wrong_rank"] += 1
                continue

            stats["counted"] += 1
            for target in targets:
                key = (barcode, target)
                if args.no_umi_dedup:
                    observed[key] = observed.get(key, 0) + 1
                else:
                    observed.setdefault(key, set()).add(umi)

    counts = {
        key: (len(value) if isinstance(value, set) else value)
        for key, value in observed.items()
    }

    # The presence call. CSI-Microbes (Robinson et al., Sci Adv 2024) does not
    # test for presence at all - it thresholds, requiring at least two UMIs of a
    # taxon in a cell. On their positive control that threshold removed 100% of
    # contaminant genera at single-cell level, which is a stronger result than
    # any of the correlation tests here achieve. It is also, in their own words,
    # arbitrary (their stated limitation 2), so --sweep writes the same table at
    # thresholds 1-5 and no conclusion should rest on one of them.
    kept = {key: count for key, count in counts.items() if count >= args.min_umis}

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "barcode", "taxid", "rank", "name", "count"])
        for (barcode, taxid), count in sorted(kept.items(), key=lambda item: (item[0][0], -item[1])):
            writer.writerow([args.sample, barcode, taxid, ranks.get(taxid, ""), names.get(taxid, ""), count])

    if args.sweep:
        with open(args.sweep, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["sample", "min_umis", "cells", "taxa", "cell_taxon_pairs"])
            for threshold in range(1, 6):
                surviving = {key for key, count in counts.items() if count >= threshold}
                writer.writerow([
                    args.sample, threshold,
                    len({barcode for barcode, _taxid in surviving}),
                    len({taxid for _barcode, taxid in surviving}),
                    len(surviving),
                ])

    barcodes = {barcode for barcode, _taxid in kept}
    with open(args.summary, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "metric", "value"])
        for metric, value in stats.items():
            writer.writerow([args.sample, metric, value])
        writer.writerow([args.sample, "barcodes_with_taxa", len(barcodes)])
        writer.writerow([args.sample, "taxa", len({taxid for _barcode, taxid in kept})])
        writer.writerow([args.sample, "min_umis", args.min_umis])
        writer.writerow([args.sample, "pairs_below_min_umis", len(counts) - len(kept)])

    print(
        f"[sc_taxa_counts] {args.sample}: {stats['reads']} classified reads -> "
        f"{stats['counted']} counted over {len(barcodes)} barcode(s); dropped "
        f"{stats['no_barcode']} no-barcode, {stats['host_kmer']} host-k-mer, "
        f"{stats['below_min_frac']} below min-frac, {stats['low_complexity']} low-complexity.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

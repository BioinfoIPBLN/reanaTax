#!/usr/bin/env python3
"""
sc_kmer_denoise.py -- does a taxon's evidence scale across CELLS?

The bulk pipeline already asks this across samples (--minimizer_correlation).
SAHMI's single-cell contribution is to ask it across barcodes within one sample,
which is a much finer grain: a library has a handful of samples but thousands of
cells, so the correlation has real power where the sample-level version needs at
least five libraries to say anything at all.

The logic is the same. For a real organism, a barcode carrying more of its reads
carries proportionally more DISTINCT k-mers too, because the reads sample more
of its genome. For ambient contamination smeared over every droplet, or for
reads piling on one conserved locus, the total rises and the distinct count
saturates. SAHMI requires the Spearman correlation between the two, taken over
barcodes, to be significant, on more than three barcodes, with p adjusted by
HOLM (FWER) rather than BH.

Size the expectation correctly. On SAHMI's own shipped example
(SRR9713132.sckmer.txt) 98.5% of barcode-taxon rows have total k-mers EXACTLY
equal to distinct k-mers, and in their published PDAC scoring 94.8% of testable
rows clear p < 0.05. The correlation is therefore close to a tautology on that
data, and almost all of the filtering power comes from the requirement that a
taxon be seen on at least four barcodes carrying more than one k-mer. Treat this
as a prevalence filter with a significance test attached, not as an independent
statistical test, and do not read a high pass rate as validation.

Distinct k-mers cannot be read out of the Kraken2 report here - that is a
per-sample number, and the question is per barcode. They are recovered from the
reads themselves, which is what SAHMI's sckmer.r does:

  Kraken2's --output records, per read, a run-length encoding of which taxon
  each consecutive k-mer matched: `taxid:count taxid:count ...` in read order.
  A read of length L yields L-k+1 k-mers, so run i covers a known span of
  positions and the k-mers assigned to a taxon can be pulled straight out of
  the sequence.

Each k-mer is stored as a 2-bit-per-base integer rather than a string: exact
(no hash collisions), deterministic across resumes, and about an order of
magnitude smaller in memory than keeping the sequences. k-mers containing any
non-ACGT base are skipped, as Kraken2 could not have matched them either.
"""
import argparse
import csv
import gzip
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import benjamini_hochberg, holm, spearman

BARCODE_IN_NAME = re.compile(r"\|CB:(?P<cb>[^|]*)\|UB:(?P<ub>[^|]*)$")
NAMED_TAXID = re.compile(r"\(taxid\s+(?P<taxid>\d+)\)\s*$")
BASES = {"A": 0, "C": 1, "G": 2, "T": 3}


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path, encoding="utf-8", errors="replace")


def pack_read(sequence):
    """(2-bit packing of the whole read, index of the last non-ACGT base at or
    before each position).

    Encoding each k-mer on its own re-walks all 35 of its bases, so a read
    whose runs cover it end to end pays thousands of base lookups where the
    read holds only 90 bases. Walking it once turns a k-mer into a bit slice,
    and a window is valid exactly when the last non-ACGT base falls before it.
    """
    packed = 0
    prev_bad = []
    last_bad = -1
    for index, base in enumerate(sequence):
        code = BASES.get(base)
        if code is None:
            code = 0
            last_bad = index
        packed = (packed << 2) | code
        prev_bad.append(last_bad)
    return packed, prev_bad


def parse_assignment(field):
    text = field.strip()
    match = NAMED_TAXID.search(text)
    if match:
        return int(match.group("taxid"))
    try:
        return int(text)
    except ValueError:
        return None


def kmer_runs(kmers):
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


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reads", required=True, help="Kraken2 --output file")
    parser.add_argument("--fastq", required=True, help="the reads that were classified")
    parser.add_argument("--sample", default="sample")
    parser.add_argument("--evidence", required=True, help="per-taxon evidence table")
    parser.add_argument("--drop-list", required=True, help="taxids that failed, one per line")
    parser.add_argument("--mqc", help="MultiQC custom-content section")
    parser.add_argument("--kmer-len", type=int, default=35, help="Kraken2's k (default database: 35)")
    parser.add_argument("--min-barcodes", type=int, default=4,
                        help="a taxon needs more than three barcodes to be tested (SAHMI's rule)")
    parser.add_argument("--min-kmers", type=int, default=2,
                        help="a barcode counts towards the test only above this many k-mers")
    parser.add_argument("--max-barcodes-per-taxon", type=int, default=1000,
                        help="cap on barcodes held per taxon, SAHMI's nsample; 0 removes the cap")
    parser.add_argument("--correlation-p", type=float, default=0.05,
                        help="adjusted p the correlation must clear")
    parser.add_argument("--adjust", choices=["holm", "bh"], default="holm",
                        help="multiple-testing correction; SAHMI uses holm (FWER)")
    parser.add_argument("--host-taxid", default="", help="comma-separated taxids to treat as host")
    args = parser.parse_args()

    host_taxids = {int(part) for part in args.host_taxid.split(",") if part.strip()}
    sequences = SequenceReader(args.fastq)

    # (taxid, barcode) -> [total k-mers, {distinct encoded k-mers}]
    per_barcode = {}
    held_barcodes = {}
    dropped_host = 0
    kmer_mask = (1 << (2 * args.kmer_len)) - 1

    with open_maybe_gzip(args.reads) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t", 4)
            if len(fields) < 5:
                continue
            taxid = parse_assignment(fields[2])
            if not taxid:
                continue
            read_id = fields[1]
            match = BARCODE_IN_NAME.search(read_id)
            if not match or match.group("cb") == "NA":
                continue
            barcode = match.group("cb")

            runs = kmer_runs(fields[4])
            if host_taxids and any(run_taxid in host_taxids for run_taxid, _count in runs):
                dropped_host += 1
                continue

            sequence = sequences.get(read_id)
            if not sequence:
                continue

            # Packed lazily: a read assigned to an ancestor by LCA can have no
            # run matching its own taxid, and then nothing here is encoded.
            packed = None
            prev_bad = None
            length = len(sequence)

            # Walk the runs in read order; run i starts at k-mer index `offset`.
            offset = 0
            for run_taxid, count in runs:
                if run_taxid == taxid:
                    key = (taxid, barcode)
                    entry = per_barcode.get(key)
                    if entry is None:
                        # SAHMI's nsample: hold at most this many barcodes per
                        # taxon. Without a cap an abundant taxon in a deep
                        # library would hold a k-mer set for every droplet.
                        # Counted, not recounted: entries are never removed,
                        # so the running tally is exactly what scanning the
                        # dict used to produce - at O(1) instead of O(n) per
                        # new barcode, which was quadratic over the file.
                        if args.max_barcodes_per_taxon and held_barcodes.get(taxid, 0) >= args.max_barcodes_per_taxon:
                            offset += count
                            continue
                        held_barcodes[taxid] = held_barcodes.get(taxid, 0) + 1
                        entry = per_barcode[key] = [0, set()]
                    entry[0] += count
                    if packed is None:
                        packed, prev_bad = pack_read(sequence)
                    for start in range(offset, offset + count):
                        stop = start + args.kmer_len
                        if stop > length:
                            break
                        if prev_bad[stop - 1] >= start:
                            continue
                        entry[1].add((packed >> (2 * (length - stop))) & kmer_mask)
                offset += count

    # Regroup by taxon, then correlate total against distinct over barcodes.
    by_taxon = {}
    for (taxid, barcode), (total, distinct) in per_barcode.items():
        if total <= args.min_kmers or len(distinct) <= 1:
            continue
        by_taxon.setdefault(taxid, []).append((total, len(distinct)))

    tested, raw, saturated = [], [], []
    for taxid, observations in sorted(by_taxon.items()):
        if len(observations) < args.min_barcodes:
            continue
        totals = [total for total, _distinct in observations]
        distincts = [distinct for _total, distinct in observations]
        rho, pvalue = spearman(totals, distincts)
        if rho is None:
            # Spearman is undefined when a vector never varies. Which vector
            # decides what that means, and the two cases are opposite:
            #
            #   distinct constant while total rises  the taxon accumulates reads
            #       across cells without ever covering anything new. That is not
            #       an untestable taxon, it is the contamination signature in its
            #       purest form, and it fails outright.
            #   total constant  nothing varied, so nothing was measured. Left
            #       alone, like the sample-level test does.
            if len(set(distincts)) == 1 and len(set(totals)) > 1:
                saturated.append((taxid, len(observations)))
            continue
        tested.append((taxid, len(observations), rho))
        raw.append(pvalue)

    adjusted = holm(raw) if args.adjust == 'holm' else benjamini_hochberg(raw)
    dropped = []
    with open(args.evidence, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", "taxid", "barcodes", "rho_kmers_distinct", "p_adj", "verdict"])
        for taxid, barcodes in saturated:
            dropped.append(taxid)
            writer.writerow([args.sample, taxid, barcodes, "NA", "NA", "saturated_distinct"])
        for (taxid, barcodes, rho), qvalue in zip(tested, adjusted):
            passed = rho > 0 and qvalue is not None and qvalue <= args.correlation_p
            if not passed:
                dropped.append(taxid)
            writer.writerow([
                args.sample, taxid, barcodes, f"{rho:.4f}",
                "NA" if qvalue is None else f"{qvalue:.6g}",
                "kept" if passed else "no_barcode_correlation",
            ])

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in dropped:
            handle.write(f"{taxid}\n")

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_sc_kmer'\n"
                "# section_name: 'Barcode-level k-mer evidence'\n"
                "# description: 'Taxa judged on whether their evidence SCALES across cells: a barcode\n"
                "#     carrying more of a real organism carries proportionally more distinct k-mers too,\n"
                "#     while ambient contamination smeared over every droplet, or reads piling on one\n"
                "#     conserved locus, saturates. Spearman over barcodes, BH-adjusted. Taxa on three or\n"
                "#     fewer barcodes are not tested.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_sc_kmer_plot'\n"
                "#     title: 'reanaTax: barcode-level k-mer evidence'\n"
                "#     ylab: 'Taxa'\n"
                "Sample\tkept\tno_barcode_correlation\n"
                f"{args.sample}\t{len(tested) + len(saturated) - len(dropped)}\t{len(dropped)}\n"
            )

    print(
        f"[sc_kmer_denoise] {args.sample}: {len(by_taxon)} taxa over "
        f"{len({barcode for _taxid, barcode in per_barcode})} barcode(s); "
        f"{len(tested) + len(saturated)} testable (>{args.min_barcodes - 1} barcodes), "
        f"{len(tested) + len(saturated) - len(dropped)} passed, {len(dropped)} listed for removal"
        + (f"; {dropped_host} read(s) dropped for host k-mers." if host_taxids else "."),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
polya_carryover.py -- is the microbial signal in a poly(A) library real capture,
or physical carry-over?

Most bacterial transcripts have no poly(A) tail, so their presence in a
poly(A)-selected mRNA-seq library needs an explanation. There are two:

  * INTERNAL MISPRIMING - the oligo-dT primer catches an A-rich stretch inside
    the transcript. Reads captured this way are enriched for internal poly-A/T
    runs relative to what their own base composition would predict.
  * NON-SPECIFIC CARRY-OVER - the fragment came along for the ride, with no
    poly-A involvement at all. Those reads look compositionally ordinary.

Distinguishing them matters because it says how the abundances should be read:
mispriming biases towards A-rich transcripts, carry-over does not.

The comparison must be against the right null. A/T runs are common in AT-rich
genomes for trivial reasons, so counting them alone says nothing. Following
Monteleone et al. (Microbiome 2026), each read is shuffled to build an empirical
baseline that PRESERVES ITS LENGTH AND BASE COMPOSITION, and the observed rate is
compared against that. An enrichment near 1.0 means the runs are exactly what
composition predicts - carry-over. Substantially above 1.0 means capture.

Reads are streamed and sampled, so this runs in seconds on a full FASTQ.
"""
import argparse
import gzip
import io
import random
import re
import sys

DEFAULT_RUNS = (8, 10, 12)


def open_maybe_gzip(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def iter_reads(paths, max_reads):
    """Yield sequences from FASTQ files, stopping once max_reads is reached."""
    seen = 0
    for path in paths:
        try:
            with open_maybe_gzip(path) as handle:
                for index, line in enumerate(handle):
                    if index % 4 != 1:
                        continue
                    seq = line.strip().upper()
                    if not seq:
                        continue
                    yield seq
                    seen += 1
                    if max_reads and seen >= max_reads:
                        return
        except OSError as exc:
            print(f"[polya] cannot read {path}: {exc}", file=sys.stderr)


def has_run(seq, length, patterns):
    return bool(patterns[length].search(seq))


def merge(paths, output):
    """Fold the per-sample TSVs into one MultiQC line plot: enrichment against
    run length, one series per sample. A flat line at 1 is carry-over."""
    rows = []
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[0] == "sample":
                    continue
                rows.append(fields)
    lengths = sorted({int(r[1]) for r in rows})
    by_sample = {}
    for sample, run_length, _reads, _obs, _exp, enrichment in rows:
        value = float(enrichment)
        by_sample.setdefault(sample, {})[int(run_length)] = value

    with open(output, "w", encoding="utf-8") as handle:
        handle.write(
            "# id: 'reanatax_polya_carryover'\n"
            "# section_name: 'poly(A) carry-over check'\n"
            "# description: 'Internal poly-A/T runs in the non-host reads, relative to a null built by\n"
            "#     shuffling each read (preserving its length and base composition). Around 1 means the runs\n"
            "#     are what composition alone predicts, i.e. the microbial reads are non-specific carry-over.\n"
            "#     Clearly above 1 means oligo-dT internal mispriming captured them, which biases abundances\n"
            "#     towards A-rich transcripts. Only meaningful for poly(A)-selected libraries.'\n"
            "# plot_type: 'linegraph'\n"
            "# pconfig:\n"
            "#     id: 'reanatax_polya_carryover_plot'\n"
            "#     title: 'reanaTax: internal poly-A/T enrichment'\n"
            "#     xlab: 'Run length (nt)'\n"
            "#     ylab: 'Observed / expected'\n"
            "#     ymin: 0\n"
        )
        handle.write("Sample\t" + "\t".join(str(n) for n in lengths) + "\n")
        for sample in sorted(by_sample):
            handle.write(sample + "\t" + "\t".join(f"{by_sample[sample].get(n, 0):.3f}" for n in lengths) + "\n")
    print(f"[polya] merged {len(by_sample)} sample(s) into {output}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Internal poly-A/T enrichment vs a composition-matched null.")
    parser.add_argument("--fastq", nargs="*", default=[], help="non-host FASTQ file(s)")
    parser.add_argument("--merge", nargs="*", default=[], help="per-sample TSVs to fold into a MultiQC section")
    parser.add_argument("--sample", help="sample id (required unless --merge)")
    parser.add_argument("--output", required=True, help="TSV to write")
    parser.add_argument("--runs", type=int, nargs="+", default=list(DEFAULT_RUNS),
                        help="poly-A/T run lengths to test (default: 8 10 12)")
    parser.add_argument("--max-reads", type=int, default=200000,
                        help="reads to sample; 0 for all (default: 200000)")
    parser.add_argument("--permutations", type=int, default=1,
                        help="shuffles per read for the null (default: 1)")
    parser.add_argument("--seed", type=int, default=1,
                        help="fixed so the result is reproducible across runs")
    args = parser.parse_args()

    if args.merge:
        merge(args.merge, args.output)
        return 0
    if not args.fastq or not args.sample:
        parser.error("--fastq and --sample are required unless --merge is given")

    patterns = {n: re.compile(f"A{{{n},}}|T{{{n},}}") for n in args.runs}
    rng = random.Random(args.seed)

    observed = {n: 0 for n in args.runs}
    permuted = {n: 0 for n in args.runs}
    total = 0
    gc_total = 0.0

    for seq in iter_reads(args.fastq, args.max_reads):
        total += 1
        bases = list(seq)
        gc = sum(1 for b in bases if b in "GC")
        gc_total += gc / len(bases) if bases else 0.0
        for n in args.runs:
            if has_run(seq, n, patterns):
                observed[n] += 1
        # Shuffling in place preserves length and exact base composition (and so
        # GC) while destroying any positional structure -- which is the whole
        # point: only runs that survive shuffling are compositional accidents.
        for _ in range(args.permutations):
            rng.shuffle(bases)
            shuffled = "".join(bases)
            for n in args.runs:
                if has_run(shuffled, n, patterns):
                    permuted[n] += 1

    if total == 0:
        print(f"[polya] {args.sample}: no reads found; nothing written.", file=sys.stderr)
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write("sample\trun_length\treads\tobserved_pct\texpected_pct\tenrichment\n")
        return 0

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("sample\trun_length\treads\tobserved_pct\texpected_pct\tenrichment\n")
        for n in args.runs:
            draws = total * max(args.permutations, 1)
            obs = 100.0 * observed[n] / total
            exp = 100.0 * permuted[n] / draws
            # Half-count smoothing on the permuted side. A long run can be
            # absent from every shuffle, and an unsmoothed ratio would then be
            # infinite -- which has to be rendered as *something* in a plot, and
            # every choice (0, a cap) misreads as its opposite. The smoothed
            # ratio stays finite and monotonic in the observed count.
            exp_smoothed = 100.0 * (permuted[n] + 0.5) / (draws + 0.5)
            enrichment = obs / exp_smoothed
            handle.write(f"{args.sample}\t{n}\t{total}\t{obs:.4f}\t{exp:.4f}\t{enrichment:.3f}\n")

    mean_gc = 100.0 * gc_total / total
    summary = ", ".join(
        f"{n}nt obs={100.0 * observed[n] / total:.3f}% exp={100.0 * permuted[n] / (total * max(args.permutations, 1)):.3f}%"
        for n in args.runs
    )
    print(f"[polya] {args.sample}: {total} reads, mean GC {mean_gc:.1f}%; {summary}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

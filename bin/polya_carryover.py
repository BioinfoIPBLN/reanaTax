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

Two definitions of "run" are counted, because they answer different questions
and the literature uses the word for both:

  * MIXED (`[AT]{n,}`) - n consecutive positions each of which is A or T. This
    is the definition behind the ~5.5-6% background Monteleone et al. report at
    n=8, and it is the one to use when reproducing their figures.
  * HOMOPOLYMER (`A{n,}|T{n,}`) - n consecutive identical bases. Two orders of
    magnitude rarer, and the mechanistically sharper test: oligo-dT anneals to a
    stretch of A, not to an AT-rich stretch.

Reads are streamed and sampled, so this runs in seconds on a full FASTQ.
"""
import argparse
import gzip
import io
import random
import re
import sys

DEFAULT_RUNS = (8, 10, 12)
TRACTS = {"mixed": "[AT]{{{n},}}", "homopolymer": "A{{{n},}}|T{{{n},}}"}
HEADER = "sample\ttract\trun_length\treads\tobserved_pct\texpected_pct\tenrichment\n"


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


def has_run(seq, key, patterns):
    return bool(patterns[key].search(seq))


def merge(paths, output):
    """Fold the per-sample TSVs into one MultiQC line plot: enrichment against
    run length, one series per sample. A flat line at 1 is carry-over."""
    rows = []
    for path in paths:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7 or fields[0] == "sample":
                    continue
                rows.append(fields)
    lengths = sorted({int(r[2]) for r in rows})
    by_sample = {}
    for sample, tract, run_length, _reads, _obs, _exp, enrichment in rows:
        series = sample if tract == "mixed" else f"{sample} [{tract}]"
        by_sample.setdefault(series, {})[int(run_length)] = float(enrichment)

    with open(output, "w", encoding="utf-8") as handle:
        handle.write(
            "# id: 'reanatax_polya_carryover'\n"
            "# section_name: 'poly(A) carry-over check'\n"
            "# description: 'Internal poly-A/T runs in the non-host reads, relative to a null built by\n"
            "#     shuffling each read (preserving its length and base composition). Around 1 means the runs\n"
            "#     are what composition alone predicts, i.e. the microbial reads are non-specific carry-over.\n"
            "#     Unlabelled series count mixed A/T tracts; those marked [homopolymer] require identical bases,\n"
            "#     which is the stricter test for oligo-dT internal priming.\n"
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
    parser.add_argument("--tracts", nargs="+", choices=sorted(TRACTS), default=sorted(TRACTS),
                        help="tract definition(s) to count (default: both)")
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

    keys = [(tract, n) for tract in args.tracts for n in args.runs]
    patterns = {(tract, n): re.compile(TRACTS[tract].format(n=n)) for tract, n in keys}
    rng = random.Random(args.seed)

    observed = {key: 0 for key in keys}
    permuted = {key: 0 for key in keys}
    total = 0
    gc_total = 0.0

    for seq in iter_reads(args.fastq, args.max_reads):
        total += 1
        bases = list(seq)
        gc = sum(1 for b in bases if b in "GC")
        gc_total += gc / len(bases) if bases else 0.0
        for key in keys:
            if has_run(seq, key, patterns):
                observed[key] += 1
        # Shuffling in place preserves length and exact base composition (and so
        # GC) while destroying any positional structure -- which is the whole
        # point: only runs that survive shuffling are compositional accidents.
        for _ in range(args.permutations):
            rng.shuffle(bases)
            shuffled = "".join(bases)
            for key in keys:
                if has_run(shuffled, key, patterns):
                    permuted[key] += 1

    if total == 0:
        print(f"[polya] {args.sample}: no reads found; nothing written.", file=sys.stderr)
        with open(args.output, "w", encoding="utf-8") as handle:
            handle.write(HEADER)
        return 0

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(HEADER)
        for key in keys:
            tract, n = key
            draws = total * max(args.permutations, 1)
            obs = 100.0 * observed[key] / total
            exp = 100.0 * permuted[key] / draws
            # Half-count smoothing on the permuted side. A long run can be
            # absent from every shuffle, and an unsmoothed ratio would then be
            # infinite -- which has to be rendered as *something* in a plot, and
            # every choice (0, a cap) misreads as its opposite. The smoothed
            # ratio stays finite and monotonic in the observed count.
            exp_smoothed = 100.0 * (permuted[key] + 0.5) / (draws + 0.5)
            enrichment = obs / exp_smoothed
            handle.write(f"{args.sample}\t{tract}\t{n}\t{total}\t{obs:.4f}\t{exp:.4f}\t{enrichment:.3f}\n")

    mean_gc = 100.0 * gc_total / total
    draws = total * max(args.permutations, 1)
    summary = ", ".join(
        f"{tract} {n}nt obs={100.0 * observed[key] / total:.3f}% "
        f"exp={100.0 * permuted[key] / draws:.3f}%"
        for key in keys
        for tract, n in (key,)
    )
    print(f"[polya] {args.sample}: {total} reads, mean GC {mean_gc:.1f}%; {summary}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

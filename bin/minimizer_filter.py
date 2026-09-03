#!/usr/bin/env python3
"""
minimizer_filter.py -- is a taxon's read count backed by BREADTH of evidence?

An abundance filter asks how many reads landed on a taxon. That is exactly the
question a spurious call passes: reads piled on one conserved gene, one
repetitive element, one adapter-like stretch, look abundant. What separates them
from a real organism is how much of the reference those reads actually cover.

Both classifiers this pipeline can run report that, in slightly different terms:

  * Kraken2 with `--report-minimizer-data` inserts two columns - the minimizers
    observed for a clade and, of those, how many were DISTINCT. Minimizers are a
    subsample of k-mers (roughly one per five positions at the default k=35,
    l=31), so the counts are coarser than KrakenUniq's but behave the same way.
  * KrakenUniq reports `kmers` (distinct k-mers, by HyperLogLog), `dup` (their
    average duplication) and `cov` (the fraction of the reference they span)
    directly, with no extra flag.

Three quantities are derived, whichever the source:

  duplication  observed / distinct. A real organism spreads reads over the
               genome, so most of what it hits is new; a taxon called from one
               locus re-hits the same handful over and over.
  distinct     the raw count. Below a few, the call rests on one place.
  coverage     distinct observed / distinct held in the database for that clade.
               KrakenUniq gives it; for Kraken2 the denominator comes from the
               `inspect.txt` shipped with the database, whose second column is
               the distinct minimizers per clade.

Reads alone are never a reason to drop: a taxon under --min-reads is left
untouched, because a handful of reads cannot demonstrate breadth either way and
condemning them would just re-derive the abundance filter that already ran.

--correlation-filter adds a second, independent question, taken from SAHMI's
sample-level denoising step (Ghaddar, Blaser & De, Nat Comput Sci 2023). The
thresholds above ask whether a taxon's evidence is broad ENOUGH; the correlation
asks whether it SCALES. Across samples, a real organism that is more abundant in
one library contributes proportionally more k-mers and more distinct k-mers
there too, so all three of

    reads ~ observed,    reads ~ distinct,    observed ~ distinct

rise together. A reagent contaminant sitting at a flat low level in every
library, or reads piling on one conserved locus, breaks that proportionality:
the read count moves and the distinct count does not. SAHMI requires all three
Spearman correlations to be significant, at genus or species rank, for a taxon
seen in at least three samples.

The two tests catch different things and neither subsumes the other. A constant
contaminant with plenty of distinct minimizers passes the thresholds and fails
the correlation; a genuine bloom in a single library passes the correlation
(where it can be computed at all) and is judged only on its thresholds. They are
combined by union - a taxon is dropped if either condemns it.

One hard limit, enforced rather than documented: with four or fewer samples the
smallest attainable two-sided Spearman p-value is 2/4! = 0.083, so no taxon can
clear p < 0.05 however clean its evidence. Five samples is the minimum at which
this test can do anything; below that it is refused up front.
"""
import argparse
import csv
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from sahmi_stats import MIN_N_FOR_SIGNIFICANCE, benjamini_hochberg, spearman

RANK_CODE = re.compile(r"^[URDKPCOFGS][0-9]*$")
# KrakenUniq spells ranks out; these are the ones its reports use.
RANK_WORDS = {
    "no rank", "superkingdom", "kingdom", "phylum", "class", "order",
    "family", "genus", "species", "subspecies", "clade", "domain",
    "sequence", "assembly", "root",
}


def parse_report(path):
    """Return {taxid: {...}} for one Kraken2-minimizer or KrakenUniq report.

    Layouts, by position:
      Kraken2 + --report-minimizer-data
        pct, clade_reads, direct_reads, minimizers, distinct, rank, taxid, name
      KrakenUniq
        pct, reads, taxReads, kmers, dup, cov, taxid, rank, name

    They are told apart by where the taxid sits - Kraken2 puts the rank before
    it, KrakenUniq after - which is the same test Bracken and KrakenTools use.
    """
    rows = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("%"):
                continue
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if len(fields) < 8:
                continue
            try:
                clade_reads = int(fields[1])
            except ValueError:
                continue

            name = fields[-1].strip()
            krakenuniq = fields[-2].lower() in RANK_WORDS
            try:
                taxid = int(fields[-3]) if krakenuniq else int(fields[-2])
            except ValueError:
                continue
            rank = fields[-2] if krakenuniq else fields[-3]

            if krakenuniq:
                # kmers, dup, cov are columns 4, 5, 6.
                distinct = to_number(fields[3])
                duplication = to_number(fields[4], default=None)
                coverage = to_number(fields[5], default=None)
                observed = distinct * duplication if duplication else None
            else:
                observed = to_number(fields[3])
                distinct = to_number(fields[4])
                duplication = observed / distinct if distinct else None
                coverage = None

            # taxid 0 is the unclassified bucket and 1 is the root: bookkeeping
            # rows, not organisms, and neither has a reference to cover.
            if taxid in (0, 1):
                continue

            rows[taxid] = {
                "taxid": taxid,
                "rank": rank,
                "name": name,
                "reads": clade_reads,
                "observed": observed,
                "distinct": distinct,
                "duplication": duplication,
                "coverage": coverage,
            }
    return rows


def to_number(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def parse_inspect(path):
    """{taxid: distinct minimizers in the database for that clade}.

    `kraken2-inspect` is, in its own words, a wrapper that reports minimizer
    counts per taxon; its second column is the clade total. That is the
    denominator Kraken2's report does not carry.
    """
    totals = {}
    with open(path, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = [f.strip() for f in line.rstrip("\n").split("\t")]
            if len(fields) < 5:
                continue
            try:
                totals[int(fields[-2])] = float(fields[1])
            except ValueError:
                continue
    return totals


# The ranks SAHMI restricts its sample-level test to, in both spellings: Kraken2
# writes rank codes, KrakenUniq writes the words. Sub-ranks (G1, S2) count as
# their parent, which is how Kraken2 spells strain-level clades.
CORRELATION_RANKS = {"g": "genus", "s": "species"}


def rank_key(rank):
    """'S1' -> 'species', 'genus' -> 'genus', anything else -> None.

    Kraken2 writes rank codes and spells strain-level clades as sub-ranks (S1,
    S2); KrakenUniq spells the words out and calls the same thing 'subspecies'.
    Both are counted as their parent, so the two classifiers' reports are judged
    by the same rule rather than S1 qualifying where 'subspecies' does not.
    """
    text = (rank or "").strip().lower()
    if text in ("genus", "species"):
        return text
    if text == "subspecies":
        return "species"
    if text[:1] in CORRELATION_RANKS and (text[1:] == "" or text[1:].isdigit()):
        return CORRELATION_RANKS[text[:1]]
    return None


def correlation_test(per_sample, ranks, args):
    """{taxid: {...}} for taxa the three-correlation test could judge.

    `per_sample` is {taxid: [(reads, observed, distinct), ...]} with one entry
    per sample the taxon appeared in. Taxa the test cannot speak about - too few
    samples, too little evidence, wrong rank, or a vector that never varies -
    are simply absent from the result and keep whatever the thresholds said.
    """
    raw = {}
    for taxid, observations in per_sample.items():
        # SAHMI applies this at genus and species resolution only. Higher ranks
        # are clade sums whose counts move for reasons that have nothing to do
        # with any one organism, so proportionality there means little.
        if rank_key(ranks.get(taxid)) is None:
            continue
        if len(observations) < args.correlation_min_samples:
            continue
        usable = [
            (reads, observed, distinct)
            for reads, observed, distinct in observations
            if reads > args.correlation_min_reads
            and (distinct or 0) > args.correlation_min_distinct
        ]
        if len(usable) < args.correlation_min_samples:
            continue
        reads = [entry[0] for entry in usable]
        observed = [entry[1] or 0.0 for entry in usable]
        distinct = [entry[2] or 0.0 for entry in usable]
        tests = {
            "reads_observed": spearman(reads, observed),
            "reads_distinct": spearman(reads, distinct),
            "observed_distinct": spearman(observed, distinct),
        }
        # A pair whose statistic is undefined - one vector constant across every
        # sample - has not been tested, so the taxon cannot be said to have
        # failed. It is dropped from consideration rather than condemned.
        if any(rho is None for rho, _p in tests.values()):
            continue
        raw[taxid] = {"n": len(usable), "tests": tests}

    if not raw:
        return {}

    # BH within each of the three tests, over the taxa that were actually
    # tested: correcting across the three as one family would treat a taxon's
    # own three p-values as independent hypotheses, which they are not.
    taxids = sorted(raw)
    adjusted = {}
    for name in ("reads_observed", "reads_distinct", "observed_distinct"):
        qvalues = benjamini_hochberg([raw[taxid]["tests"][name][1] for taxid in taxids])
        for taxid, qvalue in zip(taxids, qvalues):
            adjusted.setdefault(taxid, {})[name] = qvalue

    results = {}
    for taxid in taxids:
        tests = raw[taxid]["tests"]
        qvalues = adjusted[taxid]
        results[taxid] = {
            "n": raw[taxid]["n"],
            "rho_reads_observed": tests["reads_observed"][0],
            "rho_reads_distinct": tests["reads_distinct"][0],
            "rho_observed_distinct": tests["observed_distinct"][0],
            "q_max": max(qvalues.values()),
            # SAHMI's rule: every one of the three has to hold.
            "passed": all(q <= args.correlation_p for q in qvalues.values()),
        }
    return results


def verdict(row, args):
    """Why a taxon is being dropped, or None if it survives."""
    if row["reads"] < args.min_reads:
        return None  # too few reads to judge breadth either way
    reasons = []
    # Union with the correlation test: either question can condemn a taxon, and
    # a taxon the correlation could not judge carries `correlation` as None.
    if row.get("correlation") is not None and not row["correlation"]["passed"]:
        reasons.append("no_correlation")
    if row["distinct"] is not None and row["distinct"] < args.min_distinct:
        reasons.append("few_distinct")
    if row["duplication"] is not None and row["duplication"] > args.max_duplication:
        reasons.append("high_duplication")
    if (
        args.min_coverage > 0
        and row["coverage"] is not None
        and row["coverage"] < args.min_coverage
    ):
        reasons.append("low_coverage")
    return "+".join(reasons) if reasons else None


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reports", nargs="+", required=True,
                        help="per-sample Kraken2 (--report-minimizer-data) or KrakenUniq reports")
    parser.add_argument("--inspect", help="kraken2-inspect output, for the coverage denominator")
    parser.add_argument("--evidence", required=True, help="per-taxon evidence table to write")
    parser.add_argument("--drop-list", required=True, help="taxids to drop, one per line")
    parser.add_argument("--mqc", help="MultiQC custom-content section to write")
    parser.add_argument("--min-reads", type=float, default=50,
                        help="only judge taxa with at least this many reads")
    parser.add_argument("--min-distinct", type=float, default=10,
                        help="a taxon needs at least this many distinct minimizers/k-mers")
    parser.add_argument("--distinct-scale", type=float, default=0.0,
                        help="scale the distinct threshold with library depth: the effective "
                             "threshold becomes max(--min-distinct, scale x classified reads). "
                             "0 disables. See the module docstring for the benchmarked value.")
    parser.add_argument("--max-duplication", type=float, default=100,
                        help="...and its observed/distinct ratio must stay below this")
    parser.add_argument("--min-coverage", type=float, default=0.0,
                        help="...and cover at least this fraction of its reference (0 disables)")
    parser.add_argument("--correlation-filter", action="store_true",
                        help="also require a taxon's reads, k-mers and distinct k-mers to rise "
                             "together across samples (SAHMI's sample-level denoising)")
    parser.add_argument("--correlation-p", type=float, default=0.05,
                        help="BH-adjusted p every one of the three correlations must clear")
    parser.add_argument("--correlation-min-samples", type=int, default=3,
                        help="a taxon must appear in at least this many samples to be tested")
    parser.add_argument("--correlation-min-reads", type=float, default=2,
                        help="a sample counts towards the test only above this many reads")
    parser.add_argument("--correlation-min-distinct", type=float, default=2,
                        help="...and above this many distinct minimizers/k-mers")
    args = parser.parse_args()

    if args.correlation_filter and len(args.reports) < MIN_N_FOR_SIGNIFICANCE:
        raise SystemExit(
            f"minimizer_filter: --correlation-filter needs at least {MIN_N_FOR_SIGNIFICANCE} samples, "
            f"but {len(args.reports)} were given. With n samples the smallest attainable two-sided "
            "Spearman p-value is 2/n!, so at n=4 nothing can reach 0.083 let alone --correlation-p; "
            "the filter would drop every taxon it tested. Run without it on cohorts this small."
        )

    db_totals = parse_inspect(args.inspect) if args.inspect else {}

    # Depth scaling. Sun et al.'s threshold sweep found the optimal unique-k-mer
    # cut is linearly associated with sequencing depth (p = 7e-7), at roughly
    # 0.002 unique k-mers per read - 200 per 100,000 reads - which reproduces
    # and refines KrakenUniq's original "about 2000 unique k-mers per million
    # reads". A flat threshold across libraries that differ 100-fold in depth,
    # which is normal in public-data reanalysis, is therefore wrong in both
    # directions at once: too strict for shallow libraries, far too lenient for
    # deep ones. Depth is taken as the classified reads in that report, so it
    # needs no external input.
    scaled_floor = 0.0
    if args.distinct_scale > 0:
        depths = []
        for path in args.reports:
            rows = parse_report(path)
            depths.append(sum(row["reads"] for row in rows.values() if row["rank"].upper().startswith("S")))
        if depths:
            scaled_floor = args.distinct_scale * (sum(depths) / len(depths))
            args.min_distinct = max(args.min_distinct, scaled_floor)
            print(
                f"[minimizer_filter] depth scaling: mean {sum(depths) / len(depths):.0f} classified "
                f"reads x {args.distinct_scale:g} -> distinct threshold {args.min_distinct:.0f}.",
                file=sys.stderr,
            )

    # A taxon is judged at its BEST showing across samples: the most distinct
    # minimizers and the lowest duplication it managed anywhere. Breadth in one
    # library is enough to show the organism is real; demanding it in every
    # library would just penalise the shallow ones.
    best = {}
    # Per-taxon vectors across samples, for the correlation test. Kept beside
    # `best` rather than derived from it: `best` is a collapse, and a collapse
    # is exactly what a correlation cannot be computed from.
    per_sample = {}
    for path in args.reports:
        sample = os.path.basename(path)
        for taxid, row in parse_report(path).items():
            per_sample.setdefault(taxid, []).append(
                (row["reads"], row["observed"], row["distinct"])
            )
            if db_totals and row["coverage"] is None and row["distinct"] is not None:
                total = db_totals.get(taxid)
                if total:
                    row["coverage"] = row["distinct"] / total
            entry = best.get(taxid)
            if entry is None:
                best[taxid] = dict(row, samples=1, best_in=sample)
                continue
            entry["samples"] += 1
            entry["reads"] = max(entry["reads"], row["reads"])
            if (row["distinct"] or 0) > (entry["distinct"] or 0):
                entry["distinct"] = row["distinct"]
                entry["observed"] = row["observed"]
                entry["best_in"] = sample
            # Best = most breadth: the lowest duplication and highest coverage
            # the taxon reached in any single library.
            if row["duplication"] is not None:
                entry["duplication"] = (
                    row["duplication"] if entry["duplication"] is None
                    else min(entry["duplication"], row["duplication"])
                )
            if row["coverage"] is not None:
                entry["coverage"] = (
                    row["coverage"] if entry["coverage"] is None
                    else max(entry["coverage"], row["coverage"])
                )

    if not best:
        raise SystemExit(
            "minimizer_filter: no minimizer columns found in any report. For Kraken2 this "
            "needs --kraken2_report_minimizer_data; KrakenUniq reports carry them already."
        )

    correlations = {}
    if args.correlation_filter:
        correlations = correlation_test(
            per_sample, {taxid: row["rank"] for taxid, row in best.items()}, args
        )
        print(
            f"[minimizer_filter] correlation test judged {len(correlations)} of {len(best)} taxa "
            f"across {len(args.reports)} samples "
            f"({sum(1 for r in correlations.values() if r['passed'])} passed).",
            file=sys.stderr,
        )

    dropped, counts = [], {}
    with open(args.evidence, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "taxid", "rank", "name", "reads", "samples",
            "distinct", "observed", "duplication", "coverage",
            "cor_samples", "rho_reads_observed", "rho_reads_distinct",
            "rho_observed_distinct", "cor_q_max", "verdict",
        ])
        for taxid in sorted(best, key=lambda t: -best[t]["reads"]):
            row = best[taxid]
            row["correlation"] = correlations.get(taxid)
            reason = verdict(row, args)
            call = reason or ("kept" if row["reads"] >= args.min_reads else "below_read_gate")
            counts[call] = counts.get(call, 0) + 1
            cor = row["correlation"]
            writer.writerow([
                taxid, row["rank"], row["name"], int(row["reads"]), row["samples"],
                fmt(row["distinct"], 0), fmt(row["observed"], 0),
                fmt(row["duplication"], 3), fmt(row["coverage"], 6),
                cor["n"] if cor else "NA",
                fmt(cor["rho_reads_observed"], 3) if cor else "NA",
                fmt(cor["rho_reads_distinct"], 3) if cor else "NA",
                fmt(cor["rho_observed_distinct"], 3) if cor else "NA",
                fmt(cor["q_max"], 6) if cor else "NA",
                call,
            ])
            if reason:
                dropped.append(taxid)

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in dropped:
            handle.write(f"{taxid}\n")

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_minimizer_filter'\n"
                "# section_name: 'Minimizer evidence'\n"
                "# description: 'Taxa judged on the BREADTH of their evidence rather than its\n"
                "#     volume: how many distinct minimizers (Kraken2) or k-mers (KrakenUniq) their\n"
                "#     reads cover, and how often the same ones are re-hit. Reads concentrated on one\n"
                "#     locus look abundant but cover almost nothing. Taxa under the read gate are not\n"
                "#     judged - too few reads to show breadth either way.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_minimizer_filter_plot'\n"
                "#     title: 'reanaTax: minimizer evidence'\n"
                "#     ylab: 'Taxa'\n"
            )
            keys = sorted(counts)
            handle.write("Sample\t" + "\t".join(keys) + "\n")
            handle.write("all taxa\t" + "\t".join(str(counts[key]) for key in keys) + "\n")

    print(
        f"[minimizer_filter] {len(best)} taxa: "
        + ", ".join(f"{count} {call}" for call, count in sorted(counts.items()))
        + f" (>= {args.min_distinct:g} distinct, duplication <= {args.max_duplication:g}"
        + (f", coverage >= {args.min_coverage:g}" if args.min_coverage > 0 else "")
        + f", gate {args.min_reads:g} reads); {len(dropped)} taxid(s) listed for removal.",
        file=sys.stderr,
    )
    return 0


def fmt(value, places=1):
    if value is None:
        return "NA"
    return f"{value:.{places}f}" if places else str(int(value))


if __name__ == "__main__":
    sys.exit(main())

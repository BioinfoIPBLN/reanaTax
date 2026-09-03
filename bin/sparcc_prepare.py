#!/usr/bin/env python3
"""
sparcc_prepare.py -- the combined Bracken table as an OTU table FastSpar reads.

Also, and more importantly, the place where the question "can this cohort
support a co-occurrence network at all?" is asked out loud.

Why SparCC and not a correlation matrix.  Taxonomic profiles are compositional:
the counts sum to a library size that has nothing to do with the biology, so
when one taxon rises every other one falls whether or not anything about them is
related. A Pearson or Spearman matrix over such data manufactures negative
correlations everywhere and a spurious positive block among the rare taxa.
SparCC (Friedman & Alm, PLoS Comput Biol 2012) works on log-ratio variances,
which are invariant to the total, and infers the correlations of the underlying
basis instead. That is the same reason ANCOM-BC2 and ALDEx2 are the pipeline's
differential-abundance methods rather than a t-test on proportions.

Why it needs a guard.  SparCC assumes the true network is sparse and that most
pairs are uncorrelated; it estimates each pair from the variance of a log ratio
across samples. With a handful of samples those variances are noise, and the
bootstrap p-values that come out will be small for some pair purely because a
bootstrap resample of five samples has very few distinct outcomes. The refusal
below is not conservatism, it is the difference between "no associations were
found" and "this cohort cannot be asked the question".

Zeros are the other reason to be careful. A log ratio is undefined at zero, so
SparCC adds a pseudocount, and a taxon that is zero in most samples has its
correlations decided almost entirely by that pseudocount. Taxa are therefore
filtered on prevalence here rather than left for the tool to cope with.
"""
import argparse
import csv
import gzip
import math
import sys

# Below this, the log-ratio variances SparCC is built on are estimated from too
# few observations for the result to mean anything, whatever p-values come out.
MIN_SAMPLES = 10


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(
        path, encoding="utf-8", errors="replace"
    )


def read_bracken(path):
    """(taxa, samples, counts) from the combined Bracken table.

    `taxa` is a list of (taxid, name); `counts[i][j]` is taxon i in sample j.
    The `_num` columns are used and the `_frac` ones ignored: SparCC is given
    counts, and a table already closed to 1 has had the very quantity it needs
    removed.
    """
    with open_maybe_gzip(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None:
            sys.exit(f"[sparcc_prepare] '{path}' is empty.")
        header = [field.strip() for field in header]
        columns = [index for index, field in enumerate(header) if field.endswith("_num")]
        if not columns:
            sys.exit(
                f"[sparcc_prepare] no `<sample>_num` column in '{path}'. This needs the combined "
                "Bracken table (bracken_combined_<level>.txt), not a relative-abundance profile: "
                "SparCC is given counts, and a table already closed to 1 has had the quantity it "
                "works on removed."
            )
        samples = [header[index][: -len("_num")] for index in columns]
        try:
            name_column = header.index("name")
        except ValueError:
            sys.exit(f"[sparcc_prepare] no `name` column in '{path}'.")
        taxid_column = header.index("taxonomy_id") if "taxonomy_id" in header else None

        taxa, counts = [], []
        for row in reader:
            if not row or len(row) <= max(columns):
                continue
            values = []
            for index in columns:
                try:
                    values.append(int(round(float(row[index]))))
                except ValueError:
                    values.append(0)
            taxa.append((row[taxid_column].strip() if taxid_column is not None else "",
                         row[name_column].strip()))
            counts.append(values)
    return taxa, samples, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True, help="combined Bracken table")
    parser.add_argument("--min-prevalence", type=float, default=0.5,
                        help="a taxon must be seen in this fraction of samples")
    parser.add_argument("--min-reads", type=int, default=10, help="...with at least this many reads there")
    parser.add_argument("--top-taxa", type=int, default=200,
                        help="keep only the N most abundant survivors (0 = all)")
    parser.add_argument("--permutations", type=int, default=1000,
                        help="bootstrap replicates the p-values will be built from")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--force", action="store_true", help="run even when the cohort is too small")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    taxa, samples, counts = read_bracken(args.counts)
    if len(samples) < MIN_SAMPLES and not args.force:
        sys.exit(
            f"[sparcc_prepare] REFUSING: {len(samples)} sample(s), and SparCC needs at least "
            f"{MIN_SAMPLES} before its log-ratio variances are estimated from enough observations "
            "to mean anything. Every pair's correlation would be decided by two or three points, "
            "and the bootstrap p-values would be small for some pairs purely because a resample of "
            f"{len(samples)} samples has very few distinct outcomes. An empty or a full network "
            "here would both be artefacts of the cohort size, not findings. Pass --force to run "
            "anyway and read the result as exploratory."
        )

    # Prevalence first: a taxon that is zero in most samples has its
    # correlations decided by SparCC's pseudocount rather than by its data.
    needed = math.ceil(args.min_prevalence * len(samples))
    keep = []
    for index, values in enumerate(counts):
        present = sum(1 for value in values if value >= args.min_reads)
        if present >= needed:
            keep.append(index)
    if args.top_taxa and len(keep) > args.top_taxa:
        keep.sort(key=lambda index: -sum(counts[index]))
        keep = keep[: args.top_taxa]
    keep.sort()

    if len(keep) < 3:
        sys.exit(
            f"[sparcc_prepare] only {len(keep)} taxon/taxa are present in >= {args.min_prevalence:g} "
            f"of samples with >= {args.min_reads} reads. A co-occurrence network needs at least "
            "three nodes to say anything a pairwise test could not."
        )

    pairs = len(keep) * (len(keep) - 1) // 2
    floor = 1.0 / args.permutations
    # BH can only reject at all when the smallest attainable p clears alpha/m.
    if floor > args.p_threshold / pairs and not args.force:
        needed_permutations = int(math.ceil(pairs / args.p_threshold))
        sys.exit(
            f"[sparcc_prepare] REFUSING: {pairs} pairs will be tested, so a Benjamini-Hochberg "
            f"threshold of {args.p_threshold:g} needs some pair to reach p <= "
            f"{args.p_threshold / pairs:.3g}. With {args.permutations} bootstrap replicates the "
            f"smallest p FastSpar can report is {floor:.3g}, so NO pair could be called "
            "significant however strong its correlation. Raise --sparcc_permutations to at least "
            f"{needed_permutations}, cut --sparcc_top_taxa, or pass --force and read the "
            "correlations without the p-values."
        )

    with open(f"{args.prefix}.sparcc_otu.tsv", "w", encoding="utf-8") as handle:
        handle.write("#OTU ID\t" + "\t".join(samples) + "\n")
        for index in keep:
            taxid, name = taxa[index]
            # FastSpar keys rows by the id string, so it has to be unique and
            # free of whitespace; the readable name is restored by the network
            # step from the map below.
            label = taxid if taxid else name.replace(" ", "_")
            handle.write(label + "\t" + "\t".join(str(value) for value in counts[index]) + "\n")

    with open(f"{args.prefix}.sparcc_taxa.tsv", "w", encoding="utf-8") as handle:
        handle.write("id\ttaxid\tname\ttotal_reads\tprevalence\n")
        for index in keep:
            taxid, name = taxa[index]
            label = taxid if taxid else name.replace(" ", "_")
            present = sum(1 for value in counts[index] if value >= args.min_reads)
            handle.write(f"{label}\t{taxid}\t{name}\t{sum(counts[index])}\t"
                         f"{present / len(samples):.4f}\n")

    print(
        f"[sparcc_prepare] {len(keep)}/{len(taxa)} taxa kept across {len(samples)} samples; "
        f"{pairs} pairs, smallest attainable p {floor:.3g} against a BH requirement of "
        f"{args.p_threshold / pairs:.3g}.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

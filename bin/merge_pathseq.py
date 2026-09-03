#!/usr/bin/env python3
"""
merge_pathseq.py -- PathSeq's per-sample scores into one table per cohort.

PathSeq writes one scores file per sample with a row per taxon and these
columns: tax_id, taxonomy, type, name, kingdom, score, score_normalized, reads,
unambiguous, reference_length. Three of them are worth understanding before
reading any output built from them.

  reads         the number of reads assigned to the taxon, INCLUDING those
                shared with other taxa. Not an integer count of distinct reads
                at that node in the Kraken sense - a read aligning equally well
                to three species contributes to all three.
  unambiguous   the reads that aligned to this taxon and nowhere else. This is
                the conservative count, and it is the one to compare against a
                Kraken2 species count.
  score         PathSeq's abundance score: each read's weight divided among the
                taxa it maps to, then summed up the tree and divided by the
                reference length. It is a length-normalised abundance, so it is
                comparable BETWEEN taxa within a sample in a way `reads` is not,
                and it is not a count.

Three matrices come out, one per quantity, taxa by sample - the same shape the
combined Bracken table has, so the same tools can read them. `score_normalized`
is PathSeq's own within-sample percentage and is emitted as given rather than
recomputed.

Deliberately NOT converted into a Kraken report. The two tools disagree about
what a read assigned to a taxon means - PathSeq divides an ambiguous read
between taxa, Kraken2 pushes it up to their common ancestor - and a conversion
would have to invent one convention or the other. Keeping PathSeq's own numbers
means a disagreement between the two routes stays visible instead of being
averaged away by the format.
"""
import argparse
import csv
import gzip
import os
import sys

QUANTITIES = ("score", "score_normalized", "reads", "unambiguous")


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(
        path, encoding="utf-8", errors="replace"
    )


def sample_name(path):
    base = os.path.basename(path)
    for suffix in (".pathseq.scores.txt", ".scores.txt", ".txt", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def read_scores(path):
    """{taxid: {quantity: value}} plus {taxid: (name, type, kingdom, taxonomy)}."""
    values, meta = {}, {}
    with open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "tax_id" not in [f.strip() for f in reader.fieldnames]:
            sys.exit(
                f"[merge_pathseq] '{path}' does not look like a PathSeq scores file: its header is "
                f"{reader.fieldnames}. Expected a tax_id column."
            )
        lookup = {name.strip(): name for name in reader.fieldnames}
        for row in reader:
            taxid = (row.get(lookup.get("tax_id", "")) or "").strip()
            if not taxid:
                continue
            entry = {}
            for quantity in QUANTITIES:
                raw = (row.get(lookup.get(quantity, ""), "") or "").strip()
                try:
                    entry[quantity] = float(raw)
                except ValueError:
                    entry[quantity] = 0.0
            values[taxid] = entry
            meta[taxid] = (
                (row.get(lookup.get("name", ""), "") or "").strip(),
                (row.get(lookup.get("type", ""), "") or "").strip(),
                (row.get(lookup.get("kingdom", ""), "") or "").strip(),
                (row.get(lookup.get("taxonomy", ""), "") or "").strip(),
            )
    return values, meta


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scores", nargs="+", required=True, help="PathSeq scores files, one per sample")
    parser.add_argument("--rank", default=None,
                        help="keep only this PathSeq `type` (e.g. species, genus); default keeps every row")
    parser.add_argument("--prefix", default="pathseq")
    args = parser.parse_args()

    samples = []
    per_sample = {}
    meta = {}
    for path in sorted(args.scores):
        sample = sample_name(path)
        if sample in per_sample:
            sys.exit(f"[merge_pathseq] two scores files resolve to the sample id '{sample}'.")
        values, sample_meta = read_scores(path)
        samples.append(sample)
        per_sample[sample] = values
        for taxid, entry in sample_meta.items():
            meta.setdefault(taxid, entry)

    taxids = sorted(meta, key=lambda taxid: -sum(
        per_sample[sample].get(taxid, {}).get("reads", 0.0) for sample in samples
    ))
    if args.rank:
        wanted = args.rank.strip().lower()
        taxids = [taxid for taxid in taxids if meta[taxid][1].lower() == wanted]
        if not taxids:
            sys.exit(
                f"[merge_pathseq] no row has type '{args.rank}'. PathSeq's `type` column holds the "
                "rank name spelled out (species, genus, family, ...), not a Kraken rank code."
            )

    for quantity in QUANTITIES:
        with open(f"{args.prefix}_{quantity}.tsv", "w", encoding="utf-8") as handle:
            handle.write("\t".join(["taxonomy_id", "name", "type", "kingdom"] + samples) + "\n")
            for taxid in taxids:
                name, rank, kingdom, _lineage = meta[taxid]
                row = [taxid, name, rank, kingdom]
                for sample in samples:
                    value = per_sample[sample].get(taxid, {}).get(quantity, 0.0)
                    row.append(f"{value:.6g}")
                handle.write("\t".join(row) + "\n")

    with open(f"{args.prefix}_lineage.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxonomy_id\tname\ttype\tkingdom\ttaxonomy\n")
        for taxid in taxids:
            name, rank, kingdom, lineage = meta[taxid]
            handle.write(f"{taxid}\t{name}\t{rank}\t{kingdom}\t{lineage}\n")

    totals = {
        sample: sum(entry.get("unambiguous", 0.0) for entry in per_sample[sample].values())
        for sample in samples
    }
    with open(f"{args.prefix}_pathseq_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            "# id: 'reanatax_pathseq'",
            "# section_name: 'PathSeq'",
            "# description: 'Reads assigned by GATK PathSeq, which ALIGNS the surviving reads to a",
            "#     microbe reference rather than matching k-mers, and divides an ambiguously",
            "#     mapping read between the taxa it hits instead of pushing it to their common",
            "#     ancestor. Only unambiguous reads are plotted - the reads that aligned to one",
            "#     taxon and nowhere else - because those are the ones comparable to a Kraken2",
            "#     species count.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            "#     id: 'reanatax_pathseq_plot'",
            "#     title: 'reanaTax: PathSeq unambiguous reads'",
            "#     ylab: 'Reads'",
            "Sample\tUnambiguous reads",
        ] + [f"{sample}\t{totals[sample]:.0f}" for sample in samples] + [""]))

    print(
        f"[merge_pathseq] {len(taxids)} taxa across {len(samples)} sample(s)"
        + (f" at type '{args.rank}'" if args.rank else "")
        + f"; {len(QUANTITIES)} matrices written.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

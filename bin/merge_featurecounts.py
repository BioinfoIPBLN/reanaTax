#!/usr/bin/env python3
"""Merge per-sample featureCounts tables into one gene x sample count matrix.

featureCounts names its count column after the BAM it was given
(`SRX8592588.host.sorted.bam`), so the column has to be renamed to the sample id
before the samples can sit side by side and be matched against --da_metadata.

Gene lengths are carried through in a second file: they are constant across
samples, and both the length-normalised expression used for the host-microbe
correlation and any downstream TPM need them.
"""

import argparse
import os
import sys


def sample_name(path):
    """`SRX8592588.host.featureCounts.tsv` -> `SRX8592588`. The `.host` label is
    the depletion pass the BAM came from (see conf/modules.config), not part of
    the sample id, so it is stripped along with the extension."""
    name = os.path.basename(path)
    for suffix in (".featureCounts.tsv", ".featureCounts.txt", ".tsv", ".txt"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    for label in (".host", ".host1", ".host2", ".host3"):
        if name.endswith(label):
            name = name[: -len(label)]
            break
    return name


def read_table(path):
    """featureCounts writes a `# Program:...` line, then a header, then one row
    per gene with six annotation columns and one count column."""
    counts = {}
    lengths = {}
    order = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Geneid":
                continue
            if len(fields) < 7:
                continue
            gene = fields[0]
            try:
                lengths[gene] = int(fields[5])
                # The count is the LAST column, not column 7: featureCounts adds
                # one column per BAM, and while this pipeline passes exactly one,
                # keying on the last column keeps the reader honest if that ever
                # changes.
                counts[gene] = int(float(fields[-1]))
            except ValueError:
                continue
            order.append(gene)
    return order, counts, lengths


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tables", nargs="+", help="per-sample featureCounts tables")
    parser.add_argument("--output", required=True, help="merged count matrix (TSV)")
    parser.add_argument("--lengths", help="gene lengths (TSV), optional")
    args = parser.parse_args()

    samples = []
    per_sample = {}
    gene_order = []
    seen = set()
    lengths = {}

    for path in sorted(args.tables):
        name = sample_name(path)
        if name in per_sample:
            sys.exit(f"[merge_featurecounts] two tables map to sample '{name}'; "
                     "the count matrix would silently keep only one of them")
        order, counts, lens = read_table(path)
        if not counts:
            sys.exit(f"[merge_featurecounts] {path} holds no gene rows")
        samples.append(name)
        per_sample[name] = counts
        for gene in order:
            if gene not in seen:
                seen.add(gene)
                gene_order.append(gene)
        # Lengths are a property of the annotation, so every sample should agree.
        # A disagreement means the tables were built against different GTFs and
        # the matrix would be quietly meaningless.
        for gene, value in lens.items():
            if gene in lengths and lengths[gene] != value:
                sys.exit(f"[merge_featurecounts] gene '{gene}' has length {lengths[gene]} "
                         f"in one table and {value} in {path}; the tables were not "
                         "counted against the same annotation")
            lengths[gene] = value

    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write("gene_id\t" + "\t".join(samples) + "\n")
        for gene in gene_order:
            row = [str(per_sample[s].get(gene, 0)) for s in samples]
            handle.write(gene + "\t" + "\t".join(row) + "\n")

    if args.lengths:
        with open(args.lengths, "w", encoding="utf-8") as handle:
            handle.write("gene_id\tlength\n")
            for gene in gene_order:
                handle.write(f"{gene}\t{lengths.get(gene, 0)}\n")

    print(f"[merge_featurecounts] {len(gene_order)} genes x {len(samples)} samples "
          f"-> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

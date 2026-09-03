#!/usr/bin/env python3
"""Merge per-sample kallisto abundances into one counts matrix.

The output is deliberately the same shape bin/merge_featurecounts.py produces -
a `gene_id` column followed by one column per sample - so the same differential
expression path serves both quantifiers and there is only one of it to keep
correct.

Two things are worth knowing about the numbers.

`est_counts` are ESTIMATED and not integers: a read that pseudoaligns to several
transcripts of the same gene is apportioned between them. DESeq2 is handed
rounded counts downstream, which is an approximation; the fully correct route is
tximport, which passes kallisto's effective lengths to the model as offsets
rather than discarding them. Summing transcripts to genes with --tx2gene removes
most of the difference, because the apportioning that rounding damages is mostly
*within* a gene.

Effective length, not length, is what kallisto normalises by, so that is what is
carried in the lengths file.
"""

import argparse
import os
import sys


def sample_name(path):
    name = os.path.basename(path)
    for suffix in (".abundance.tsv", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def read_abundance(path):
    counts, tpms, lengths = {}, {}, {}
    with open(path, encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            i_id = header.index("target_id")
            i_eff = header.index("eff_length")
            i_cnt = header.index("est_counts")
            i_tpm = header.index("tpm")
        except ValueError:
            sys.exit(f"[merge_kallisto] {path} is not a kallisto abundance.tsv "
                     f"(header: {', '.join(header)})")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(i_id, i_eff, i_cnt, i_tpm):
                continue
            try:
                counts[fields[i_id]] = float(fields[i_cnt])
                tpms[fields[i_id]] = float(fields[i_tpm])
                lengths[fields[i_id]] = float(fields[i_eff])
            except ValueError:
                continue
    return counts, tpms, lengths


def read_tx2gene(path):
    mapping = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2 and fields[0] and fields[1]:
                if fields[0].lower() in ("transcript_id", "tx", "target_id"):
                    continue  # header row
                mapping[fields[0]] = fields[1]
    if not mapping:
        sys.exit(f"[merge_kallisto] {path} yielded no transcript->gene pairs; it "
                 "should be a two-column TSV of transcript id and gene id")
    return mapping


def write_matrix(path, features, samples, table, fmt="{:.4f}"):
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("gene_id\t" + "\t".join(samples) + "\n")
        for feature in features:
            row = [fmt.format(table[s].get(feature, 0.0)) for s in samples]
            handle.write(feature + "\t" + "\t".join(row) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("abundances", nargs="+")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--tpm", required=True)
    parser.add_argument("--lengths", required=True)
    parser.add_argument("--tx2gene", help="two-column TSV: transcript id, gene id")
    args = parser.parse_args()

    mapping = read_tx2gene(args.tx2gene) if args.tx2gene else None

    samples, counts, tpms = [], {}, {}
    features, seen = [], set()
    lengths = {}

    for path in sorted(args.abundances):
        name = sample_name(path)
        if name in counts:
            sys.exit(f"[merge_kallisto] two files map to sample '{name}'")
        c, t, eff = read_abundance(path)
        if not c:
            sys.exit(f"[merge_kallisto] {path} holds no transcript rows")

        if mapping is not None:
            unmapped = [tx for tx in c if tx not in mapping]
            if len(unmapped) == len(c):
                sys.exit(f"[merge_kallisto] none of the {len(c)} transcript ids in "
                         f"{path} appear in --tx2gene; the map is for a different "
                         f"transcriptome (e.g. '{next(iter(c))}' not found)")
            agg_c, agg_t, agg_len = {}, {}, {}
            for tx, value in c.items():
                gene = mapping.get(tx)
                if gene is None:
                    continue
                agg_c[gene] = agg_c.get(gene, 0.0) + value
                agg_t[gene] = agg_t.get(gene, 0.0) + t.get(tx, 0.0)
                # Effective length of a gene is not additive over its
                # transcripts; the abundance-weighted mean is what tximport uses
                # and is the only summary that stays comparable across samples.
                agg_len.setdefault(gene, []).append((value, eff.get(tx, 0.0)))
            c, t = agg_c, agg_t
            eff = {}
            for gene, pairs in agg_len.items():
                total = sum(w for w, _ in pairs)
                eff[gene] = (sum(w * l for w, l in pairs) / total) if total else (
                    sum(l for _, l in pairs) / len(pairs))

        samples.append(name)
        counts[name] = c
        tpms[name] = t
        for feature in c:
            if feature not in seen:
                seen.add(feature)
                features.append(feature)
        for feature, value in eff.items():
            lengths.setdefault(feature, []).append(value)

    write_matrix(args.counts, features, samples, counts)
    write_matrix(args.tpm, features, samples, tpms)

    with open(args.lengths, "w", encoding="utf-8") as handle:
        handle.write("gene_id\teff_length\n")
        for feature in features:
            values = lengths.get(feature, [])
            mean = sum(values) / len(values) if values else 0.0
            handle.write(f"{feature}\t{mean:.2f}\n")

    level = "genes" if mapping else "transcripts"
    print(f"[merge_kallisto] {len(features)} {level} x {len(samples)} samples "
          f"-> {args.counts}", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
db_prepare_genome.py -- one genome, made ready for a classifier's database.

Kraken2 and KrakenUniq both place a sequence on the taxonomy by its taxid, and
both can take that taxid from the sequence itself rather than from NCBI's
accession-to-taxid tables (which run to tens of gigabytes and know nothing of a
genome that is not in GenBank). So every header is rewritten to carry it:

    >NZ_CP000001.1 Escherichia coli ...   becomes
    >ecoli__NZ_CP000001.1|kraken:taxid|562 Escherichia coli ...

`kraken:taxid|<n>` is the form Kraken2's manual documents for --add-to-library.
The genome id is prefixed to the sequence id because two genomes may share
contig names (`contig_1`), and KrakenUniq's seqid-to-taxid map is keyed on them.
The same file therefore serves both builds, plus a two-column `<file>.map`
(seqid, taxid), which is what krakenuniq-build reads next to each library file.

sylph needs nothing of this: it sketches genomes as files, and the taxonomy is
attached afterwards by file name (bin/db_lineage.py).
"""
import argparse
import gzip
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="genome FASTA, plain or gzipped")
    parser.add_argument("--id", required=True, help="genome id from the sheet (names the output)")
    parser.add_argument("--taxid", required=True, help="NCBI taxid every sequence is assigned to")
    args = parser.parse_args()

    opener = gzip.open if args.fasta.endswith(".gz") else open
    out_fasta = f"{args.id}.fna"
    sequences = 0
    bases = 0
    with opener(args.fasta, "rt", encoding="utf-8", errors="replace") as handle, \
            open(out_fasta, "w", encoding="utf-8") as fasta, \
            open(f"{out_fasta}.map", "w", encoding="utf-8") as mapping:
        for line in handle:
            if line.startswith(">"):
                fields = line[1:].rstrip("\n").split(None, 1)
                if not fields:
                    sys.exit(f"[db_prepare_genome] {args.fasta}: a header with no sequence id.")
                seqid = f"{args.id}__{fields[0]}|kraken:taxid|{args.taxid}"
                rest = f" {fields[1]}" if len(fields) > 1 else ""
                fasta.write(f">{seqid}{rest}\n")
                mapping.write(f"{seqid}\t{args.taxid}\n")
                sequences += 1
            else:
                fasta.write(line)
                bases += len(line.strip())

    if sequences == 0 or bases == 0:
        sys.exit(f"[db_prepare_genome] {args.fasta} holds no sequence: is it a FASTA file?")
    print(f"[db_prepare_genome] {args.id}: {sequences} sequence(s), {bases} bp, taxid {args.taxid}.", file=sys.stderr)


if __name__ == "__main__":
    main()

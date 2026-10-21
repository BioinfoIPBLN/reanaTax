#!/usr/bin/env python3
"""
mpa_uniqueness.py -- which genomes each MetaPhlAn marker candidate also occurs in.

MetaPhlAn 4 tests a candidate by cutting it into 150-nt fragments, the size of
a read, and mapping them with Bowtie2 against the genomes of every species; a
hit of any fragment is a hit of the gene (Blanco-Miguez et al., Nat Biotechnol
2023, Methods). Two subcommands around that Bowtie2 run:

  prepare  writes the fragments, and the target genomes plus any --screen
           sequences as one FASTA. Screen sequences (host, vector, anything a
           read could come from but that is not in the database) are renamed
           `__screen__<n>`: a marker they contain would credit their reads to a
           microbe.
  hits     reads Bowtie2's SAM on stdin and writes each (candidate, genome) pair
           once. A contig name's genome is its prefix up to `__`, as
           db_prepare_genome.py wrote it; genome ids cannot contain `__`, so a
           screen hit, whose prefix is empty, is written as `__screen__`.

A trailing piece shorter than --min-tail is dropped rather than mapped: it would
align by chance more readily than a full-length fragment.
"""
import argparse
import gzip
import sys


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def prepare(args):
    count = 0
    with open(args.fragments_out, "w", encoding="utf-8") as handle:
        for name, sequence in read_fasta(args.candidates):
            candidate = name.split()[0]
            for index, start in enumerate(range(0, len(sequence), args.fragment)):
                piece = sequence[start:start + args.fragment]
                if len(piece) < min(args.fragment, args.min_tail):
                    continue
                handle.write(f">{candidate}#{index}\n{piece}\n")
                count += 1

    with open(args.targets, encoding="utf-8") as handle:
        targets = [line.split("\t")[0] for line in handle if line.strip()]
    screened = 0
    with open(args.targets_out, "w", encoding="utf-8") as out:
        for genome in targets:
            for name, sequence in read_fasta(f"{args.genomes_dir}/{genome}.fna"):
                out.write(f">{name}\n{sequence}\n")
        for path in args.screen:
            for _name, sequence in read_fasta(path):
                screened += 1
                out.write(f">__screen__{screened}\n{sequence}\n")
    print(f"[mpa_uniqueness] {count} fragment(s); {len(targets)} target genome(s), "
          f"{screened} screen sequence(s).", file=sys.stderr)


def hits(args):
    pairs = set()
    for line in sys.stdin:
        if line.startswith("@"):
            continue
        fields = line.split("\t", 4)
        if len(fields) < 3 or int(fields[1]) & 4:
            continue
        candidate = fields[0].rsplit("#", 1)[0]
        genome = fields[2].split("__", 1)[0] or "__screen__"
        pairs.add((candidate, genome))
    with open(args.output, "w", encoding="utf-8") as handle:
        for candidate, genome in sorted(pairs):
            handle.write(f"{candidate}\t{genome}\n")
    print(f"[mpa_uniqueness] {len(pairs)} candidate-genome hit(s).", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--candidates", required=True)
    p.add_argument("--targets", required=True, help="targets.tsv: genome, taxid")
    p.add_argument("--genomes-dir", required=True, help="directory of <genome>.fna")
    p.add_argument("--screen", nargs="*", default=[])
    p.add_argument("--fragment", type=int, default=150)
    p.add_argument("--min-tail", type=int, default=75)
    p.add_argument("--fragments-out", required=True)
    p.add_argument("--targets-out", required=True)
    h = sub.add_parser("hits")
    h.add_argument("--output", required=True)
    args = parser.parse_args()
    prepare(args) if args.command == "prepare" else hits(args)


if __name__ == "__main__":
    main()

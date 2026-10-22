#!/usr/bin/env python3
"""
host_mask_kmers.py -- the protected genomes cut into overlapping k-mers, for
masking the host genome against them.

Hostile's masking (Constantinides et al., Bioinformatics 2023; `hostile mask`):
every genome whose reads must survive host depletion is cut into k-mers of a
read's length (150 nt), one every --step bases (10), the k-mers are aligned to
the host genome, and every host region they reach is replaced with N. A read of
the protected organism then has no host sequence left to align to. A k-mer is
named <sequence>_<offset>; sequences shorter than k give none.
"""
import argparse
import gzip
import sys


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = (line[1:].split() or [""])[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--targets", nargs="+", required=True, help="genomes to protect, FASTA, plain or gzipped")
    parser.add_argument("--length", type=int, default=150)
    parser.add_argument("--step", type=int, default=10)
    parser.add_argument("--output", required=True, help="k-mer FASTA, or - for standard output")
    args = parser.parse_args()
    if args.length < 1 or args.step < 1:
        sys.exit("[host_mask_kmers] --length and --step must be positive.")

    sequences = kmers = 0
    out = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    with out:
        for path in args.targets:
            for name, sequence in read_fasta(path):
                sequences += 1
                for offset in range(0, len(sequence) - args.length + 1, args.step):
                    out.write(f">{name}_{offset}\n{sequence[offset:offset + args.length]}\n")
                    kmers += 1
    if not kmers:
        sys.exit(f"[host_mask_kmers] no {args.length}-mer from {sequences} sequence(s): are the --host_mask "
                 "files FASTA, and is any sequence at least that long?")
    print(f"[host_mask_kmers] {kmers} {args.length}-mers from {sequences} sequence(s), step {args.step}.", file=sys.stderr)


if __name__ == "__main__":
    main()

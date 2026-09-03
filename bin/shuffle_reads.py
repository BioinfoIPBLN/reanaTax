#!/usr/bin/env python3
"""
shuffle_reads.py -- destroy the sequence, keep the composition.

Every other filter in this pipeline asks whether a taxon's evidence is strong.
This one asks a different question: how much evidence would this database have
produced from reads that contain no biological sequence at all?

A k-mer classifier matches exact substrings, so its false positives are driven
by base composition. An AT-rich read finds AT-rich genomes; a poly-G tail from
two-colour chemistry finds whatever GC-rich organism happens to be in the
database. Shuffling a read leaves its length, its GC content and its
dinucleotide frequencies exactly where they were and destroys every 35-mer it
contained, so the taxa that still collect reads afterwards are being called on
composition alone. Their real counts are not evidence of anything.

This is a NEGATIVE CONTROL, not a filter of the usual kind. It does not judge a
taxon against a threshold that someone chose; it measures how many reads that
taxon attracts from noise with the same composition as the real library, and
that number is directly comparable to its real count.

Three shuffles, in decreasing order of how conservative they are:

  dinuc    Altschul-Erickson: an Euler-path walk that preserves the exact
           dinucleotide frequency as well as the base composition. The default,
           and the honest one. Most compositional bias in a genome is
           dinucleotide-level (CpG depletion above all), so a shuffle that
           discards it understates how many chance matches real data produces.
  mono     a plain permutation of the bases. Preserves GC content and nothing
           else, so it is the weakest null and will flag fewer taxa.
  reverse  the read read backwards - NOT reverse-complemented, which would be a
           real sequence and would map. Deterministic, free, and preserves
           every compositional statistic there is, including the trinucleotide
           spectrum. Its weakness is that a palindromic or low-complexity
           stretch survives it unchanged, so a homopolymer-driven false
           positive stays a false positive - which is a reason to use it, if
           that is the artefact being chased.

Quality strings are carried through unchanged for the two permutation methods
(Kraken2 does not read them) and reversed alongside the sequence for `reverse`,
so the record stays internally consistent for anything downstream that does.

A fixed seed is required, not optional. A control whose result changes between
runs of the same pipeline is not a control.
"""
import argparse
import gzip
import random
import sys


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode, errors="replace") if "t" in mode else gzip.open(path, mode)
    return open(path, mode, encoding="utf-8", errors="replace")


def fastq_records(path):
    """Yield (name, sequence, plus, quality). Raises on a truncated record."""
    with open_maybe_gzip(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            if not quality:
                raise SystemExit(
                    f"[shuffle_reads] '{path}' ends mid-record. A truncated FASTQ would make the "
                    "control disagree with the real classification for reasons that have nothing "
                    "to do with composition."
                )
            yield name.rstrip("\n"), sequence.rstrip("\n"), plus.rstrip("\n"), quality.rstrip("\n")


def mono_shuffle(sequence, rng):
    bases = list(sequence)
    rng.shuffle(bases)
    return "".join(bases)


def dinuc_shuffle(sequence, rng):
    """Altschul-Erickson dinucleotide-preserving shuffle.

    The sequence is an Eulerian path through a multigraph whose vertices are the
    symbols and whose edges are the observed dinucleotides. Any other Eulerian
    path with the same first and last symbol has, by construction, the identical
    dinucleotide count - so the problem is to draw one uniformly. The standard
    construction: pick a random spanning tree of the vertices oriented towards
    the last symbol, shuffle each vertex's outgoing edges freely, then move the
    edge that belongs to the tree to the end of that vertex's list. That last
    step is what guarantees the walk cannot strand itself, and is the whole
    content of the algorithm.
    """
    if len(sequence) < 3:
        return sequence

    first, last = sequence[0], sequence[-1]
    edges = {}
    for source, target in zip(sequence, sequence[1:]):
        edges.setdefault(source, []).append(target)

    vertices = sorted(set(sequence))
    if len(vertices) == 1:
        return sequence

    # A random arborescence towards `last`: for every other vertex, keep drawing
    # one of its outgoing edges until the walk it starts reaches `last`.
    for _attempt in range(100):
        tree_edge = {}
        for vertex in vertices:
            if vertex == last or not edges.get(vertex):
                continue
            tree_edge[vertex] = rng.choice(edges[vertex])
        if _reaches(tree_edge, vertices, last):
            break
    else:
        # Vanishingly unlikely, and a wrong answer here would be a silent one.
        return mono_shuffle(sequence, rng)

    walk_edges = {}
    for vertex, targets in edges.items():
        remaining = list(targets)
        chosen = tree_edge.get(vertex)
        if chosen is not None:
            remaining.remove(chosen)
        rng.shuffle(remaining)
        if chosen is not None:
            remaining.append(chosen)
        walk_edges[vertex] = remaining

    positions = dict.fromkeys(walk_edges, 0)
    out = [first]
    current = first
    for _step in range(len(sequence) - 1):
        index = positions[current]
        following = walk_edges[current][index]
        positions[current] = index + 1
        out.append(following)
        current = following
    return "".join(out)


def _reaches(tree_edge, vertices, last):
    """True when every vertex with a tree edge walks into `last`."""
    for vertex in vertices:
        if vertex == last or vertex not in tree_edge:
            continue
        seen = set()
        current = vertex
        while current != last:
            if current in seen or current not in tree_edge:
                return False
            seen.add(current)
            current = tree_edge[current]
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--reads", nargs="+", required=True, help="one or two FASTQ files")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--method", default="dinuc", choices=("dinuc", "mono", "reverse"))
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument(
        "--max-reads",
        type=int,
        default=0,
        help="shuffle only the first N reads (pairs); 0 takes the whole library",
    )
    parser.add_argument("--stats", default=None, help="write a one-line accounting TSV here")
    args = parser.parse_args()

    if len(args.reads) > 2:
        sys.exit(f"[shuffle_reads] expected one or two FASTQ files, got {len(args.reads)}.")

    rng = random.Random(args.seed)
    shuffler = {"mono": mono_shuffle, "dinuc": dinuc_shuffle}.get(args.method)

    written = 0
    bases = 0
    gc = 0
    handles = []
    try:
        for index, _path in enumerate(args.reads, start=1):
            suffix = f"_{index}" if len(args.reads) == 2 else ""
            handles.append(gzip.open(f"{args.prefix}_shuffled{suffix}.fastq.gz", "wt"))

        # Mates are stepped together so that truncation by --max-reads keeps the
        # pairing: Kraken2 reads a pair as one fragment and a half-truncated
        # file would silently change what the control is counting.
        streams = [fastq_records(path) for path in args.reads]
        while True:
            records = []
            for stream in streams:
                records.append(next(stream, None))
            if records[0] is None:
                break
            if any(record is None for record in records):
                sys.exit("[shuffle_reads] the mate files hold different numbers of reads.")
            if args.max_reads and written >= args.max_reads:
                break
            for handle, (name, sequence, plus, quality) in zip(handles, records):
                upper = sequence.upper()
                bases += len(upper)
                gc += upper.count("G") + upper.count("C")
                if args.method == "reverse":
                    new_sequence = sequence[::-1]
                    new_quality = quality[::-1]
                else:
                    new_sequence = shuffler(upper, rng)
                    new_quality = quality
                handle.write(f"{name}\n{new_sequence}\n{plus}\n{new_quality}\n")
            written += 1
    finally:
        for handle in handles:
            handle.close()

    if args.stats:
        with open(args.stats, "w", encoding="utf-8") as handle:
            handle.write("sample\tmethod\tseed\tfragments_shuffled\tbases\tgc_fraction\n")
            handle.write(
                f"{args.prefix}\t{args.method}\t{args.seed}\t{written}\t{bases}\t"
                f"{gc / bases if bases else 0:.6f}\n"
            )

    print(
        f"[shuffle_reads] {args.prefix}: {written} fragment(s) shuffled by '{args.method}' "
        f"(seed {args.seed}), GC {100 * gc / bases if bases else 0:.2f}%.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

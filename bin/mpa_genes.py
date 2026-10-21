#!/usr/bin/env python3
"""
mpa_genes.py -- one genome's protein-coding genes, as candidate MetaPhlAn markers.

MetaPhlAn 4 markers are core protein-coding genes of 450 to 4,500 nt
(Blanco-Miguez et al., Nat Biotechnol 2023, Methods), so that is what is written
here, as nucleotide sequences for the markers themselves and as proteins for the
clustering that decides which genes a species holds in common.

Genes come from one of two places:

  * --cds, a coding-sequence FASTA from the genome's own annotation. Required
    for eukaryotes: pyrodigal is a prokaryotic gene finder and splits every
    intron-containing gene. VEuPathDB, Ensembl and NCBI all distribute one.
  * pyrodigal otherwise. Trained on the genome when it holds 100 kb or more,
    run in metagenomic mode below that, as Prodigal's authors advise. Genetic
    code 4 is chosen over 11 when it codes markedly more of the genome - the
    rule CheckM applies (density 5 points higher, and code 11 under 70%) -
    because Mycoplasma and Spiroplasma read TGA as tryptophan, and under code
    11 their genes come out in fragments. Genes running off a contig edge are
    left out: a truncated gene makes a truncated marker.

Writes <id>.genes.fna, <id>.genes.faa, <id>.genes.tsv (gene, genome, length) and
<id>.genome.tsv (genome, taxid, length, contigs, genes, source, genetic code).
"""
import argparse
import gzip
import sys

CODONS = {
    a + b + c: aa
    for (a, b, c), aa in zip(
        ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
        "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    )
}


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, chunks = None, []
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split(None, 1)[0] if len(line) > 1 else "", []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def translate(sequence):
    """Standard-code translation; None for a sequence that is not one clean ORF."""
    if len(sequence) % 3:
        return None
    protein = "".join(CODONS.get(sequence[i:i + 3], "X") for i in range(0, len(sequence), 3))
    if protein.endswith("*"):
        protein = protein[:-1]
    return None if "*" in protein or not protein else protein


def genes_from_cds(path):
    for _name, sequence in read_fasta(path):
        protein = translate(sequence)
        if protein is not None:
            yield sequence, protein


def genes_from_pyrodigal(contigs):
    import pyrodigal

    total = sum(len(sequence) for _name, sequence in contigs)
    if total >= 100000:
        sequences = [sequence.encode() for _name, sequence in contigs]
        finders = {}
        density = {}
        for code in (11, 4):
            finder = pyrodigal.GeneFinder()
            finder.train(*sequences, translation_table=code)
            coding = sum(len(gene.sequence()) for sequence in sequences for gene in finder.find_genes(sequence))
            finders[code], density[code] = finder, coding / total
        code = 4 if density[4] > density[11] + 0.05 and density[11] < 0.7 else 11
        finder = finders[code]
    else:
        code = 11
        finder = pyrodigal.GeneFinder(meta=True)
    genes = []
    for _name, sequence in contigs:
        for gene in finder.find_genes(sequence.encode()):
            if gene.partial_begin or gene.partial_end:
                continue
            genes.append((gene.sequence(), gene.translate(include_stop=False)))
    return genes, code


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True, help="the genome")
    parser.add_argument("--id", required=True, help="genome id from the sheet")
    parser.add_argument("--taxid", required=True)
    parser.add_argument("--cds", default=None, help="coding sequences from the genome's annotation")
    parser.add_argument("--min-length", type=int, default=450)
    parser.add_argument("--max-length", type=int, default=4500)
    args = parser.parse_args()

    contigs = list(read_fasta(args.fasta))
    length = sum(len(sequence) for _name, sequence in contigs)
    if not contigs or not length:
        sys.exit(f"[mpa_genes] {args.fasta} holds no sequence.")

    if args.cds:
        genes, source, code = list(genes_from_cds(args.cds)), "cds", 1
        if not genes:
            sys.exit(f"[mpa_genes] {args.cds} holds no complete coding sequence (a length divisible by three, "
                     "no internal stop codon). Is it a nucleotide CDS FASTA?")
    else:
        genes, code = genes_from_pyrodigal(contigs)
        source = "pyrodigal"

    kept = 0
    with open(f"{args.id}.genes.fna", "w", encoding="utf-8") as fna, \
            open(f"{args.id}.genes.faa", "w", encoding="utf-8") as faa, \
            open(f"{args.id}.genes.tsv", "w", encoding="utf-8") as table:
        for nucleotides, protein in genes:
            if not args.min_length <= len(nucleotides) <= args.max_length:
                continue
            kept += 1
            gene = f"{args.id}__g{kept}"
            fna.write(f">{gene}\n{nucleotides}\n")
            faa.write(f">{gene}\n{protein}\n")
            table.write(f"{gene}\t{args.id}\t{len(nucleotides)}\n")

    with open(f"{args.id}.genome.tsv", "w", encoding="utf-8") as handle:
        handle.write(f"{args.id}\t{args.taxid}\t{length}\t{len(contigs)}\t{kept}\t{source}\t{code}\n")
    print(f"[mpa_genes] {args.id}: {kept} of {len(genes)} {source} gene(s) within "
          f"{args.min_length}-{args.max_length} nt.", file=sys.stderr)


if __name__ == "__main__":
    main()

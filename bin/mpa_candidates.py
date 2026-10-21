#!/usr/bin/env python3
"""
mpa_candidates.py -- the core genes of each species: MetaPhlAn marker candidates.

Follows the MetaPhlAn 4 recipe (Blanco-Miguez et al., Nat Biotechnol 2023,
Methods), as HROM reimplemented it (github.com/netbiolab/HROM,
7.Marker_database). A species is every genome of the sheet that shares a taxid.

  * A protein cluster (MMseqs2, 90% identity, 80% coverage) is core to a species
    when it holds a gene from at least the coreness threshold of the species'
    genomes. The threshold is the highest that still leaves 800 core clusters,
    but never below 0.6 for a species of fewer than 100 genomes or 0.5 for a
    larger one. A species of one genome has every gene at coreness 1.
  * A cluster found in more than 10 species is nobody's marker (HROM's
    prefilter): it is a conserved family, and the fragment mapping that follows
    would reject it anyway at a far higher cost.
  * Each surviving (species, cluster) pair becomes one candidate, represented
    by the species' member gene of median length.

Up to --max-genomes genomes per species are the mapping targets of the
uniqueness test, as MetaPhlAn maps against at most 100 per SGB. The most
contiguous ones are taken: fewest contigs, then longest.

Species absent from the MetaPhlAn taxonomy table (viruses; see db_lineage.py)
are left out.
"""
import argparse
import csv
import sys

MAX_SPECIES_PER_CLUSTER = 10
TARGET_CORE = 800


def coreness_threshold(values, n_genomes):
    """MetaPhlAn 4's rule: the highest threshold leaving TARGET_CORE clusters, floored."""
    floor = 0.6 if n_genomes < 100 else 0.5
    ranked = sorted(values, reverse=True)
    best = ranked[TARGET_CORE - 1] if len(ranked) >= TARGET_CORE else 0.0
    return max(floor, best)


def read_fasta_ids(paths, wanted):
    sequences = {}
    for path in paths:
        name, chunks = None, []
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if line.startswith(">"):
                    if name in wanted:
                        sequences[name] = "".join(chunks)
                    name, chunks = line[1:].split()[0], []
                elif name in wanted:
                    chunks.append(line)
        if name in wanted:
            sequences[name] = "".join(chunks)
    return sequences


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--clusters", required=True, help="MMseqs2 cluster TSV: representative, member")
    parser.add_argument("--genes", nargs="+", required=True, help="*.genes.tsv: gene, genome, length")
    parser.add_argument("--genomes", nargs="+", required=True, help="*.genome.tsv from mpa_genes.py")
    parser.add_argument("--sequences", nargs="+", required=True, help="*.genes.fna")
    parser.add_argument("--taxonomy", required=True, help="MetaPhlAn taxonomy table from db_lineage.py")
    parser.add_argument("--max-genomes", type=int, default=100)
    args = parser.parse_args()

    with open(args.taxonomy, encoding="utf-8") as handle:
        allowed = {row["taxid"] for row in csv.DictReader(handle, delimiter="\t")}

    genomes = {}
    for path in args.genomes:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                genome, taxid, length, contigs, n_genes, _source, _code = line.rstrip("\n").split("\t")
                genomes[genome] = (taxid, int(length), int(contigs), int(n_genes))
    species = {}
    for genome, (taxid, _length, _contigs, _genes) in genomes.items():
        if taxid in allowed:
            species.setdefault(taxid, []).append(genome)
    if not species:
        sys.exit("[mpa_candidates] no genome belongs to a species MetaPhlAn can profile (viruses are left out).")

    gene_genome, gene_length = {}, {}
    for path in args.genes:
        with open(path, encoding="utf-8") as handle:
            for line in handle:
                gene, genome, length = line.rstrip("\n").split("\t")
                gene_genome[gene], gene_length[gene] = genome, int(length)

    # cluster -> species -> {genome: [genes]}
    clusters = {}
    with open(args.clusters, encoding="utf-8") as handle:
        for line in handle:
            representative, member = line.rstrip("\n").split("\t")[:2]
            genome = gene_genome.get(member)
            if genome is None:
                continue
            taxid = genomes[genome][0]
            if taxid not in species:
                continue
            clusters.setdefault(representative, {}).setdefault(taxid, {}).setdefault(genome, []).append(member)

    per_species = {taxid: [] for taxid in species}
    for representative, by_species in clusters.items():
        if len(by_species) > MAX_SPECIES_PER_CLUSTER:
            continue
        for taxid, by_genome in by_species.items():
            coreness = len(by_genome) / len(species[taxid])
            per_species[taxid].append((representative, coreness, [g for genes in by_genome.values() for g in genes]))

    chosen = {}
    bins = []
    for taxid in sorted(species, key=lambda t: int(t) if t.isdigit() else t):
        members = species[taxid]
        threshold = coreness_threshold([coreness for _rep, coreness, _genes in per_species[taxid]], len(members))
        core = [(rep, coreness, genes) for rep, coreness, genes in per_species[taxid] if coreness >= threshold]
        for index, (rep, coreness, genes) in enumerate(sorted(core), start=1):
            genes = sorted(genes, key=lambda gene: (gene_length[gene], gene))
            chosen[f"TAX{taxid}_C{index}"] = (genes[len(genes) // 2], taxid, rep, coreness, threshold)
        targets = sorted(members, key=lambda genome: (genomes[genome][2], -genomes[genome][1], genome))[: args.max_genomes]
        bins.append((taxid, len(members), targets, threshold, len(core)))

    sequences = read_fasta_ids(args.sequences, {gene for gene, *_rest in chosen.values()})
    with open("candidates.fna", "w", encoding="utf-8") as fasta, \
            open("candidates.tsv", "w", encoding="utf-8") as table:
        table.write("candidate\ttaxid\tcluster\tgene\tlength\tcoreness\tthreshold\n")
        for candidate, (gene, taxid, rep, coreness, threshold) in chosen.items():
            fasta.write(f">{candidate}\n{sequences[gene]}\n")
            table.write(f"{candidate}\t{taxid}\t{rep}\t{gene}\t{gene_length[gene]}\t{coreness:.4f}\t{threshold:.4f}\n")
    with open("targets.tsv", "w", encoding="utf-8") as handle:
        for taxid, _n, targets, _threshold, _core in bins:
            for genome in targets:
                handle.write(f"{genome}\t{taxid}\n")
    with open("species.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxid\tgenomes\tmapping_targets\tmean_genome_length\tcoreness_threshold\tcore_clusters\n")
        for taxid, n_genomes, targets, threshold, n_core in bins:
            mean_length = sum(genomes[genome][1] for genome in species[taxid]) / n_genomes
            handle.write(f"{taxid}\t{n_genomes}\t{len(targets)}\t{mean_length:.0f}\t{threshold:.4f}\t{n_core}\n")
    print(f"[mpa_candidates] {len(chosen)} candidate(s) across {len(bins)} species.", file=sys.stderr)


if __name__ == "__main__":
    main()

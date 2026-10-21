#!/usr/bin/env python3
"""
mpa_select.py -- choose each species' markers and write the MetaPhlAn database.

A candidate (mpa_candidates.py) becomes a marker when its fragments
(mpa_uniqueness.py):

  * hit at least the species' coreness threshold of its own mapping targets -
    MetaPhlAn's check that the gene really is in almost every genome of the
    species, at the nucleotide level the reads will be mapped at;
  * hit no screen sequence (host, vector);
  * hit fewer than --max-other of the targets of every other species. With up to
    100 targets per species that means none, unless a species has more; a hit
    below the limit makes a quasi-marker, its species listed in `ext` so that
    MetaPhlAn discounts the marker when that species is present.

Up to --max-markers per species are kept, perfectly unique ones first and then
the longest, as in MetaPhlAn 4. A species left with fewer than --min-markers is
not written into the database: MetaPhlAn discards those too, and a species
judged on a handful of genes is easily called from a single shared locus.

The database is what MetaPhlAn 4.2 reads (metaphlan/utils/database_controller.py):
<name>.pkl, a bz2-compressed pickle of

    {"taxonomy": {lineage: (taxid lineage, genome length)},
     "markers": {marker: {"clade": "t__...", "ext": [...], "len": n, "taxon": lineage}},
     "merged_taxon": {}}

beside <name>.fna, the markers, which the caller indexes with Bowtie2.
"""
import argparse
import bz2
import csv
import pickle
import sys


def read_fasta(path):
    name, chunks = None, []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def read_tsv(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--candidates", required=True, help="candidates.tsv")
    parser.add_argument("--sequences", required=True, help="candidates.fna")
    parser.add_argument("--targets", required=True, help="targets.tsv: genome, taxid")
    parser.add_argument("--species", required=True, help="species.tsv")
    parser.add_argument("--hits", required=True, help="candidate, genome pairs")
    parser.add_argument("--taxonomy", required=True, help="MetaPhlAn taxonomy table from db_lineage.py")
    parser.add_argument("--name", required=True, help="database name (the MetaPhlAn --index)")
    parser.add_argument("--max-markers", type=int, default=200)
    parser.add_argument("--min-markers", type=int, default=10)
    parser.add_argument("--max-other", type=float, default=0.01,
                        help="share of another species' targets a marker may hit (quasi-marker)")
    args = parser.parse_args()

    taxonomy = {row["taxid"]: row for row in read_tsv(args.taxonomy)}
    species = {row["taxid"]: row for row in read_tsv(args.species)}
    targets = {}
    with open(args.targets, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                genome, taxid = line.rstrip("\n").split("\t")[:2]
                targets.setdefault(taxid, set()).add(genome)
    genome_species = {genome: taxid for taxid, genomes in targets.items() for genome in genomes}
    hits = {}
    with open(args.hits, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                candidate, genome = line.rstrip("\n").split("\t")[:2]
                hits.setdefault(candidate, set()).add(genome)

    verdicts = {taxid: {"lack_coreness": 0, "shared": 0, "screen": 0} for taxid in species}
    accepted = {taxid: [] for taxid in species}
    for row in read_tsv(args.candidates):
        candidate, taxid = row["candidate"], row["taxid"]
        hit = hits.get(candidate, set())
        if "__screen__" in hit:
            verdicts[taxid]["screen"] += 1
            continue
        own = len(hit & targets[taxid]) / len(targets[taxid])
        if own < float(row["threshold"]) - 1e-9:
            verdicts[taxid]["lack_coreness"] += 1
            continue
        ext, shared = [], False
        for other in sorted({genome_species[g] for g in hit if g in genome_species} - {taxid}):
            share = len(hit & targets[other]) / len(targets[other])
            if share >= args.max_other:
                shared = True
                break
            ext.append(other)
        if shared:
            verdicts[taxid]["shared"] += 1
            continue
        accepted[taxid].append((len(ext), -int(row["length"]), candidate, ext))

    sequences = dict(read_fasta(args.sequences))
    label = {taxid: f"TAX{taxid}" for taxid in species}
    chosen_by = {taxid: sorted(accepted[taxid])[: args.max_markers] for taxid in species}
    kept_species = {taxid for taxid, chosen in chosen_by.items() if len(chosen) >= args.min_markers}
    database = {"taxonomy": {}, "markers": {}, "merged_taxon": {}}
    summary = []
    with open(f"{args.name}.fna", "w", encoding="utf-8") as fasta:
        for taxid in sorted(species, key=lambda t: int(t) if t.isdigit() else t):
            chosen = chosen_by[taxid]
            status = "kept" if taxid in kept_species else f"dropped: {len(chosen)} marker(s), --min-markers {args.min_markers}"
            if status == "kept":
                lineage = taxonomy[taxid]["lineage"]
                database["taxonomy"][lineage] = (taxonomy[taxid]["taxid_lineage"], int(float(species[taxid]["mean_genome_length"])))
                for number, (_n_ext, _length, candidate, ext) in enumerate(chosen, start=1):
                    marker = f"{label[taxid]}_M{number}"
                    database["markers"][marker] = {
                        "clade": f"t__{label[taxid]}",
                        # MetaPhlAn looks every ext species up in its tree, so
                        # one left out of the database cannot be named here.
                        "ext": [label[other] for other in ext if other in kept_species],
                        "len": len(sequences[candidate]),
                        "taxon": lineage,
                    }
                    fasta.write(f">{marker}\n{sequences[candidate]}\n")
            counts = verdicts[taxid]
            summary.append([taxid, taxonomy[taxid]["lineage"], species[taxid]["genomes"],
                            species[taxid]["mapping_targets"], species[taxid]["coreness_threshold"],
                            species[taxid]["core_clusters"], counts["lack_coreness"], counts["shared"],
                            counts["screen"], len(chosen), status])

    with open(f"{args.name}.mpa_species.tsv", "w", encoding="utf-8") as handle:
        handle.write("taxid\tlineage\tgenomes\tmapping_targets\tcoreness_threshold\tcore_clusters\t"
                     "rejected_coreness\trejected_shared\trejected_screen\tmarkers\tstatus\n")
        for entry in summary:
            handle.write("\t".join(str(value) for value in entry) + "\n")

    if not database["markers"]:
        sys.exit(f"[mpa_select] no species kept {args.min_markers} or more markers; see {args.name}.mpa_species.tsv.")
    with bz2.BZ2File(f"{args.name}.pkl", "w") as handle:
        pickle.dump(database, handle, protocol=2)
    kept = sum(1 for entry in summary if entry[-1] == "kept")
    print(f"[mpa_select] {len(database['markers'])} marker(s) for {kept} of {len(summary)} species.", file=sys.stderr)


if __name__ == "__main__":
    main()

"""The MetaPhlAn database build: mpa_genes, mpa_candidates, mpa_uniqueness,
mpa_select and db_lineage --mpa-output."""

import bz2
import csv
import pickle


def write(path, text):
    path.write_text(text, encoding="utf-8")
    return path


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_genes_from_cds_keep_complete_orfs_of_marker_length(tmp_path, run_script):
    """With --cds, a gene is kept when it is one clean ORF (length divisible by
    three, no internal stop) of 450-4500 nt; its protein drops the stop codon."""
    good = "ATG" + "GCT" * 150 + "TAA"
    short = "ATG" + "GCT" * 98 + "TAA"
    internal_stop = "ATG" + "GCT" * 74 + "TAA" + "GCT" * 75 + "TAA"
    frameshift = good + "A"
    cds = write(tmp_path / "g1.cds.fna", "".join(
        f">{name}\n{seq}\n" for name, seq in
        [("good", good), ("short", short), ("stop", internal_stop), ("shift", frameshift)]))
    genome = write(tmp_path / "g1.fna", ">contig1\n" + "ACGT" * 100 + "\n")

    result = run_script("mpa_genes.py", "--fasta", genome, "--id", "g1", "--taxid", 562, "--cds", cds, cwd=tmp_path)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "g1.genes.tsv").read_text(encoding="utf-8").splitlines() == ["g1__g1\tg1\t456"]
    assert (tmp_path / "g1.genes.faa").read_text(encoding="utf-8").splitlines() == [">g1__g1", "M" + "A" * 150]
    assert (tmp_path / "g1.genes.fna").read_text(encoding="utf-8").splitlines() == [">g1__g1", good]
    assert (tmp_path / "g1.genome.tsv").read_text(encoding="utf-8").strip() == "g1\t562\t400\t1\t1\tcds\t1"


def test_candidates_are_core_families_of_each_species(tmp_path, run_script):
    """A family is core to a species at the coreness floor (0.6 under 100
    genomes), a family in more than 10 species is nobody's, the candidate is the
    species' member of median length, and the most contiguous genomes are the
    mapping targets."""
    genomes = {"a1": ("100", 3000, 3), "a2": ("100", 3000, 1), "a3": ("100", 4000, 1), "b1": ("200", 2000, 1)}
    genomes.update({f"x{i}": (str(300 + i), 1000, 1) for i in range(11)})
    genes = {
        "a1__g1": 500, "a2__g1": 510, "a3__g1": 520, "b1__g1": 600,
        "a1__g2": 700, "a2__g2": 700,
        "a1__g3": 900,
    }
    genes.update({f"x{i}__g1": 800 for i in range(11)})
    clusters = (
        [("a1__g1", member) for member in ("a1__g1", "a2__g1", "a3__g1", "b1__g1")]
        + [("a1__g2", "a1__g2"), ("a1__g2", "a2__g2"), ("a1__g3", "a1__g3")]
        + [("x0__g1", f"x{i}__g1") for i in range(11)]
    )
    write(tmp_path / "clusters.tsv", "".join(f"{rep}\t{member}\n" for rep, member in clusters))
    write(tmp_path / "all.genes.tsv", "".join(f"{gene}\t{gene.split('__')[0]}\t{length}\n" for gene, length in genes.items()))
    write(tmp_path / "all.genome.tsv", "".join(
        f"{genome}\t{taxid}\t{length}\t{contigs}\t1\tpyrodigal\t11\n" for genome, (taxid, length, contigs) in genomes.items()))
    write(tmp_path / "all.genes.fna", "".join(f">{gene}\n{'A' * length}\n" for gene, length in genes.items()))
    taxids = sorted({taxid for taxid, _length, _contigs in genomes.values()})
    write(tmp_path / "taxonomy.tsv", "taxid\tlineage\ttaxid_lineage\tdomain\n"
          + "".join(f"{taxid}\tk__Bacteria|t__TAX{taxid}\t2||||||\tBacteria\n" for taxid in taxids))

    result = run_script(
        "mpa_candidates.py",
        "--clusters", tmp_path / "clusters.tsv",
        "--genes", tmp_path / "all.genes.tsv",
        "--genomes", tmp_path / "all.genome.tsv",
        "--sequences", tmp_path / "all.genes.fna",
        "--taxonomy", tmp_path / "taxonomy.tsv",
        "--max-genomes", 2,
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr
    candidates = {(row["taxid"], row["cluster"]): row for row in read_tsv(tmp_path / "candidates.tsv")}
    assert set(candidates) == {("100", "a1__g1"), ("100", "a1__g2"), ("200", "a1__g1")}
    assert candidates[("100", "a1__g1")]["gene"] == "a2__g1"
    assert candidates[("100", "a1__g2")]["coreness"] == "0.6667"
    assert candidates[("100", "a1__g2")]["threshold"] == "0.6000"
    targets = [line.split("\t") for line in (tmp_path / "targets.tsv").read_text(encoding="utf-8").splitlines()]
    assert [genome for genome, taxid in targets if taxid == "100"] == ["a3", "a2"]
    species = {row["taxid"]: row for row in read_tsv(tmp_path / "species.tsv")}
    assert species["100"]["core_clusters"] == "2"
    assert species["300"]["core_clusters"] == "0"
    assert len((tmp_path / "candidates.fna").read_text(encoding="utf-8").splitlines()) == 6


def test_uniqueness_fragments_and_hits(tmp_path, run_script):
    """Candidates are cut into 150-nt fragments (a tail under 75 nt is dropped),
    screen sequences are renamed so no genome id can match them, and each
    candidate-genome hit is written once whatever the number of alignments."""
    write(tmp_path / "candidates.fna", f">TAX100_C1\n{'A' * 400}\n>TAX100_C2\n{'C' * 160}\n")
    write(tmp_path / "targets.tsv", "a1\t100\n")
    (tmp_path / "genomes").mkdir()
    write(tmp_path / "genomes" / "a1.fna", ">a1__contig1|kraken:taxid|562\nACGTACGT\n")
    write(tmp_path / "host.fa", ">chr1\nGGGG\n>chr2\nTTTT\n")

    prepared = run_script(
        "mpa_uniqueness.py", "prepare",
        "--candidates", tmp_path / "candidates.fna",
        "--targets", tmp_path / "targets.tsv",
        "--genomes-dir", tmp_path / "genomes",
        "--screen", tmp_path / "host.fa",
        "--fragments-out", tmp_path / "fragments.fa",
        "--targets-out", tmp_path / "targets.fa",
    )
    assert prepared.returncode == 0, prepared.stderr
    fragments = (tmp_path / "fragments.fa").read_text(encoding="utf-8").splitlines()
    assert fragments[0::2] == [">TAX100_C1#0", ">TAX100_C1#1", ">TAX100_C1#2", ">TAX100_C2#0"]
    assert [len(line) for line in fragments[1::2]] == [150, 150, 100, 150]
    headers = [line for line in (tmp_path / "targets.fa").read_text(encoding="utf-8").splitlines() if line.startswith(">")]
    assert headers == [">a1__contig1|kraken:taxid|562", ">__screen__1", ">__screen__2"]

    sam = "\n".join([
        "@HD\tVN:1.0",
        "TAX100_C1#0\t0\ta1__contig1|kraken:taxid|562\t1\t42\t150M\t*\t0\t0\tA\tI",
        "TAX100_C1#1\t256\ta1__contig1|kraken:taxid|562\t9\t1\t150M\t*\t0\t0\tA\tI",
        "TAX100_C2#0\t0\t__screen__2\t1\t42\t150M\t*\t0\t0\tC\tI",
        "TAX100_C2#1\t4\t*\t0\t0\t*\t*\t0\t0\tC\tI",
    ]) + "\n"
    hits = run_script("mpa_uniqueness.py", "hits", "--output", tmp_path / "hits.tsv", stdin=sam)
    assert hits.returncode == 0, hits.stderr
    assert (tmp_path / "hits.tsv").read_text(encoding="utf-8").splitlines() == [
        "TAX100_C1\ta1", "TAX100_C2\t__screen__"]


def write_select_inputs(tmp_path):
    lengths = {"C1": 600, "C2": 900, "C3": 700, "C4": 650, "C5": 500, "D1": 800}
    owner = {"C1": "100", "C2": "100", "C3": "100", "C4": "100", "C5": "100", "D1": "200"}
    threshold = {"100": 0.6, "200": 1.0}
    write(tmp_path / "candidates.fna", "".join(f">{c}\n{'A' * n}\n" for c, n in lengths.items()))
    write(tmp_path / "candidates.tsv", "candidate\ttaxid\tcluster\tgene\tlength\tcoreness\tthreshold\n" + "".join(
        f"{c}\t{owner[c]}\t{c}\t{c}\t{n}\t1.0\t{threshold[owner[c]]}\n" for c, n in lengths.items()))
    write(tmp_path / "targets.tsv", "a1\t100\na2\t100\nb1\t200\n")
    write(tmp_path / "species.tsv",
          "taxid\tgenomes\tmapping_targets\tmean_genome_length\tcoreness_threshold\tcore_clusters\n"
          "100\t2\t2\t2500000\t0.6\t5\n200\t1\t1\t1800000\t1.0\t1\n")
    hits = {
        "C1": ["a1", "a2"],           # unique
        "C2": ["a1"],                 # half of its own species: lacks coreness
        "C3": ["a1", "a2", "b1"],     # in the other species too
        "C4": ["a1", "a2", "__screen__"],
        "C5": ["a1", "a2"],           # unique
        "D1": ["b1"],
    }
    write(tmp_path / "hits.tsv", "".join(f"{c}\t{g}\n" for c, genomes in hits.items() for g in genomes))
    write(tmp_path / "taxonomy.tsv", "taxid\tlineage\ttaxid_lineage\tdomain\n"
          "100\tk__Bacteria|p__P|c__C|o__O|f__F|g__G|s__G_one|t__TAX100\t2|1|1|1|1|1|100|\tBacteria\n"
          "200\tk__Bacteria|p__P|c__C|o__O|f__F|g__G|s__G_two|t__TAX200\t2|1|1|1|1|1|200|\tBacteria\n")


def run_select(tmp_path, run_script, *extra):
    return run_script(
        "mpa_select.py",
        "--candidates", tmp_path / "candidates.tsv",
        "--sequences", tmp_path / "candidates.fna",
        "--targets", tmp_path / "targets.tsv",
        "--species", tmp_path / "species.tsv",
        "--hits", tmp_path / "hits.tsv",
        "--taxonomy", tmp_path / "taxonomy.tsv",
        "--name", "testdb",
        *extra,
        cwd=tmp_path,
    )


def load_db(tmp_path):
    with bz2.BZ2File(tmp_path / "testdb.pkl", "r") as handle:
        return pickle.load(handle)


def test_select_keeps_unique_core_markers_and_writes_the_metaphlan_pickle(tmp_path, run_script):
    """Markers lacking coreness, shared with another species or found in a screen
    sequence are rejected; the longest unique ones come first; a species under
    --min-markers is left out; and the pickle has MetaPhlAn's three tables."""
    write_select_inputs(tmp_path)
    result = run_select(tmp_path, run_script, "--min-markers", 2)
    assert result.returncode == 0, result.stderr

    db = load_db(tmp_path)
    assert set(db) == {"taxonomy", "markers", "merged_taxon"}
    lineage = "k__Bacteria|p__P|c__C|o__O|f__F|g__G|s__G_one|t__TAX100"
    assert db["taxonomy"] == {lineage: ("2|1|1|1|1|1|100|", 2500000)}
    assert db["markers"] == {
        "TAX100_M1": {"clade": "t__TAX100", "ext": [], "len": 600, "taxon": lineage},
        "TAX100_M2": {"clade": "t__TAX100", "ext": [], "len": 500, "taxon": lineage},
    }
    assert db["merged_taxon"] == {}
    assert (tmp_path / "testdb.fna").read_text(encoding="utf-8").splitlines()[0::2] == [">TAX100_M1", ">TAX100_M2"]

    species = {row["taxid"]: row for row in read_tsv(tmp_path / "testdb.mpa_species.tsv")}
    assert (species["100"]["rejected_coreness"], species["100"]["rejected_shared"],
            species["100"]["rejected_screen"], species["100"]["markers"]) == ("1", "1", "1", "2")
    assert species["200"]["status"].startswith("dropped")


def test_quasi_markers_name_only_species_in_the_database(tmp_path, run_script):
    """Below --max-other a shared candidate is a quasi-marker ranked after the
    unique ones, its other species in `ext` - but only while that species is in
    the database, because MetaPhlAn looks every `ext` entry up."""
    write_select_inputs(tmp_path)
    result = run_select(tmp_path, run_script, "--max-other", 2, "--min-markers", 1)
    assert result.returncode == 0, result.stderr
    markers = load_db(tmp_path)["markers"]
    assert markers["TAX100_M3"]["ext"] == ["TAX200"]
    assert markers["TAX100_M3"]["len"] == 700
    assert "TAX200_M1" in markers

    result = run_select(tmp_path, run_script, "--max-other", 2, "--min-markers", 2)
    assert result.returncode == 0, result.stderr
    markers = load_db(tmp_path)["markers"]
    assert markers["TAX100_M3"]["ext"] == []
    assert not any(marker.startswith("TAX200") for marker in markers)


def write_taxonomy(directory):
    directory.mkdir()
    nodes = [("1", "1", "no rank"), ("2", "1", "superkingdom"), ("1224", "2", "phylum"),
             ("1236", "1224", "class"), ("543", "1236", "family"), ("561", "543", "genus"),
             ("562", "561", "species"), ("83333", "562", "strain"),
             ("2759", "1", "superkingdom"), ("5794", "2759", "phylum"), ("5820", "5794", "genus"),
             ("5833", "5820", "species"),
             ("10239", "1", "superkingdom"), ("10710", "10239", "species")]
    names = [("1", "root"), ("2", "Bacteria"), ("1224", "Pseudomonadota"), ("1236", "Gammaproteobacteria"),
             ("543", "Enterobacteriaceae"), ("561", "Escherichia"), ("562", "Escherichia coli"),
             ("83333", "Escherichia coli K-12"), ("2759", "Eukaryota"), ("5794", "Apicomplexa"),
             ("5820", "Plasmodium"), ("5833", "Plasmodium falciparum"), ("10239", "Viruses"),
             ("10710", "Lambdavirus lambda")]
    (directory / "nodes.dmp").write_text(
        "".join(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n" for taxid, parent, rank in nodes), encoding="utf-8")
    (directory / "names.dmp").write_text(
        "".join(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for taxid, name in names), encoding="utf-8")


def test_metaphlan_lineages_name_missing_ranks_and_skip_viruses(tmp_path, run_script):
    """One row per taxid in MetaPhlAn's form: k__ for the domain, a missing rank
    named after the nearest one above with an empty taxid, t__TAX<taxid> as the
    bin; viruses are left out with a note."""
    write_taxonomy(tmp_path / "taxonomy")
    genomes = write(tmp_path / "genomes.tsv", "ecoli.fna\t83333\t0\npfal.fna\t5833\t1\nlambda.fna\t10710\t0\n")

    result = run_script("db_lineage.py", "--taxonomy", tmp_path / "taxonomy", "--genomes", genomes,
                        "--output", tmp_path / "sylph.tsv", "--mpa-output", tmp_path / "mpa.tsv")

    assert result.returncode == 0, result.stderr
    assert "lambda.fna" in result.stderr
    rows = {row["taxid"]: row for row in read_tsv(tmp_path / "mpa.tsv")}
    assert set(rows) == {"83333", "5833"}
    assert rows["83333"]["lineage"] == (
        "k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Gammaproteobacteria_unclassified|"
        "f__Enterobacteriaceae|g__Escherichia|s__Escherichia_coli|t__TAX83333")
    assert rows["83333"]["taxid_lineage"] == "2|1224|1236||543|561|562|"
    assert rows["5833"]["lineage"] == (
        "k__Eukaryota|p__Apicomplexa|c__Apicomplexa_unclassified|o__Apicomplexa_unclassified|"
        "f__Apicomplexa_unclassified|g__Plasmodium|s__Plasmodium_falciparum|t__TAX5833")
    assert rows["5833"]["domain"] == "Eukaryota"


def test_eukaryote_without_cds_stops_the_metaphlan_build(tmp_path, run_script):
    """A eukaryote with no coding sequences is refused before any gene is called."""
    write_taxonomy(tmp_path / "taxonomy")
    genomes = write(tmp_path / "genomes.tsv", "pfal.fna\t5833\t0\n")

    result = run_script("db_lineage.py", "--taxonomy", tmp_path / "taxonomy", "--genomes", genomes,
                        "--output", tmp_path / "sylph.tsv", "--mpa-output", tmp_path / "mpa.tsv")

    assert result.returncode != 0
    assert "pfal.fna" in result.stderr
    assert not (tmp_path / "mpa.tsv").exists()

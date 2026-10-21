"""bin/db_prepare_genome.py and bin/db_lineage.py: genomes made ready for a database build."""

import gzip


def test_headers_carry_the_taxid_and_a_map_is_written(tmp_path, run_script):
    """Every header becomes `<id>__<seqid>|kraken:taxid|<taxid>` (the description
    kept), and `<id>.fna.map` maps that id to the taxid, from gzipped input."""
    fasta = tmp_path / "genome.fa.gz"
    with gzip.open(fasta, "wt", encoding="utf-8") as handle:
        handle.write(">chr1 first contig\nACGT\nACGT\n>plasmid\nGGCC\n")

    result = run_script("db_prepare_genome.py", "--fasta", fasta, "--id", "g1", "--taxid", "562", cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    assert (tmp_path / "g1.fna").read_text(encoding="utf-8").splitlines() == [
        ">g1__chr1|kraken:taxid|562 first contig", "ACGT", "ACGT", ">g1__plasmid|kraken:taxid|562", "GGCC",
    ]
    assert (tmp_path / "g1.fna.map").read_text(encoding="utf-8").splitlines() == [
        "g1__chr1|kraken:taxid|562\t562", "g1__plasmid|kraken:taxid|562\t562",
    ]


def _taxonomy(directory):
    directory.mkdir()
    nodes = [("1", "1", "no rank"), ("2", "1", "domain"), ("1224", "2", "phylum"), ("1236", "1224", "class"),
             ("91347", "1236", "order"), ("543", "91347", "family"), ("561", "543", "genus"),
             ("562", "561", "species")]
    names = [("1", "root"), ("2", "Bacteria"), ("1224", "Pseudomonadota"), ("1236", "Gammaproteobacteria"),
             ("91347", "Enterobacterales"), ("543", "Enterobacteriaceae"), ("561", "Escherichia"),
             ("562", "Escherichia coli")]
    (directory / "nodes.dmp").write_text(
        "".join(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n" for taxid, parent, rank in nodes), encoding="utf-8")
    (directory / "names.dmp").write_text(
        "".join(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n" for taxid, name in names), encoding="utf-8")
    (directory / "merged.dmp").write_text("999\t|\t562\t|\n", encoding="utf-8")


def test_lineages_follow_merged_taxids(tmp_path, run_script):
    """Each genome file gets a GTDB-style NCBI lineage; a merged taxid resolves
    to the taxon it was merged into."""
    _taxonomy(tmp_path / "taxonomy")
    genomes = tmp_path / "genomes.tsv"
    genomes.write_text("g1.fna\t562\ng2.fna\t999\n", encoding="utf-8")
    output = tmp_path / "lineages.tsv"

    result = run_script("db_lineage.py", "--taxonomy", tmp_path / "taxonomy", "--genomes", genomes, "--output", output)
    assert result.returncode == 0, result.stderr

    lineage = ("d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;o__Enterobacterales;"
               "f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli")
    assert output.read_text(encoding="utf-8").splitlines() == [f"g1.fna\t{lineage}", f"g2.fna\t{lineage}"]


def test_unknown_taxid_stops_the_build(tmp_path, run_script):
    """A taxid the taxonomy does not know fails the step, naming it, and nothing is written."""
    _taxonomy(tmp_path / "taxonomy")
    genomes = tmp_path / "genomes.tsv"
    genomes.write_text("g1.fna\t562\ng3.fna\t12345\n", encoding="utf-8")
    output = tmp_path / "lineages.tsv"

    result = run_script("db_lineage.py", "--taxonomy", tmp_path / "taxonomy", "--genomes", genomes, "--output", output)
    assert result.returncode != 0
    assert "12345" in result.stderr
    assert not output.exists()

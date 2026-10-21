"""bin/sc_pseudobulk.py: per-library sums set side by side for the donor-level test."""

import csv


def _table(path, sample, groups):
    lines = ["sample\ttaxid\ttaxon\tcell_type\tstate\tn_cells\tgene\tcount"]
    for taxid, cell_type, state, gene, count in groups:
        lines.append(f"{sample}\t{taxid}\tTaxon {taxid}\t{cell_type}\t{state}\t12\t{gene}\t{count}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_groups_with_enough_paired_donors_are_written(tmp_path, run_script):
    """A taxon x cell type with both states in two donors becomes one count matrix
    (columns <donor>.<state>); one present in a single donor is listed as skipped."""
    both = [(100, "T cell", state, gene, 20) for state in ("infected", "bystander") for gene in ("G1", "G2")]
    _table(tmp_path / "D1.sc_host_de_pseudobulk.tsv", "D1", both + [(200, "T cell", "infected", "G1", 50)])
    _table(tmp_path / "D2.sc_host_de_pseudobulk.tsv", "D2", both)

    result = run_script("sc_pseudobulk.py", "--tables", tmp_path / "D1.sc_host_de_pseudobulk.tsv",
                        tmp_path / "D2.sc_host_de_pseudobulk.tsv", "--min-donors", "2", "--prefix", "t",
                        cwd=tmp_path)
    assert result.returncode == 0, result.stderr

    counts = (tmp_path / "t.100.T_cell.pseudobulk_counts.tsv").read_text(encoding="utf-8").splitlines()
    assert counts[0].split("\t") == ["gene", "D1.infected", "D1.bystander", "D2.infected", "D2.bystander"]
    assert counts[1].split("\t") == ["G1", "20", "20", "20", "20"]
    assert not list(tmp_path.glob("t.200.*.pseudobulk_counts.tsv"))

    with open(tmp_path / "t.pseudobulk_groups.tsv", encoding="utf-8") as handle:
        status = {row["taxid"]: row["status"] for row in csv.DictReader(handle, delimiter="\t")}
    assert status["100"] == "tested"
    assert status["200"].startswith("skipped")

    with open(tmp_path / "t.pseudobulk_metadata.tsv", encoding="utf-8") as handle:
        metadata = {row["column"]: (row["donor"], row["state"]) for row in csv.DictReader(handle, delimiter="\t")}
    assert metadata["D2.bystander"] == ("D2", "bystander")

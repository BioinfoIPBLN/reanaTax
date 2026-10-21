"""bin/filter_abundance.py: drop-list application and renormalisation of a combined Bracken table."""

import csv

import pytest

HEADER = ["name", "taxonomy_id", "taxonomy_lvl", "S1_num", "S1_frac", "S2_num", "S2_frac"]
ROWS = [
    ["Homo sapiens", "9606", "S", "500", "0.5", "100", "0.1"],
    ["Alpha bacterium", "1001", "S", "300", "0.3", "450", "0.45"],
    ["Beta bacterium", "1002", "S", "100", "0.1", "300", "0.3"],
    ["Gamma bacterium", "1003", "S", "100", "0.1", "150", "0.15"],
]


def test_drop_lists_remove_taxa_and_renormalise_fractions(tmp_path, run_script):
    """Taxa named in any drop list (taken as a union) leave the table, and every
    `_frac` column is rescaled to sum to 1 over the surviving rows, so abundances
    are fractions of what is left; `_num` counts are untouched."""
    table = tmp_path / "bracken_combined.tsv"
    table.write_text("\n".join("\t".join(row) for row in [HEADER] + ROWS) + "\n", encoding="utf-8")
    host_list = tmp_path / "host.txt"
    host_list.write_text("# host taxon\n9606\n", encoding="utf-8")
    minimizer_list = tmp_path / "minimizer_drop.txt"
    minimizer_list.write_text("1003\n", encoding="utf-8")

    result = run_script(
        "filter_abundance.py",
        "--input", table,
        "--output", tmp_path / "filtered.tsv",
        "--removed", tmp_path / "removed.tsv",
        "--drop-taxids", host_list, minimizer_list,
        "--min-rel-abundance", 0.001,
    )

    assert result.returncode == 0, result.stderr
    with open(tmp_path / "filtered.tsv", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == HEADER
        kept = {row["taxonomy_id"]: row for row in reader}
    assert set(kept) == {"1001", "1002"}
    assert float(kept["1001"]["S1_frac"]) == pytest.approx(0.75)
    assert float(kept["1002"]["S1_frac"]) == pytest.approx(0.25)
    assert float(kept["1001"]["S2_frac"]) == pytest.approx(0.6)
    assert float(kept["1002"]["S2_frac"]) == pytest.approx(0.4)
    for column in ("S1_frac", "S2_frac"):
        assert sum(float(row[column]) for row in kept.values()) == pytest.approx(1.0)
    assert (kept["1001"]["S1_num"], kept["1001"]["S2_num"]) == ("300", "450")
    assert (kept["1002"]["S1_num"], kept["1002"]["S2_num"]) == ("100", "300")

    with open(tmp_path / "removed.tsv", encoding="utf-8", newline="") as handle:
        removed = {row["taxonomy_id"] for row in csv.DictReader(handle, delimiter="\t")}
    assert removed == {"9606", "1003"}


def test_missing_drop_list_is_an_error(tmp_path, run_script):
    """A drop list that is not there stops the run instead of filtering nothing."""
    table = tmp_path / "bracken_combined.tsv"
    table.write_text("\n".join("\t".join(row) for row in [HEADER] + ROWS) + "\n", encoding="utf-8")

    result = run_script(
        "filter_abundance.py",
        "--input", table,
        "--output", tmp_path / "filtered.tsv",
        "--drop-taxids", tmp_path / "no_such_list.txt",
    )

    assert result.returncode != 0
    assert "no_such_list.txt" in result.stderr
    assert not (tmp_path / "filtered.tsv").exists()


def test_fractions_keep_one_denominator_across_passes(tmp_path, run_script):
    """The thresholds remove rows without rescaling, and a second pass on the
    filtered table (as the PRISM and gene-diversity filters run) rescales against
    the first pass's denominator rather than re-summing what the thresholds left:
    a kept taxon's fraction is its reads over the sample's reads less the
    dropped ones, whichever pass dropped them."""
    header = ["name", "taxonomy_id", "taxonomy_lvl", "S1_num", "S1_frac"]
    rows = [
        ["Alpha bacterium", "1001", "S", "6000", "0.6"],
        ["Beta bacterium", "1002", "S", "3000", "0.3"],
        ["Sparse bacterium", "1003", "S", "5", "0.0005"],
        ["Gamma bacterium", "1004", "S", "995", "0.0995"],
    ]
    table = tmp_path / "bracken_combined.tsv"
    table.write_text("\n".join("\t".join(row) for row in [header] + rows) + "\n", encoding="utf-8")

    first = run_script(
        "filter_abundance.py",
        "--input", table,
        "--output", tmp_path / "pass1.tsv",
        "--min-rel-abundance", 0.001,
    )
    assert first.returncode == 0, first.stderr
    with open(tmp_path / "pass1.tsv", encoding="utf-8", newline="") as handle:
        pass1 = {row["taxonomy_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert set(pass1) == {"1001", "1002", "1004"}
    assert pass1["1001"]["S1_frac"] == "0.6"

    drop = tmp_path / "prism_drop.txt"
    drop.write_text("1004\n", encoding="utf-8")
    second = run_script(
        "filter_abundance.py",
        "--input", tmp_path / "pass1.tsv",
        "--output", tmp_path / "pass2.tsv",
        "--drop-taxids", drop,
        "--min-rel-abundance", 0,
    )
    assert second.returncode == 0, second.stderr
    with open(tmp_path / "pass2.tsv", encoding="utf-8", newline="") as handle:
        pass2 = {row["taxonomy_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    # 10,000 reads, 995 dropped: the sparse taxon's 5 stay in the denominator.
    assert float(pass2["1001"]["S1_frac"]) == pytest.approx(6000 / 9005)
    assert float(pass2["1002"]["S1_frac"]) == pytest.approx(3000 / 9005)

"""bin/control_filter.py: per-sample comparison against negative-control libraries."""

import csv

import pytest

SAMPLES = ["S1", "S2", "NC1", "NC2"]

# Real libraries hold 1,000,000 classified reads, the blanks 10,000.
TAXA = [
    # taxid, name, {sample: reads}
    ("100", "Realbacter verus", {"S1": 500000, "S2": 400000, "NC1": 1, "NC2": 0}),
    # Same CPM (5000) in the samples as in the blanks.
    ("200", "Contaminobacter kitensis", {"S1": 5000, "S2": 5000, "NC1": 50, "NC2": 50}),
    # The reagent background that dominates the blanks.
    ("300", "Reagentella ubiquita", {"S1": 445000, "S2": 594950, "NC1": 9948, "NC2": 9950}),
    # Real in S1, carry-over level in S2.
    ("400", "Fusobacterium hoppens", {"S1": 50000, "S2": 50, "NC1": 1, "NC2": 0}),
]


def write_bracken_table(path, samples=SAMPLES):
    totals = {sample: sum(counts[sample] for _taxid, _name, counts in TAXA) for sample in samples}
    header = ["name", "taxonomy_id", "taxonomy_lvl"]
    for sample in samples:
        header += [f"{sample}_num", f"{sample}_frac"]
    lines = ["\t".join(header)]
    for taxid, name, counts in TAXA:
        fields = [name, taxid, "S"]
        for sample in samples:
            fields += [str(counts[sample]), f"{counts[sample] / totals[sample]:.6f}"]
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return totals


def read_tsv(path, key):
    with open(path, encoding="utf-8", newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def test_taxa_are_judged_per_sample_against_the_controls(tmp_path, run_script):
    """A taxon far above the blanks (in CPM) is kept; one at the blanks' level is
    zeroed - counts and fractions together - and listed for removal when that
    happens in every real sample; a taxon real in one sample and at control level
    in another is zeroed only where it fails."""
    totals = write_bracken_table(tmp_path / "bracken.tsv")
    assert totals == {"S1": 1000000, "S2": 1000000, "NC1": 10000, "NC2": 10000}

    result = run_script(
        "control_filter.py",
        "--input", tmp_path / "bracken.tsv",
        "--output", tmp_path / "filtered.tsv",
        "--evidence", tmp_path / "evidence.tsv",
        "--drop-list", tmp_path / "drop.txt",
        "--drop-cells", tmp_path / "cells.tsv",
        "--controls", "NC1,NC2",
        "--ratio", 10,
        "--statistic", "mean",
        "--floor-reads", 1,
        "--min-reads", 2,
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr

    evidence = read_tsv(tmp_path / "evidence.tsv", "taxid")
    assert evidence["100"]["verdict"] == "kept"
    assert evidence["200"]["verdict"] == "control_level_everywhere"
    assert evidence["300"]["verdict"] == "control_level_everywhere"
    assert evidence["400"]["verdict"] == "control_level_somewhere"
    assert (evidence["400"]["samples_kept"], evidence["400"]["samples_zeroed"]) == ("1", "1")

    filtered = read_tsv(tmp_path / "filtered.tsv", "taxonomy_id")
    assert filtered["100"]["S1_num"] == "500000"
    assert filtered["100"]["S2_num"] == "400000"
    for sample in ("S1", "S2"):
        assert filtered["200"][f"{sample}_num"] == "0"
        assert filtered["200"][f"{sample}_frac"] == "0"
    # Controls themselves are never zeroed.
    assert filtered["200"]["NC1_num"] == "50"
    assert filtered["400"]["S1_num"] == "50000"
    assert filtered["400"]["S2_num"] == "0"
    assert filtered["400"]["S2_frac"] == "0"
    # The zeroed reads leave the library's denominator: S1 lost 450,000 of its
    # million, S2 600,000.
    assert float(filtered["100"]["S1_frac"]) == pytest.approx(500000 / 550000)
    assert float(filtered["400"]["S1_frac"]) == pytest.approx(50000 / 550000)
    assert float(filtered["100"]["S2_frac"]) == pytest.approx(1.0)
    # A control is never zeroed, so its fractions keep their denominator.
    assert float(filtered["300"]["NC1_frac"]) == pytest.approx(0.9948)

    assert sorted((tmp_path / "drop.txt").read_text(encoding="utf-8").split()) == ["200", "300"]

    with open(tmp_path / "cells.tsv", encoding="utf-8", newline="") as handle:
        cells = {(row["taxid"], row["sample"]) for row in csv.DictReader(handle, delimiter="\t")}
    assert cells == {("200", "S1"), ("200", "S2"), ("300", "S1"), ("300", "S2"), ("400", "S2")}


def run_filter(run_script, tmp_path, *extra):
    return run_script(
        "control_filter.py",
        "--input", tmp_path / "bracken.tsv",
        "--output", tmp_path / "filtered.tsv",
        "--evidence", tmp_path / "evidence.tsv",
        "--drop-list", tmp_path / "drop.txt",
        "--drop-cells", tmp_path / "cells.tsv",
        "--ratio", 10,
        "--statistic", "mean",
        "--floor-reads", 1,
        "--min-reads", 2,
        *extra,
        cwd=tmp_path,
    )


def test_a_control_without_a_column_stops_the_filter(tmp_path, run_script):
    """A control named in --controls but missing from the table is an error."""
    write_bracken_table(tmp_path / "bracken.tsv", samples=["S1", "S2", "NC1"])

    result = run_filter(run_script, tmp_path, "--controls", "NC1,NC2")

    assert result.returncode != 0
    assert "no column for control sample(s) NC2" in result.stderr


def test_controls_left_out_upstream_are_reported_not_scored(tmp_path, run_script):
    """A control in --absent-controls is reported; the table is judged against the controls it holds."""
    write_bracken_table(tmp_path / "bracken.tsv", samples=["S1", "S2", "NC1"])

    result = run_filter(run_script, tmp_path, "--controls", "NC1", "--absent-controls", "NC2")

    assert result.returncode == 0, result.stderr
    assert "left out upstream: NC2" in result.stderr
    evidence = read_tsv(tmp_path / "evidence.tsv", "taxid")
    assert {row["controls"] for row in evidence.values()} == {"1"}
    assert evidence["100"]["verdict"] == "kept"
    assert evidence["200"]["verdict"] == "control_level_everywhere"


def test_no_control_left_passes_the_table_through(tmp_path, run_script):
    """With every control left out upstream and no --prevalence, no cell is zeroed and nothing is listed."""
    write_bracken_table(tmp_path / "bracken.tsv", samples=["S1", "S2"])

    result = run_filter(run_script, tmp_path, "--absent-controls", "NC1,NC2")

    assert result.returncode == 0, result.stderr
    assert "no control left to judge against" in result.stderr
    assert (tmp_path / "filtered.tsv").read_text(encoding="utf-8") == (tmp_path / "bracken.tsv").read_text(encoding="utf-8")
    assert (tmp_path / "drop.txt").read_text(encoding="utf-8") == ""
    assert (tmp_path / "cells.tsv").read_text(encoding="utf-8") == "taxid\tsample\n"
    evidence = read_tsv(tmp_path / "evidence.tsv", "taxid")
    assert {row["verdict"] for row in evidence.values()} == {"kept"}

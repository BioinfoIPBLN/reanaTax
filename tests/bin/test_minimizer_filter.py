"""bin/minimizer_filter.py: judge taxa on the breadth of their minimizer evidence."""

import csv

# Kraken2 --report-minimizer-data layout:
# pct, clade_reads, direct_reads, minimizers, distinct, rank, taxid, name
REPORT_ROWS = [
    ("1.00", 1000, 1000, 0, 0, "U", 0, "unclassified"),
    ("99.00", 101710, 0, 1590090, 7207, "R", 1, "root"),
    # Under the read gate: never judged, although it would fail few_distinct.
    ("0.01", 10, 10, 50, 2, "S", 1001, "Below gate"),
    # Plenty of reads on five distinct minimizers.
    ("0.20", 200, 200, 40, 5, "S", 1002, "Few distinct"),
    # Duplication 250 over the flat limit, and too few distinct (200 < 1000)
    # for the limit to rise with its 20% coverage.
    ("0.49", 500, 500, 50000, 200, "S", 1003, "Single locus"),
    # Duplication 300 over the flat limit, but 5000 distinct at 50% coverage:
    # the limit rises to 100 x 0.5 / 0.01 = 5000.
    ("97.34", 100000, 100000, 1500000, 5000, "S", 1004, "Small genome deep"),
    ("0.97", 1000, 1000, 20000, 2000, "S", 1005, "Clean"),
]

# kraken2-inspect: pct, clade minimizers, taxon minimizers, rank, taxid, name
INSPECT_ROWS = [
    ("100.00", 2000000, 0, "R", 1, "root"),
    ("0.05", 1000, 1000, "S", 1001, "Below gate"),
    ("5.00", 100000, 100000, "S", 1002, "Few distinct"),
    ("0.05", 1000, 1000, "S", 1003, "Single locus"),
    ("0.50", 10000, 10000, "S", 1004, "Small genome deep"),
    ("5.00", 100000, 100000, "S", 1005, "Clean"),
]


def write_inputs(tmp_path):
    report = tmp_path / "S1.kraken2.report.txt"
    report.write_text(
        "".join(
            f"{pct}\t{clade}\t{direct}\t{minimizers}\t{distinct}\t{rank}\t{taxid}\t    {name}\n"
            for pct, clade, direct, minimizers, distinct, rank, taxid, name in REPORT_ROWS
        ),
        encoding="utf-8",
    )
    inspect = tmp_path / "inspect.txt"
    inspect.write_text(
        "# Database options: nucleotide db, k = 35, l = 31\n"
        "# Spaced mask = 11111111111111111111111111111111110011001100110011001100110011\n"
        "# Toggle mask = 1110001101111110001010001100010000100111000110110101101000101101\n"
        "# Total taxonomy nodes: 6\n"
        "# Table size: 2000000\n"
        "# Table capacity: 3000000\n"
        "# Min clear hash value = 0\n"
        + "".join(
            f"{pct}\t{clade}\t{direct}\t{rank}\t{taxid}\t  {name}\n"
            for pct, clade, direct, rank, taxid, name in INSPECT_ROWS
        ),
        encoding="utf-8",
    )
    return report, inspect


def run_filter(run_script, tmp_path, *extra):
    report, inspect = write_inputs(tmp_path)
    evidence = tmp_path / "evidence.tsv"
    drop = tmp_path / "drop.txt"
    result = run_script(
        "minimizer_filter.py",
        "--reports", report,
        "--inspect", inspect,
        "--evidence", evidence,
        "--drop-list", drop,
        "--min-reads", 50,
        "--min-distinct", 10,
        "--max-duplication", 100,
        "--duplication-min-distinct", 1000,
        *extra,
    )
    assert result.returncode == 0, result.stderr
    with open(evidence, encoding="utf-8", newline="") as handle:
        rows = {row["taxid"]: row for row in csv.DictReader(handle, delimiter="\t")}
    dropped = drop.read_text(encoding="utf-8").split()
    return rows, dropped


def test_verdicts_and_coverage_scaled_duplication_limit(tmp_path, run_script):
    """Few reads are left unjudged; few distinct minimizers or high duplication
    drop a taxon; but duplication over the flat limit is tolerated once breadth is
    genome-scale, because the limit scales with coverage (the RSV/rubella fix)."""
    rows, dropped = run_filter(run_script, tmp_path, "--duplication-coverage", 0.01)

    assert "0" not in rows and "1" not in rows
    assert rows["1001"]["verdict"] == "below_read_gate"
    assert rows["1002"]["verdict"] == "few_distinct"
    assert rows["1003"]["verdict"] == "high_duplication"
    assert float(rows["1003"]["duplication_limit"]) == 100.0
    assert rows["1004"]["verdict"] == "kept"
    assert float(rows["1004"]["duplication"]) == 300.0
    assert float(rows["1004"]["coverage"]) == 0.5
    assert float(rows["1004"]["duplication_limit"]) == 5000.0
    assert rows["1005"]["verdict"] == "kept"
    assert sorted(dropped) == ["1002", "1003"]


def test_flat_duplication_limit_condemns_the_genome_scale_taxon(tmp_path, run_script):
    """With the coverage scaling switched off the same deep small-genome taxon is
    dropped for duplication, so it is the scaling that keeps it above."""
    rows, dropped = run_filter(run_script, tmp_path, "--duplication-coverage", 0)

    assert rows["1004"]["verdict"] == "high_duplication"
    assert float(rows["1004"]["duplication_limit"]) == 100.0
    assert sorted(dropped) == ["1002", "1003", "1004"]


def test_depth_scaling_counts_each_classified_read_once(tmp_path, run_script):
    """--distinct-scale takes depth from the root's clade count. A species and the
    strain under it hold the same 6,000 reads; summing every S-rank row counted
    them twice and raised the threshold from 1,000 to 1,600 distinct minimizers,
    which would have dropped a species with 1,200."""
    report = tmp_path / "S1.kraken2.report.txt"
    rows = [
        ("4.76", 500, 500, 0, 0, "U", 0, "unclassified"),
        ("95.24", 10000, 0, 70000, 4200, "R", 1, "root"),
        ("57.14", 6000, 0, 50000, 1200, "S", 1005, "Clean"),
        ("57.14", 6000, 6000, 50000, 1200, "S1", 2001, "Clean strain"),
        ("38.10", 4000, 4000, 20000, 3000, "S", 1006, "Other"),
    ]
    report.write_text(
        "".join(
            f"{pct}\t{clade}\t{direct}\t{minimizers}\t{distinct}\t{rank}\t{taxid}\t    {name}\n"
            for pct, clade, direct, minimizers, distinct, rank, taxid, name in rows
        ),
        encoding="utf-8",
    )
    evidence = tmp_path / "evidence.tsv"
    result = run_script(
        "minimizer_filter.py",
        "--reports", report,
        "--evidence", evidence,
        "--drop-list", tmp_path / "drop.txt",
        "--min-reads", 50,
        "--min-distinct", 10,
        "--distinct-scale", 0.1,
    )

    assert result.returncode == 0, result.stderr
    assert "mean 10000 classified reads x 0.1 -> distinct threshold 1000" in result.stderr
    with open(evidence, encoding="utf-8", newline="") as handle:
        verdicts = {row["taxid"]: row["verdict"] for row in csv.DictReader(handle, delimiter="\t")}
    assert verdicts["1005"] == "kept"
    assert verdicts["1006"] == "kept"

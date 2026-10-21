"""bin/rpm_total.py: counts restated as reads per million of each whole library."""

import csv


def test_counts_become_reads_per_million_of_the_library(tmp_path, run_script):
    """Each `<s>_num` count is divided by the sample's trimmed reads from the read
    accounting and multiplied by a million; name, taxid and rank are carried, the
    fractions are not, and a sample with no library size is NA, not a guess."""
    table = tmp_path / "bracken_combined_S.filtered.tsv"
    table.write_text(
        "name\ttaxonomy_id\ttaxonomy_lvl\tS1_num\tS1_frac\tS2_num\tS2_frac\n"
        "Alpha bacterium\t1001\tS\t500\t0.5\t30\t0.3\n"
        "Beta bacterium\t1002\tS\t500\t0.5\t70\t0.7\n",
        encoding="utf-8",
    )
    accounting = tmp_path / "reanatax.read_accounting.tsv"
    accounting.write_text("sample\traw_reads\ttrimmed_reads\nS1\t2500000\t2000000\n", encoding="utf-8")
    output = tmp_path / "rpm.tsv"

    result = run_script("rpm_total.py", "--table", table, "--accounting", accounting, "--output", output)
    assert result.returncode == 0, result.stderr

    with open(output, encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0].keys()) == ["name", "taxonomy_id", "taxonomy_lvl", "S1_rpm", "S2_rpm"]
    assert rows[0]["S1_rpm"] == "250"
    assert rows[1]["S1_rpm"] == "250"
    assert rows[0]["S2_rpm"] == "NA"

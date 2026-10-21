"""bin/filter_abundance.py: a separate relative-abundance floor for viruses."""

HEADER = ["name", "taxonomy_id", "taxonomy_lvl", "S1_num", "S1_frac"]
ROWS = [
    ["Alpha bacterium", "1001", "S", "9995", "0.9995"],
    ["Small virus", "2001", "S", "5", "0.0005"],
]


def _filter(tmp_path, run_script, *extra):
    table = tmp_path / "bracken_combined.tsv"
    table.write_text("\n".join("\t".join(row) for row in [HEADER] + ROWS) + "\n", encoding="utf-8")
    output = tmp_path / "filtered.tsv"
    result = run_script("filter_abundance.py", "--input", table, "--output", output,
                        "--min-rel-abundance", "0.001", *extra)
    assert result.returncode == 0, result.stderr
    return [line.split("\t")[1] for line in output.read_text(encoding="utf-8").splitlines()[1:]]


def test_one_floor_removes_the_virus(tmp_path, run_script):
    """At 0.05% the virus is under the 0.1% floor that applies to everything."""
    assert _filter(tmp_path, run_script) == ["1001"]


def test_virus_floor_zero_keeps_it(tmp_path, run_script):
    """Named in --virus-taxids, with --virus-min-rel-abundance 0, the virus is
    judged against its own floor and kept; the bacterium is judged as before."""
    viruses = tmp_path / "viruses.txt"
    viruses.write_text("10239\n2001\n", encoding="utf-8")
    kept = _filter(tmp_path, run_script, "--virus-taxids", viruses, "--virus-min-rel-abundance", "0")
    assert kept == ["1001", "2001"]

"""bin/bracken_unique_names.py: one name per taxid before combine_bracken_outputs.py."""

HEADER = "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n"


def test_a_name_shared_by_two_taxids_is_labelled_in_every_file(tmp_path, run_script):
    """SILVA's "uncultured" genera: one sample carries two of them, the other a
    third. All three rows get their taxid in the name, in both files; a name
    with a single taxid is copied unchanged."""
    a = tmp_path / "A.bracken_G.tsv"
    a.write_text(
        HEADER
        + "Wolbachia\t2770\tG\t6745\t3261\t10006\t0.10468\n"
        + "uncultured\t45238\tG\t13\t142\t155\t0.00162\n"
        + "uncultured\t45557\tG\t22\t386\t408\t0.00427\n",
        encoding="utf-8",
    )
    b = tmp_path / "B.bracken_G.tsv"
    b.write_text(
        HEADER
        + "Wolbachia\t2770\tG\t10\t0\t10\t0.5\n"
        + "uncultured\t45010\tG\t10\t0\t10\t0.5\n",
        encoding="utf-8",
    )
    out = tmp_path / "unique"

    result = run_script("bracken_unique_names.py", "--outdir", out, a, b)
    assert result.returncode == 0, result.stderr
    assert "'uncultured' names 3 taxa (45010, 45238, 45557)" in result.stderr

    a_rows = [line.split("\t") for line in (out / a.name).read_text(encoding="utf-8").splitlines()]
    assert a_rows[0] == HEADER.rstrip("\n").split("\t")
    assert [row[0] for row in a_rows[1:]] == ["Wolbachia", "uncultured (taxid 45238)", "uncultured (taxid 45557)"]
    assert a_rows[2][1:] == ["45238", "G", "13", "142", "155", "0.00162"]

    b_names = [line.split("\t")[0] for line in (out / b.name).read_text(encoding="utf-8").splitlines()[1:]]
    assert b_names == ["Wolbachia", "uncultured (taxid 45010)"]


def test_a_cohort_without_shared_names_is_copied_byte_for_byte(tmp_path, run_script):
    text = HEADER + "Escherichia coli\t562\tS\t100\t5\t105\t1.00000\n"
    a = tmp_path / "A.bracken_S.tsv"
    a.write_text(text, encoding="utf-8")
    out = tmp_path / "unique"

    result = run_script("bracken_unique_names.py", "--outdir", out, a)
    assert result.returncode == 0, result.stderr
    assert result.stderr == ""
    assert (out / a.name).read_text(encoding="utf-8") == text


def test_a_file_that_is_not_bracken_output_is_refused(tmp_path, run_script):
    report = tmp_path / "A.kreport.txt"
    report.write_text("100.00\t10\t0\tR\t1\troot\n", encoding="utf-8")

    result = run_script("bracken_unique_names.py", "--outdir", tmp_path / "unique", report)
    assert result.returncode != 0
    assert "not a Bracken output" in result.stderr

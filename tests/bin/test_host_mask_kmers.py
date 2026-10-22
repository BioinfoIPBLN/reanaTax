"""bin/host_mask_kmers.py: the protected genomes cut into k-mers for --host_mask."""

import gzip


def _records(text):
    lines = text.splitlines()
    return list(zip(lines[0::2], lines[1::2]))


def test_kmers_from_plain_and_gzipped_targets(tmp_path, run_script):
    """Every k-mer at each step, named <sequence>_<offset>, from every file; a
    sequence shorter than k gives none."""
    plain = tmp_path / "a.fa"
    plain.write_text(">seqA description\nACGTAC\nGTAC\n>tiny\nACG\n", encoding="utf-8")
    zipped = tmp_path / "b.fa.gz"
    with gzip.open(zipped, "wt", encoding="utf-8") as handle:
        handle.write(">seqB\nTTTTTGGGGG\n")
    out = tmp_path / "kmers.fa"

    result = run_script(
        "host_mask_kmers.py",
        "--targets", plain, zipped,
        "--length", 5,
        "--step", 2,
        "--output", out,
    )

    assert result.returncode == 0, result.stderr
    assert _records(out.read_text(encoding="utf-8")) == [
        (">seqA_0", "ACGTA"),
        (">seqA_2", "GTACG"),
        (">seqA_4", "ACGTA"),
        (">seqB_0", "TTTTT"),
        (">seqB_2", "TTTGG"),
        (">seqB_4", "TGGGG"),
    ]


def test_kmers_to_standard_output(tmp_path, run_script):
    """--output - streams the k-mers, as HOST_MASK pipes them into minimap2."""
    target = tmp_path / "a.fa"
    target.write_text(">s\nACGTACGT\n", encoding="utf-8")

    result = run_script("host_mask_kmers.py", "--targets", target, "--length", 8, "--step", 1, "--output", "-")

    assert result.returncode == 0, result.stderr
    assert _records(result.stdout) == [(">s_0", "ACGTACGT")]


def test_no_kmer_at_all_is_an_error(tmp_path, run_script):
    """Targets with no sequence as long as k would mask nothing; say so."""
    target = tmp_path / "a.fa"
    target.write_text(">s\nACGT\n", encoding="utf-8")

    result = run_script("host_mask_kmers.py", "--targets", target, "--output", tmp_path / "k.fa")

    assert result.returncode != 0
    assert "no 150-mer" in result.stderr

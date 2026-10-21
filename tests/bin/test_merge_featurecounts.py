"""bin/merge_featurecounts.py: per-sample featureCounts tables -> one gene x sample matrix."""

GENES = [
    ("GENE1", "chr1", "100", "200", "+", 101),
    ("GENE2", "chr1;chr1", "300;500", "400;600", "+;+", 202),
    ("GENE3", "chr2", "1000", "1999", "-", 1000),
]


def write_featurecounts(path, bam_name, counts):
    lines = [
        '# Program:featureCounts v2.0.6; Command:"featureCounts" "-a" "genes.gtf" "-o" "out" "' + bam_name + '"',
        "\t".join(["Geneid", "Chr", "Start", "End", "Strand", "Length", bam_name]),
    ]
    for (gene, chrom, start, end, strand, length), count in zip(GENES, counts):
        lines.append("\t".join([gene, chrom, start, end, strand, str(length), str(count)]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_matrix(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        fields = line.split("\t")
        rows[fields[0]] = dict(zip(header[1:], fields[1:]))
    return header, rows


def test_depletion_and_target_labels_are_stripped_from_sample_names(tmp_path, run_script):
    """The `.host`, `.hostN` and `.target` pass labels are not part of the sample id:
    left on, the matrix columns never match --da_metadata."""
    write_featurecounts(tmp_path / "S1.host.featureCounts.tsv", "S1.host.sorted.bam", [5, 0, 12])
    write_featurecounts(tmp_path / "S2.host1.featureCounts.tsv", "S2.host1.sorted.bam", [7, 3, 0])
    write_featurecounts(tmp_path / "S3.target.featureCounts.tsv", "S3.target.sorted.bam", [1, 2, 3])
    out = tmp_path / "merged.tsv"
    lengths = tmp_path / "lengths.tsv"

    result = run_script(
        "merge_featurecounts.py",
        tmp_path / "S1.host.featureCounts.tsv",
        tmp_path / "S2.host1.featureCounts.tsv",
        tmp_path / "S3.target.featureCounts.tsv",
        "--output", out,
        "--lengths", lengths,
    )

    assert result.returncode == 0, result.stderr
    header, rows = read_matrix(out)
    assert header == ["gene_id", "S1", "S2", "S3"]
    assert rows["GENE1"] == {"S1": "5", "S2": "7", "S3": "1"}
    assert rows["GENE2"] == {"S1": "0", "S2": "3", "S3": "2"}
    assert rows["GENE3"] == {"S1": "12", "S2": "0", "S3": "3"}
    assert lengths.read_text(encoding="utf-8").splitlines() == [
        "gene_id\tlength", "GENE1\t101", "GENE2\t202", "GENE3\t1000",
    ]


def test_two_tables_for_one_sample_are_refused(tmp_path, run_script):
    """Two tables resolving to the same sample must fail loudly instead of one
    silently overwriting the other in the matrix."""
    write_featurecounts(tmp_path / "S1.host.featureCounts.tsv", "S1.host.sorted.bam", [5, 0, 12])
    write_featurecounts(tmp_path / "S1.host1.featureCounts.tsv", "S1.host1.sorted.bam", [6, 1, 11])

    result = run_script(
        "merge_featurecounts.py",
        tmp_path / "S1.host.featureCounts.tsv",
        tmp_path / "S1.host1.featureCounts.tsv",
        "--output", tmp_path / "merged.tsv",
    )

    assert result.returncode != 0
    assert "two tables map to sample 'S1'" in result.stderr
    assert not (tmp_path / "merged.tsv").exists()

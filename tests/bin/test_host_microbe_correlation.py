"""bin/host_microbe_correlation.py: host gene x taxon Spearman grid with a feasibility guard."""

import itertools

HOST_SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6", "S7"]
HOST = {
    "GENE_A": [100, 200, 300, 400, 500, 600, 700],
    "GENE_B": [600, 500, 400, 300, 200, 100, 50],
    "GENE_C": [50, 10, 80, 20, 90, 30, 60],
    "GENE_ZERO": [0, 0, 0, 0, 0, 0, 0],
}

# Fractions per sample S1..S6; S7 is host-only and must be ignored.
MICROBES = {
    "Alpha bacterium": ("11", [0.10, 0.15, 0.20, 0.25, 0.30, 0.35]),
    "Beta bacterium": ("12", [0.60, 0.50, 0.45, 0.40, 0.30, 0.25]),
    "Rare bacterium": ("13", [0.00, 0.00, 0.00, 0.00, 0.05, 0.05]),
    "Gamma bacterium": ("14", [0.30, 0.35, 0.35, 0.35, 0.35, 0.35]),
}


def write_host(path, samples, genes):
    lines = ["gene_id\t" + "\t".join(samples)]
    for gene, counts in genes.items():
        lines.append(gene + "\t" + "\t".join(str(count) for count in counts[: len(samples)]))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_bracken(path, samples, microbes):
    """Combined Bracken layout; columns deliberately in reverse sample order."""
    order = list(reversed(range(len(samples))))
    header = ["name", "taxonomy_id", "taxonomy_lvl"]
    for i in order:
        header += [f"{samples[i]}_num", f"{samples[i]}_frac"]
    lines = ["\t".join(header)]
    for name, (taxid, fractions) in microbes.items():
        fields = [name, taxid, "S"]
        for i in order:
            fields += [str(round(fractions[i] * 1000)), f"{fractions[i]:.4f}"]
        lines.append("\t".join(fields))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_results(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    comments = [line for line in lines if line.startswith("#")]
    table = [line.split("\t") for line in lines if not line.startswith("#")]
    header, rows = table[0], table[1:]
    return comments, header, [dict(zip(header, row)) for row in rows]


def test_one_row_per_tested_gene_taxon_pair(tmp_path, run_script):
    """Every gene and taxon that passes the non-zero/variability filter is tested
    against every other exactly once, samples are matched by NAME (not column
    position), and genes/taxa seen in too few samples are left out."""
    write_host(tmp_path / "host.tsv", HOST_SAMPLES, HOST)
    write_bracken(tmp_path / "bracken.tsv", HOST_SAMPLES[:6], MICROBES)
    out = tmp_path / "correlation.tsv"

    result = run_script(
        "host_microbe_correlation.py",
        "--host", tmp_path / "host.tsv",
        "--microbial", tmp_path / "bracken.tsv",
        "--output", out,
        "--mqc", tmp_path / "correlation_mqc.tsv",
    )

    assert result.returncode == 0, result.stderr
    comments, header, rows = read_results(out)
    assert header == ["gene_id", "taxon", "rho", "p_value", "q_value", "significant"]
    assert any("over 6 samples: S1, S2, S3, S4, S5, S6" in line for line in comments)

    pairs = [(row["gene_id"], row["taxon"]) for row in rows]
    expected = set(itertools.product(
        ["GENE_A", "GENE_B", "GENE_C"],
        ["Alpha bacterium", "Beta bacterium", "Gamma bacterium"],
    ))
    assert len(pairs) == len(set(pairs)) == 9
    assert set(pairs) == expected

    by_pair = {(row["gene_id"], row["taxon"]): row for row in rows}
    assert float(by_pair[("GENE_A", "Alpha bacterium")]["rho"]) == 1.0
    assert float(by_pair[("GENE_B", "Alpha bacterium")]["rho"]) == -1.0
    for row in rows:
        assert row["significant"] in ("yes", "no")
        assert 0.0 <= float(row["q_value"]) <= 1.0


def test_cohort_too_small_for_any_significant_pair_is_refused(tmp_path, run_script):
    """With four samples the smallest attainable Spearman p (2/4! = 0.083) cannot
    clear 0.05, so the run is refused instead of writing an empty-looking
    'no associations' table; --force overrides."""
    samples = ["S1", "S2", "S3", "S4"]
    write_host(tmp_path / "host.tsv", samples, {k: v for k, v in HOST.items() if k != "GENE_ZERO"})
    write_bracken(
        tmp_path / "bracken.tsv",
        samples,
        {name: (taxid, fractions[:4]) for name, (taxid, fractions) in MICROBES.items() if name.startswith(("Alpha", "Beta"))},
    )
    out = tmp_path / "correlation.tsv"
    args = ["--host", tmp_path / "host.tsv", "--microbial", tmp_path / "bracken.tsv", "--output", out]

    refused = run_script("host_microbe_correlation.py", *args)
    assert refused.returncode != 0
    assert "REFUSING" in refused.stderr
    assert not out.exists()

    forced = run_script("host_microbe_correlation.py", *args, "--force")
    assert forced.returncode == 0, forced.stderr
    _comments, _header, rows = read_results(out)
    assert len(rows) == 3 * 2


def test_combine_kreports_table_is_refused(tmp_path, run_script):
    """A hierarchical kreport table is refused by name instead of failing on
    its preamble; the correlation wants the flat Bracken table."""
    write_host(tmp_path / "host.tsv", HOST_SAMPLES, HOST)
    kreport = tmp_path / "kraken2_combined_report.txt"
    kreport.write_text(
        "#Number of Samples: 2\n#Total Number of Reads: 2000\n"
        "#perc\ttot_all\ttot_lvl\tS1_all\tS1_lvl\tS2_all\tS2_lvl\tlvl_type\ttaxid\tname\n"
        "100.00\t2000\t0\t1000\t0\t1000\t0\tR\t1\troot\n",
        encoding="utf-8",
    )

    result = run_script(
        "host_microbe_correlation.py",
        "--host", tmp_path / "host.tsv",
        "--microbial", kreport,
        "--output", tmp_path / "correlation.tsv",
    )

    assert result.returncode != 0
    assert "combine_kreports" in result.stderr
    assert not (tmp_path / "correlation.tsv").exists()

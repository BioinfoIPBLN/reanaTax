"""bin/runinfo_to_runsheet.py: fastq-dl run-info TSV -> slim, line-safe run sheet."""

import csv

# A realistic subset of ENA's `fields=all`, in fastq-dl's order, without
# `base_count`: the SRA fallback returns a smaller field set.
FIELDS = [
    "study_accession", "sample_accession", "experiment_accession", "run_accession",
    "scientific_name", "instrument_platform", "instrument_model", "library_layout",
    "library_strategy", "library_source", "read_count", "sample_title", "fastq_ftp",
    "description",
]


def run_record(run, **overrides):
    record = {
        "study_accession": "PRJNA000001",
        "sample_accession": "SAMN00000001",
        "experiment_accession": "SRX0000001",
        "run_accession": run,
        "scientific_name": "Homo sapiens",
        "instrument_platform": "ILLUMINA",
        "instrument_model": "Illumina NovaSeq 6000",
        "library_layout": "PAIRED",
        "library_strategy": "RNA-Seq",
        "library_source": "TRANSCRIPTOMIC",
        "read_count": "1000",
        "sample_title": "plain title",
        "fastq_ftp": f"ftp.sra.ebi.ac.uk/vol1/fastq/{run}_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/{run}_2.fastq.gz",
        "description": "",
    }
    record.update(overrides)
    return record


def write_runinfo(path, records):
    """Written the way fastq-dl writes it: csv.DictWriter, tab-separated, so
    free-text fields with tabs, quotes or newlines arrive quoted."""
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)


def test_run_sheet_is_one_clean_line_per_unique_run(tmp_path, run_script):
    """Whitespace inside free text is collapsed so every run is one CSV line,
    duplicate runs are kept once, a missing experiment/sample accession falls back
    to the run accession, absent columns come out empty, and the query is recorded."""
    write_runinfo(tmp_path / "PRJNA000001-run-info.tsv", [
        run_record("SRR1000001", sample_title="Patient 1\tbaseline\nvisit  two",
                   description='free "quoted"\ttext'),
        run_record("ERR2000002", experiment_accession="", sample_accession=""),
        run_record("SRR1000001"),
    ])
    out = tmp_path / "runsheet.csv"

    result = run_script(
        "runinfo_to_runsheet.py",
        tmp_path / "PRJNA000001-run-info.tsv",
        out,
        "--query", "PRJNA000001",
    )

    assert result.returncode == 0, result.stderr
    raw_lines = out.read_text(encoding="utf-8").splitlines()
    assert len(raw_lines) == 3
    with open(out, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        assert reader.fieldnames[:2] == ["query", "run_accession"]
        rows = list(reader)
    assert [row["run_accession"] for row in rows] == ["SRR1000001", "ERR2000002"]
    assert all(row["query"] == "PRJNA000001" for row in rows)
    assert rows[0]["sample_title"] == "Patient 1 baseline visit two"
    assert rows[0]["instrument_model"] == "Illumina NovaSeq 6000"
    assert rows[0]["experiment_accession"] == "SRX0000001"
    assert rows[1]["experiment_accession"] == "ERR2000002"
    assert rows[1]["sample_accession"] == "ERR2000002"
    assert all(row["base_count"] == "" for row in rows)
    assert "description" not in rows[0]


def test_ena_error_line_as_run_accession_fails(tmp_path, run_script):
    """ENA sometimes returns its error message as a data line; fastq-dl writes it
    as a run. Anything that is not an [SED]RR accession must exit non-zero (so
    the task is retried) rather than be skipped, and no run sheet is written."""
    write_runinfo(tmp_path / "PRJNA000001-run-info.tsv", [
        run_record("SRR1000001"),
        run_record("ERROR occurred. Not all results may have been written.Query: accession=PRJNA000001",
                   experiment_accession="", sample_accession=""),
    ])
    out = tmp_path / "runsheet.csv"

    result = run_script(
        "runinfo_to_runsheet.py",
        tmp_path / "PRJNA000001-run-info.tsv",
        out,
        "--query", "PRJNA000001",
    )

    assert result.returncode != 0
    assert "instead of a run accession" in result.stderr
    assert not out.exists()

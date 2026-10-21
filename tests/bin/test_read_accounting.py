"""bin/read_accounting.py: per-sample read fate across trimming, depletion and classification."""

import csv
import json

# HISAT2's default (Bowtie2-style) summary for a paired-end library.
HOST_LOG = """\
1000 reads; of these:
  1000 (100.00%) were paired; of these:
    600 (60.00%) aligned concordantly 0 times
    350 (35.00%) aligned concordantly exactly 1 time
    50 (5.00%) aligned concordantly >1 times
    ----
    600 pairs aligned concordantly 0 times; of these:
      20 (3.33%) aligned discordantly 1 time
    ----
    580 pairs aligned 0 times concordantly or discordantly; of these:
      1160 mates make up the pairs; of these:
        1100 (94.83%) aligned 0 times
        40 (3.45%) aligned exactly 1 time
        20 (1.72%) aligned >1 times
45.00% overall alignment rate
"""

UNIVEC_LOG = """\
600 reads; of these:
  600 (100.00%) were paired; of these:
    590 (98.33%) aligned concordantly 0 times
    8 (1.33%) aligned concordantly exactly 1 time
    2 (0.33%) aligned concordantly >1 times
    ----
    590 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    590 pairs aligned 0 times concordantly or discordantly; of these:
      1180 mates make up the pairs; of these:
        1178 (99.83%) aligned 0 times
        2 (0.17%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
1.83% overall alignment rate
"""


def test_dotted_sample_ids_keep_their_depletion_passes_together(tmp_path, run_script):
    """`<sample>.<label>.hisat2.summary.log` with a dotted sample id must stay one
    sample (`Pt01.Visit2`) carrying a `host` and a `univec` pass, not be split at
    the first dot into a sample `Pt01` with a pass called `Visit2`."""
    sample = "Pt01.Visit2"
    host = tmp_path / f"{sample}.host.hisat2.summary.log"
    univec = tmp_path / f"{sample}.univec.hisat2.summary.log"
    fastp = tmp_path / f"{sample}.fastp.json"
    host.write_text(HOST_LOG, encoding="utf-8")
    univec.write_text(UNIVEC_LOG, encoding="utf-8")
    fastp.write_text(json.dumps({
        "summary": {
            "before_filtering": {"total_reads": 2200},
            "after_filtering": {"total_reads": 2000, "read1_mean_length": 100, "read2_mean_length": 100},
        }
    }), encoding="utf-8")
    prefix = tmp_path / "run"

    result = run_script(
        "read_accounting.py",
        "--fastp", fastp,
        "--hisat2", host, univec,
        "--prefix", prefix,
        cwd=tmp_path,
    )

    assert result.returncode == 0, result.stderr
    with open(f"{prefix}.read_accounting.tsv", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["sample"] for row in rows] == [sample]
    row = rows[0]
    assert not any(column.startswith("Visit2_") for column in row)
    # Paired: HISAT2 counts pairs, the table counts reads.
    assert int(row["host_input_reads"]) == 2000
    assert int(row["host_aligned_reads"]) == 800
    assert int(row["host_unaligned_reads"]) == 1200
    assert int(row["univec_input_reads"]) == 1200
    assert int(row["univec_aligned_reads"]) == 20
    assert int(row["univec_unaligned_reads"]) == 1180
    # univec runs after the host pass, so it is the one whose survivors are non-host.
    assert int(row["host_removed_reads"]) == 820
    assert int(row["nonhost_reads"]) == 1180
    assert int(row["raw_reads"]) == 2200
    assert float(row["nonhost_pct_of_raw"]) == 53.64

    fate = (tmp_path / "run_read_fate_mqc.tsv").read_text(encoding="utf-8").splitlines()
    data = [line for line in fate if not line.startswith("#")]
    assert data[1].split("\t")[0] == sample
    assert len(data) == 2

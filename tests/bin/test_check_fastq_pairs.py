"""bin/check_fastq_pairs.sh: a paired samtools fastq run that lost its pairs fails."""

import subprocess
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"


def run_check(tmp_path, log_text):
    log = tmp_path / "x.fastq.log"
    log.write_text(log_text, encoding="utf-8")
    return subprocess.run(["bash", str(BIN / "check_fastq_pairs.sh"), str(log)],
                          capture_output=True, text=True, timeout=60)


def test_every_read_a_singleton_fails(tmp_path):
    """PRJNA1367057 SRX31239398: '@SRR36203002.1/1 1 length=100' headers, every
    read set aside, empty non-host FASTQs."""
    result = run_check(tmp_path, "[M::bam2fq_mainloop] discarded 8781870 singletons\n"
                                 "[M::bam2fq_mainloop] processed 8781870 reads\n")
    assert result.returncode == 1
    assert "set aside 8781870 of the 8781870 non-host reads" in result.stderr
    assert "processed 8781870 reads" in result.stderr


def test_paired_run_passes_and_keeps_the_log(tmp_path):
    result = run_check(tmp_path, "[M::bam2fq_mainloop] discarded 0 singletons\n"
                                 "[M::bam2fq_mainloop] processed 2000 reads\n")
    assert result.returncode == 0, result.stderr
    assert "processed 2000 reads" in result.stderr


def test_a_handful_of_singletons_passes(tmp_path):
    result = run_check(tmp_path, "[M::bam2fq_mainloop] discarded 3 singletons\n"
                                 "[M::bam2fq_mainloop] processed 2000 reads\n")
    assert result.returncode == 0, result.stderr


def test_no_non_host_reads_passes(tmp_path):
    result = run_check(tmp_path, "[M::bam2fq_mainloop] discarded 0 singletons\n"
                                 "[M::bam2fq_mainloop] processed 0 reads\n")
    assert result.returncode == 0, result.stderr

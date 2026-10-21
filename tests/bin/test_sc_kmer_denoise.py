"""bin/sc_kmer_denoise.py: barcode-level k-mer correlation (SAHMI's single-cell denoising)."""

import csv
import random

import pytest

KMER_LEN = 8
READ_LEN = 30
KMERS_PER_READ = READ_LEN - KMER_LEN + 1
BARCODES = [f"CELL{i}" for i in range(1, 6)]


def random_sequence(rng):
    return "".join(rng.choice("ACGT") for _ in range(READ_LEN))


def write_inputs(tmp_path):
    """Taxon 100 scales: barcode i holds i DIFFERENT reads, so its distinct k-mers
    rise with its total. Taxon 200 saturates: barcode i holds i copies of ONE
    read, so the total rises and the distinct count never moves."""
    rng = random.Random(7)
    real_reads = [random_sequence(rng) for _ in BARCODES]
    ambient_read = random_sequence(rng)

    kraken, fastq = [], []
    serial = 0
    for index, barcode in enumerate(BARCODES, start=1):
        for taxid, sequences in ((100, real_reads[:index]), (200, [ambient_read] * index)):
            for sequence in sequences:
                serial += 1
                read_id = f"r{serial}|CB:{barcode}|UB:UMI{serial}"
                kraken.append(f"C\t{read_id}\t{taxid}\t{READ_LEN}\t{taxid}:{KMERS_PER_READ}\n")
                fastq.append(f"@{read_id}\n{sequence}\n+\n{'I' * READ_LEN}\n")

    reads = tmp_path / "sc1.kraken2.output.txt"
    reads.write_text("".join(kraken), encoding="utf-8")
    fq = tmp_path / "sc1.fastq"
    fq.write_text("".join(fastq), encoding="utf-8")
    return reads, fq


def run_denoise(run_script, tmp_path, adjust):
    reads, fq = write_inputs(tmp_path)
    result = run_script(
        "sc_kmer_denoise.py",
        "--reads", reads,
        "--fastq", fq,
        "--sample", "sc1",
        "--evidence", tmp_path / "evidence.tsv",
        "--drop-list", tmp_path / "drop.txt",
        "--mqc", tmp_path / "sc_kmer_mqc.tsv",
        "--kmer-len", KMER_LEN,
        "--adjust", adjust,
    )
    assert result.returncode == 0, result.stderr
    return result


@pytest.mark.parametrize(
    "adjust, named, not_named",
    [("holm", "Holm-adjusted", "BH-adjusted"), ("bh", "BH-adjusted", "Holm-adjusted")],
)
def test_mqc_header_names_the_correction_used(tmp_path, run_script, adjust, named, not_named):
    """The MultiQC section must describe the multiple-testing correction that was
    actually applied, not a fixed one."""
    run_denoise(run_script, tmp_path, adjust)

    text = (tmp_path / "sc_kmer_mqc.tsv").read_text(encoding="utf-8")
    header = "\n".join(line for line in text.splitlines() if line.startswith("#"))
    assert named in header
    assert not_named not in header
    data = [line for line in text.splitlines() if not line.startswith("#")]
    assert data == ["Sample\tkept\tno_barcode_correlation", "sc1\t1\t1"]


def test_scaling_taxon_kept_and_saturating_taxon_dropped(tmp_path, run_script):
    """A taxon whose distinct k-mers rise with its total across barcodes is kept;
    one whose total rises while its distinct k-mers stay flat (ambient /
    single-locus signature) is dropped as saturated."""
    run_denoise(run_script, tmp_path, "holm")

    with open(tmp_path / "evidence.tsv", encoding="utf-8", newline="") as handle:
        evidence = {row["taxid"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert evidence["100"]["verdict"] == "kept"
    assert evidence["100"]["barcodes"] == "5"
    assert float(evidence["100"]["rho_kmers_distinct"]) == 1.0
    assert evidence["200"]["verdict"] == "saturated_distinct"
    assert (tmp_path / "drop.txt").read_text(encoding="utf-8").split() == ["200"]

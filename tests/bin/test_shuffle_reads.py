"""bin/shuffle_reads.py: composition-preserving shuffles for the shuffled-read negative control."""

import gzip
import random
from collections import Counter


def make_reads(seed, count=20, length=75):
    rng = random.Random(seed)
    reads = []
    for index in range(count):
        sequence = "".join(rng.choice("ACGT") for _ in range(length))
        quality = "".join(rng.choice("#,:FI") for _ in range(length))
        reads.append((f"@read{index} extra=comment", sequence, "+", quality))
    return reads


def write_fastq(path, reads, gz=False):
    text = "".join(f"{name}\n{sequence}\n{plus}\n{quality}\n" for name, sequence, plus, quality in reads)
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def read_fastq_gz(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [tuple(lines[i:i + 4]) for i in range(0, len(lines), 4)]


def dinucleotides(sequence):
    return Counter(zip(sequence, sequence[1:]))


def test_dinuc_preserves_each_reads_dinucleotide_counts(tmp_path, run_script):
    """Altschul-Erickson: every shuffled read has exactly the dinucleotide counts
    (and so the base composition) of its source read, mates stay paired, and the
    sequence is actually rearranged."""
    mate1, mate2 = make_reads(1), make_reads(2)
    write_fastq(tmp_path / "in_1.fastq.gz", mate1, gz=True)
    write_fastq(tmp_path / "in_2.fastq.gz", mate2, gz=True)
    prefix = tmp_path / "S1"

    result = run_script(
        "shuffle_reads.py",
        "--reads", tmp_path / "in_1.fastq.gz", tmp_path / "in_2.fastq.gz",
        "--prefix", prefix,
        "--method", "dinuc",
        "--seed", 42,
    )

    assert result.returncode == 0, result.stderr
    for source, shuffled_path in ((mate1, f"{prefix}_shuffled_1.fastq.gz"), (mate2, f"{prefix}_shuffled_2.fastq.gz")):
        shuffled = read_fastq_gz(shuffled_path)
        assert len(shuffled) == len(source)
        for (name, sequence, plus, quality), (new_name, new_sequence, new_plus, new_quality) in zip(source, shuffled):
            assert new_name == name
            assert new_plus == plus
            assert new_quality == quality
            assert len(new_sequence) == len(sequence)
            assert dinucleotides(new_sequence) == dinucleotides(sequence)
            assert Counter(new_sequence) == Counter(sequence)
        assert any(new[1] != old[1] for old, new in zip(source, shuffled))


def test_mono_preserves_base_composition(tmp_path, run_script):
    """The plain permutation keeps each read's length and base counts (so its GC)
    and carries names and qualities through unchanged."""
    reads = make_reads(3)
    write_fastq(tmp_path / "in.fastq", reads)
    prefix = tmp_path / "S1"
    stats = tmp_path / "stats.tsv"

    result = run_script(
        "shuffle_reads.py",
        "--reads", tmp_path / "in.fastq",
        "--prefix", prefix,
        "--method", "mono",
        "--seed", 7,
        "--stats", stats,
    )

    assert result.returncode == 0, result.stderr
    shuffled = read_fastq_gz(f"{prefix}_shuffled.fastq.gz")
    assert len(shuffled) == len(reads)
    for (name, sequence, _plus, quality), (new_name, new_sequence, _new_plus, new_quality) in zip(reads, shuffled):
        assert new_name == name
        assert new_quality == quality
        assert Counter(new_sequence) == Counter(sequence)
    assert any(new[1] != old[1] for old, new in zip(reads, shuffled))

    header, values = stats.read_text(encoding="utf-8").splitlines()
    record = dict(zip(header.split("\t"), values.split("\t")))
    assert record["method"] == "mono"
    assert int(record["fragments_shuffled"]) == len(reads)
    assert int(record["bases"]) == sum(len(read[1]) for read in reads)


def test_same_seed_gives_the_same_shuffle(tmp_path, run_script):
    """A control that changes between runs is not a control: the same seed must
    reproduce the output exactly, and a different seed must not."""
    reads = make_reads(4)
    write_fastq(tmp_path / "in.fastq", reads)

    outputs = {}
    for label, seed in (("a", 11), ("b", 11), ("c", 12)):
        prefix = tmp_path / label
        result = run_script(
            "shuffle_reads.py",
            "--reads", tmp_path / "in.fastq",
            "--prefix", prefix,
            "--method", "dinuc",
            "--seed", seed,
        )
        assert result.returncode == 0, result.stderr
        outputs[label] = read_fastq_gz(f"{prefix}_shuffled.fastq.gz")

    assert outputs["a"] == outputs["b"]
    assert outputs["a"] != outputs["c"]

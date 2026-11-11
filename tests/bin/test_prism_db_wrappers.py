"""bin/prism_blastdbcmd.py and bin/prism_makeblastdb.py: PRISM's per-library BLAST
database extracted volume by volume and built chunk by chunk, with fake tools."""

import stat


def fake_tool(path, body):
    path.write_text("#!/usr/bin/env bash\n" + body)
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


def fake_blastdbcmd(tmp_path, missing=(), broken=()):
    """Writes its volume's name to -out; `missing` volumes report no taxids, `broken` ones fail."""
    return fake_tool(tmp_path / "blastdbcmd", f"""
db=; out=
while [ $# -gt 0 ]; do case $1 in -db) db=$2; shift;; -out) out=$2; shift;; esac; shift; done
name=$(basename "$db")
case " {' '.join(missing)} " in *" $name "*) : > "$out"; echo "Error: [blastdbcmd] Taxonomy ID(s) not found in the $name database." >&2; exit 1;; esac
case " {' '.join(broken)} " in *" $name "*) echo "Error: [blastdbcmd] disk on fire" >&2; exit 2;; esac
echo "$name" > "$out"
""")


def fake_makeblastdb(tmp_path, fail=()):
    """Copies its input (-in, or standard input for -) to <out>.nsq and its -title to
    <out>.title; `fail` database names exit 1."""
    return fake_tool(tmp_path / "makeblastdb", f"""
in=; out=; title=
while [ $# -gt 0 ]; do case $1 in -in) in=$2; shift;; -out) out=$2; shift;; -title) title=$2; shift;; esac; shift; done
case " {' '.join(fail)} " in *" $(basename "$out") "*) cat > /dev/null; exit 1;; esac
if [ "$in" = - ]; then cat > "$out.nsq"; else cat "$in" > "$out.nsq"; fi
printf '%s' "$title" > "$out.title"
""")


def argv_recorder(tmp_path, name):
    return fake_tool(tmp_path / name, 'printf "%s\\n" "$@" > "$(dirname "$0")/argv.txt"\n')


def write_fasta(path, records, length=60):
    path.write_text("".join(f">seq{i} title {i}\n{'ACGT' * (length // 4)}\n" for i in range(records)))


def test_blastdbcmd_extracts_volumes_in_alias_order(tmp_path, run_script):
    """One blastdbcmd per volume, outputs joined in DBLIST order; a volume without the
    taxids contributes nothing and does not fail the call."""
    db = tmp_path / "nt"
    db.with_suffix(".nal").write_text("TITLE nt\nDBLIST nt.00 nt.01 nt.02 nt.03\n")
    tool = fake_blastdbcmd(tmp_path, missing=("nt.02",))
    out = tmp_path / "custom.fasta"
    result = run_script("prism_blastdbcmd.py", tool, 2, "-db", db, "-dbtype", "nucl", "-outfmt", "%f", "-taxidlist", "ids.txt", "-out", out)
    assert result.returncode == 0, result.stderr
    assert out.read_text() == "nt.00\nnt.01\nnt.03\n"
    assert "4 volumes, 3 with sequences" in result.stderr
    assert not list(tmp_path.glob("custom.fasta.volumes.*"))


def test_blastdbcmd_fails_as_blastdbcmd_would_when_no_volume_has_the_taxids(tmp_path, run_script):
    db = tmp_path / "nt"
    db.with_suffix(".nal").write_text("DBLIST nt.00 nt.01\n")
    tool = fake_blastdbcmd(tmp_path, missing=("nt.00", "nt.01"))
    result = run_script("prism_blastdbcmd.py", tool, 2, "-db", db, "-taxids", "1", "-out", tmp_path / "x.fasta")
    assert result.returncode == 1
    assert "Taxonomy ID(s) not found" in result.stderr


def test_blastdbcmd_reports_a_real_failure(tmp_path, run_script):
    db = tmp_path / "nt"
    db.with_suffix(".nal").write_text("DBLIST nt.00 nt.01\n")
    tool = fake_blastdbcmd(tmp_path, broken=("nt.01",))
    result = run_script("prism_blastdbcmd.py", tool, 2, "-db", db, "-taxids", "1", "-out", tmp_path / "x.fasta")
    assert result.returncode == 2
    assert "disk on fire" in result.stderr and "failed on nt.01" in result.stderr


def test_blastdbcmd_passes_other_calls_through(tmp_path, run_script):
    """No alias file, an -entry lookup, or another -outfmt: the call reaches blastdbcmd as given."""
    tool = argv_recorder(tmp_path, "blastdbcmd")
    db = tmp_path / "nt"
    db.with_suffix(".nal").write_text("DBLIST nt.00 nt.01\n")
    cases = [
        ["-db", str(tmp_path / "plain"), "-taxids", "1", "-out", "o"],
        ["-db", str(db), "-entry", "ACC1", "-out", "o"],
        ["-db", str(db), "-taxids", "1", "-outfmt", "%a", "-out", "o"],
    ]
    for args in cases:
        result = run_script("prism_blastdbcmd.py", tool, 4, *args)
        assert result.returncode == 0, result.stderr
        assert (tmp_path / "argv.txt").read_text().split("\n")[:-1] == args


def test_makeblastdb_streams_record_aligned_ranges_and_an_alias_in_order(tmp_path, run_script, monkeypatch):
    """Each database gets a run of whole records, the alias lists them in file order,
    and together they hold the FASTA byte for byte."""
    monkeypatch.setenv("PRISM_MAKEBLASTDB_MIN_BYTES", "1")
    monkeypatch.setenv("PRISM_MAKEBLASTDB_CHUNK_BYTES", "2000")
    fasta = tmp_path / "customdb.fasta"
    write_fasta(fasta, 100)
    tool = fake_makeblastdb(tmp_path)
    out = tmp_path / "db" / "custom_taxa_db"
    out.parent.mkdir()
    result = run_script("prism_makeblastdb.py", tool, 3, "-in", fasta, "-dbtype", "nucl", "-parse_seqids", "-taxid_map", "map.txt", "-out", out)
    assert result.returncode == 0, result.stderr
    alias = out.with_suffix(".nal").read_text().splitlines()
    names = alias[1].split()[1:]
    assert alias[0] == f"TITLE {fasta}"
    assert names == [f"custom_taxa_db.p{k:03d}" for k in range(len(names))]
    assert len(names) >= 3
    parts = [(out.parent / (name + ".nsq")).read_bytes() for name in names]
    assert b"".join(parts) == fasta.read_bytes()
    assert all(part.startswith(b">") and part.endswith(b"\n") for part in parts)
    assert {(out.parent / (name + ".title")).read_text() for name in names} == {str(fasta)}
    assert list(tmp_path.rglob("*.fasta")) == [fasta], "no chunk files are written"


def test_makeblastdb_ranges_cover_long_records(tmp_path, run_script, monkeypatch):
    """Records longer than a share of the file: no range starts inside one, and none is empty."""
    monkeypatch.setenv("PRISM_MAKEBLASTDB_MIN_BYTES", "1")
    fasta = tmp_path / "customdb.fasta"
    fasta.write_text(">a\n" + "ACGT" * 500 + "\n>b\nAC\n>c\n" + "\n".join(["ACGTACGT"] * 300) + "\n")
    tool = fake_makeblastdb(tmp_path)
    result = run_script("prism_makeblastdb.py", tool, 8, "-in", fasta, "-dbtype", "nucl", "-out", tmp_path / "custom_taxa_db")
    assert result.returncode == 0, result.stderr
    names = (tmp_path / "custom_taxa_db.nal").read_text().splitlines()[1].split()[1:]
    parts = [(tmp_path / (name + ".nsq")).read_bytes() for name in names]
    assert b"".join(parts) == fasta.read_bytes()
    assert all(part.startswith(b">") for part in parts)
    assert 2 <= len(parts) <= 3


def test_makeblastdb_small_input_passes_through(tmp_path, run_script):
    fasta = tmp_path / "small.fasta"
    write_fasta(fasta, 3)
    tool = argv_recorder(tmp_path, "makeblastdb")
    args = ["-in", str(fasta), "-dbtype", "nucl", "-out", str(tmp_path / "db")]
    result = run_script("prism_makeblastdb.py", tool, 4, *args)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "argv.txt").read_text().split("\n")[:-1] == args
    assert not (tmp_path / "db.nal").exists()


def test_makeblastdb_reports_a_failed_chunk(tmp_path, run_script, monkeypatch):
    monkeypatch.setenv("PRISM_MAKEBLASTDB_MIN_BYTES", "1")
    fasta = tmp_path / "customdb.fasta"
    write_fasta(fasta, 20)
    tool = fake_makeblastdb(tmp_path, fail=("custom_taxa_db.p001",))
    result = run_script("prism_makeblastdb.py", tool, 2, "-in", fasta, "-dbtype", "nucl", "-out", tmp_path / "custom_taxa_db")
    assert result.returncode == 1
    assert "failed on 1 of 2 range(s)" in result.stderr
    assert not (tmp_path / "custom_taxa_db.nal").exists()

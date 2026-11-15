"""bin/exploremetatax_bundle.py: HUMAnN tables renamed for the exploreMetaTax app."""

import gzip
import tarfile
from pathlib import Path

import pytest

PATHABUNDANCE = "# Pathway\tS_Abundance\nPWY-101: example pathway\t12.5\nUNMAPPED\t100.0\n"
GENEFAMILIES = "# Gene Family\tS_Abundance-RPKs\nUniRef90_A0A000\t3.0\nUNMAPPED\t50.0\n"
REGROUP = "# Gene Family\tS_Abundance-RPKs\nK00001\t7.0\nUNGROUPED\t46.0\n"


def write_humann(directory):
    """S1 arrives gzipped as the pipeline publishes it, S2 as plain text."""
    directory.mkdir()
    for name, text in (("S1_pathabundance.tsv", PATHABUNDANCE),
                       ("S1_genefamilies.tsv", GENEFAMILIES),
                       ("S1_regroup.tsv", REGROUP)):
        with gzip.open(directory / f"{name}.gz", "wt") as handle:
            handle.write(text)
    for name, text in (("S2_pathabundance.tsv", PATHABUNDANCE),
                       ("S2_genefamilies.tsv", GENEFAMILIES),
                       ("S2_regroup.tsv", REGROUP)):
        (directory / name).write_text(text, encoding="utf-8")


def bundle(run_script, tmp_path, group):
    write_humann(tmp_path / "humann")
    output = tmp_path / "reanatax.exploreMetaTax.tar.gz"
    args = [
        "--humann", tmp_path / "humann",
        "--outdir", tmp_path / "reanatax_exploreMetaTax",
        "--output", output,
    ]
    if group is not None:
        args += ["--humann-regroup", group]
    result = run_script("exploremetatax_bundle.py", *args, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    with tarfile.open(output, "r:gz") as tar:
        members = {Path(m.name).name: tar.extractfile(m).read().decode("utf-8")
                   for m in tar.getmembers() if m.isfile()}
        top = {Path(m.name).parts[0] for m in tar.getmembers()}
    assert top == {"reanatax_exploreMetaTax"}
    return members


BASE_FILES = {"README.txt", "S1_pathabundance.tsv", "S1_genefamilies.tsv",
              "S2_pathabundance.tsv", "S2_genefamilies.tsv"}


@pytest.mark.parametrize(
    "group, suffix",
    [
        (None, "_ko.tsv"),
        ("uniref90_ko", "_ko.tsv"),
        ("uniref50_ko", "_ko.tsv"),
        ("uniref90_rxn", "_reactions.tsv"),
    ],
)
def test_regroup_table_is_named_after_its_grouping(tmp_path, run_script, group, suffix):
    """`<sample>_regroup.tsv(.gz)` is bundled as `<sample>_ko.tsv` for a *_ko
    grouping (the default is uniref90_ko) and `<sample>_reactions.tsv` for *_rxn,
    decompressed; path abundance and gene families keep their names."""
    members = bundle(run_script, tmp_path, group)

    assert set(members) == BASE_FILES | {f"S1{suffix}", f"S2{suffix}"}
    assert members[f"S1{suffix}"] == REGROUP
    assert members[f"S2{suffix}"] == REGROUP
    assert members["S1_pathabundance.tsv"] == PATHABUNDANCE
    assert members["S1_genefamilies.tsv"] == GENEFAMILIES


@pytest.mark.parametrize("group", ["uniref90_go", "uniref90_level4ec", "uniref90_pfam"])
def test_regroup_table_without_an_app_slot_is_left_out(tmp_path, run_script, group):
    """A grouping the app has no slot for (GO, EC, Pfam, ...) is left out of the
    bundle rather than shown under the wrong name, and that is not an error."""
    members = bundle(run_script, tmp_path, group)

    assert set(members) == BASE_FILES
    assert not any(name.endswith(("_ko.tsv", "_reactions.tsv", "_regroup.tsv")) for name in members)


ANCOMBC = "taxon\tlfc\tpvalue\tqvalue\tsignificant\tcomparison\nStaphylococcus aureus\t2.5\t0.001\t0.01\tTRUE\tcase_over_control\n"
EMU_HEADER = "tax_id\tspecies\tgenus\tsuperkingdom\tS1\tS2\n"


def test_da_tables_and_one_emu_table_are_bundled_as_side_files(tmp_path, run_script):
    """The ANCOM-BC2 / ALDEx2 result tables keep their names, which the app loads
    on its own and which match no report format; of the per-rank Emu tables only
    the deepest rank goes in, its count variant preferred, and metadata.tsv is
    labelled as loaded by itself."""
    da = tmp_path / "da"
    da.mkdir()
    (da / "bracken_combined.ancombc.results.tsv").write_text(ANCOMBC, encoding="utf-8")
    (da / "bracken_combined.aldex2.results.tsv").write_text(ANCOMBC, encoding="utf-8")
    emu = tmp_path / "emu"
    emu.mkdir()
    for name in ("emu-combined-species.tsv", "emu-combined-species-counts.tsv",
                 "emu-combined-genus.tsv", "emu-combined-genus-counts.tsv"):
        (emu / name).write_text(EMU_HEADER + "1280\tStaphylococcus aureus\tStaphylococcus\tBacteria\t1\t1\n",
                                encoding="utf-8")
    metadata = tmp_path / "samples.tsv"
    metadata.write_text("SampleID\tgroup\nS1\tcase\nS2\tcontrol\n", encoding="utf-8")
    output = tmp_path / "reanatax.exploreMetaTax.tar.gz"
    result = run_script(
        "exploremetatax_bundle.py",
        "--da", da, "--emu", emu, "--metadata", metadata,
        "--outdir", tmp_path / "reanatax_exploreMetaTax", "--output", output,
        cwd=tmp_path,
    )
    assert result.returncode == 0, result.stderr
    with tarfile.open(output, "r:gz") as tar:
        members = {Path(m.name).name: tar.extractfile(m).read().decode("utf-8")
                   for m in tar.getmembers() if m.isfile()}
    assert set(members) == {"README.txt", "bracken_combined.ancombc.results.tsv",
                            "bracken_combined.aldex2.results.tsv",
                            "emu-combined-species-counts.tsv", "metadata.tsv"}
    assert members["bracken_combined.ancombc.results.tsv"] == ANCOMBC
    assert "Upload this .tar.gz as it is" in members["README.txt"]
    assert "emu-combined-genus" in result.stderr


def test_a_da_table_with_an_unrecognised_name_fails_the_bundle(tmp_path, run_script):
    """A result table the app would not recognise is an error, not a silent
    passenger in the archive."""
    da = tmp_path / "da"
    da.mkdir()
    (da / "results.tsv").write_text(ANCOMBC, encoding="utf-8")
    result = run_script(
        "exploremetatax_bundle.py", "--da", da,
        "--outdir", tmp_path / "reanatax_exploreMetaTax",
        "--output", tmp_path / "out.tar.gz", cwd=tmp_path,
    )
    assert result.returncode != 0
    assert "differential-abundance" in result.stderr

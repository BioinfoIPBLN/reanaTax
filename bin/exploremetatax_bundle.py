#!/usr/bin/env python3
"""
exploremetatax_bundle.py -- package the taxonomic tables for the exploreMetaTax
Shiny app (https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax).

The app has no archive ingest. You extract this tarball and point its "Browse a
folder" upload at the extracted directory; it then uploads every file in the
folder and keeps only those whose NAME matches the format picked on the radio
button. So the bundle is a flat directory whose filenames are chosen to satisfy
that filter, not a copy of the results tree.

That distinction is not cosmetic. Taking reanaTax's published names as they are:

  * `<sample>.bracken_S.tsv` -- the Bracken abundance table the app actually
    wants -- matches NO format pattern and would be silently skipped.
  * `<sample>.bracken_S.kraken2.report_bracken.txt` matches the `bracken`
    pattern, but it is a kreport, not an abundance table, so the app would hand
    it to a parser looking for `new_est_reads` columns.
  * `bracken_combined_S.txt` misses `combined_bracken`, whose regex wants
    "combined" BEFORE "bracken".
  * `kraken2_combined_report.txt` false-matches the per-sample `kraken2`
    pattern and would be loaded as if it were one more sample.

So files are renamed on the way in, the two that can only mislead are left out,
and every name written is checked against every pattern before the tarball is
sealed: a name must match the one format it is meant for and no other. Change a
name here without checking and the build fails rather than shipping an archive
that loads wrongly.
"""
import argparse
import gzip
import os
import re
import shutil
import sys
import tarfile

# Verbatim from format_filename_pattern() in apps/exploreMetaTax/app/app.R.
# Kept as a literal copy so a divergence shows up as a failed bundle rather
# than as a file the app quietly ignores.
APP_PATTERNS = {
    "kraken2": r"\.report\.txt$|_report\.txt$|\.kreport2?$|_kraken2\.txt$",
    "krakenuniq": r"_krakenuniq.*\.txt$|\.krakenuniq$|_krakenuniq_report\.txt$",
    "bracken": r"\.bracken$|_bracken[^/]*\.(txt|tsv)$",
    "combined_bracken": r"combined.*bracken.*\.(txt|tsv)$",
    "metaphlan": r"_profile\.txt$|_metaphlan.*\.(txt|tsv)$|_metaphlan_bugs_list\.tsv$",
    "merged_metaphlan": r"merged_abundance.*\.(txt|tsv)$|merged_metaphlan.*\.(txt|tsv)$",
    "humann": (
        r"pathabund(ance)?.*\.tsv$|reactions?.*\.tsv$|gene[_-]?famil(ies|y).*\.tsv$"
        r"|(^|[_-])ko(_cpm)?\.tsv$|kegg[_-]?kos?.*\.tsv$"
    ),
}

# Which radio button each of our files is meant to be loaded under. `None`
# means the file is documentation or metadata and must match no format at all.
LAYOUT_NOTE = {
    "kraken2": "Kraken2 Reports",
    "krakenuniq": "KrakenUniq Reports",
    "bracken": "Bracken Per-Sample",
    "combined_bracken": "Combined Bracken Table",
    "metaphlan": "MetaPhlAn Profiles",
    "merged_metaphlan": "Merged MetaPhlAn Table",
    "humann": "HUMAnN Tables",
}


def matching_formats(name):
    return {fmt for fmt, pattern in APP_PATTERNS.items() if re.search(pattern, name, re.I)}


def sample_of(name, *suffixes):
    """Strip a reanaTax suffix to recover the sample id."""
    for suffix in sorted(suffixes, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return re.sub(r"\.(tsv|txt|gz)$", "", name)


def files_in(directory):
    if not directory or not os.path.isdir(directory):
        return []
    return sorted(
        os.path.join(directory, entry)
        for entry in os.listdir(directory)
        if os.path.isfile(os.path.join(directory, entry))
    )


def place(src, dest_dir, dest_name, wanted_format, planned):
    """Copy one file under its app-facing name, decompressing .gz on the way."""
    dest = os.path.join(dest_dir, dest_name)
    if src.endswith(".gz") and not dest_name.endswith(".gz"):
        with gzip.open(src, "rb") as fin, open(dest, "wb") as fout:
            shutil.copyfileobj(fin, fout)
    else:
        shutil.copyfile(src, dest)
    planned.append((dest_name, wanted_format))
    return dest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kraken2", help="directory of per-sample Kraken2 reports")
    parser.add_argument("--krakenuniq", help="directory of per-sample KrakenUniq reports")
    parser.add_argument("--bracken", help="directory of per-sample Bracken tables")
    parser.add_argument("--bracken-combined", help="directory holding the combined Bracken table")
    parser.add_argument("--metaphlan", help="directory of per-sample MetaPhlAn profiles")
    parser.add_argument("--metaphlan-merged", help="directory holding the merged MetaPhlAn table")
    parser.add_argument("--humann", help="directory of HUMAnN per-sample tables")
    parser.add_argument("--metadata", help="sample metadata TSV")
    parser.add_argument("--outdir", default="exploreMetaTax", help="staging directory name")
    parser.add_argument("--output", required=True, help="tarball to write")
    args = parser.parse_args()

    stage = args.outdir
    os.makedirs(stage, exist_ok=True)
    planned = []

    # Kraken2 per-sample reports already end `.report.txt`, which the app's
    # kraken2 pattern accepts and no other pattern touches.
    for path in files_in(args.kraken2):
        place(path, stage, os.path.basename(path), "kraken2", planned)

    # `<sample>.krakenuniq.report.txt` -> `<sample>.krakenuniq`. The obvious
    # `<sample>_krakenuniq_report.txt` would satisfy the app's kraken2 pattern
    # too (`_report\.txt$`) and load as an extra Kraken2 sample.
    for path in files_in(args.krakenuniq):
        name = os.path.basename(path)
        sample = re.sub(r"\.krakenuniq\.report\.txt$", "", name)
        place(path, stage, f"{sample}.krakenuniq", "krakenuniq", planned)

    # `<sample>.bracken_<lvl>.tsv` -> `<sample>_bracken.tsv`, so the `_bracken`
    # the app looks for is actually present.
    for path in files_in(args.bracken):
        name = os.path.basename(path)
        sample = re.sub(r"\.bracken_[A-Z][0-9]*\.tsv$", "", name)
        place(path, stage, f"{sample}_bracken.tsv", "bracken", planned)

    # Hyphens, not underscores: `combined_bracken.tsv` would ALSO satisfy the
    # per-sample `_bracken` pattern and be loaded as a sample.
    for path in files_in(args.bracken_combined):
        place(path, stage, "combined-bracken-table.tsv", "combined_bracken", planned)

    for path in files_in(args.metaphlan):
        place(path, stage, os.path.basename(path), "metaphlan", planned)

    # Likewise `merged_metaphlan.tsv` would satisfy the per-sample `_metaphlan`
    # pattern; `merged_abundance_table.tsv` satisfies only the merged one.
    for path in files_in(args.metaphlan_merged):
        place(path, stage, "merged_abundance_table.tsv", "merged_metaphlan", planned)

    # HUMAnN tables arrive gzipped and are decompressed by place(). Only the
    # three the app can name are taken: `_regroup` is the KO table (the
    # pipeline regroups onto UniRef90->KO by default), while `_renorm` is
    # deliberately NOT bundled - it is the renormalisation of whichever table
    # came before it, so its contents depend on --humann_regroup and its name
    # alone cannot say whether it holds KOs or gene families.
    for path in files_in(args.humann):
        name = os.path.basename(path)
        base = name[:-3] if name.endswith(".gz") else name
        if base.endswith("_regroup.tsv"):
            base = f"{sample_of(base, '_regroup.tsv')}_ko.tsv"
        elif not (base.endswith("_pathabundance.tsv") or base.endswith("_genefamilies.tsv")):
            raise SystemExit(
                f"exploremetatax_bundle: unrecognised HUMAnN table '{name}'. Expected "
                "*_pathabundance.tsv, *_genefamilies.tsv or *_regroup.tsv."
            )
        place(path, stage, base, "humann", planned)

    if args.metadata and os.path.isfile(args.metadata):
        place(args.metadata, stage, "metadata.tsv", None, planned)

    if not planned:
        raise SystemExit("exploremetatax_bundle: nothing to bundle.")

    # The guard. Every name must resolve to exactly the one format it is for -
    # a name matching two formats gets loaded under both, which is how a
    # combined table ends up masquerading as a sample.
    problems = []
    for name, wanted in planned:
        matched = matching_formats(name)
        expected = {wanted} if wanted else set()
        if matched != expected:
            problems.append(
                f"  {name}: matches {sorted(matched) or ['nothing']}, expected {sorted(expected) or ['nothing']}"
            )
    if problems:
        raise SystemExit(
            "exploremetatax_bundle: filenames would be mis-loaded by exploreMetaTax:\n"
            + "\n".join(problems)
        )

    counts = {}
    for _name, wanted in planned:
        counts[wanted] = counts.get(wanted, 0) + 1

    with open(os.path.join(stage, "README.txt"), "w", encoding="utf-8") as handle:
        handle.write(
            "exploreMetaTax bundle\n"
            "=====================\n\n"
            "https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax\n\n"
            "1. Extract this archive.\n"
            "2. In the app, open 'Upload and Filter Data'.\n"
            "3. Tick 'Browse a folder (uploads all matching reports inside)'.\n"
            "4. Pick the input format on the radio button, then Browse to the\n"
            "   extracted folder and press 'Load Data'. Only the files belonging\n"
            "   to the chosen format are loaded, so the same folder serves every\n"
            "   format below - switch the radio button and load again.\n\n"
            "Contents:\n"
        )
        for fmt, count in sorted(counts.items(), key=lambda item: str(item[0])):
            if fmt is None:
                handle.write("  metadata.tsv - load via the separate 'Metadata' box\n")
            else:
                handle.write(f"  {count} file(s) -> radio button '{LAYOUT_NOTE[fmt]}'\n")

    with tarfile.open(args.output, "w:gz") as tar:
        tar.add(stage, arcname=os.path.basename(stage))

    print(
        f"[exploremetatax_bundle] wrote {args.output} with {len(planned)} file(s): "
        + ", ".join(f"{count} {fmt or 'metadata'}" for fmt, count in sorted(counts.items(), key=lambda i: str(i[0]))),
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

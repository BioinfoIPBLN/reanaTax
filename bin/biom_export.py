#!/usr/bin/env python3
"""
biom_export.py -- finish the BIOM table kraken-biom started.

kraken-biom turns per-sample kreports into a BIOM table and stops there. Three
things have to happen afterwards for the result to be usable and, more
importantly, for it to agree with everything else this pipeline publishes.

Filtering is the one that matters. kraken-biom reads the RAW per-sample
reports, so a BIOM built straight from them carries every taxon the evidence
filters, the negative-control filter and the abundance thresholds removed from
the combined tables. Shipping that next to the filtered tables would hand a
reader two different answers from the same run and no way to tell which is
which. --keep-table takes the filtered combined table and restricts the BIOM to
the taxids that survived it, so one filtering decision is expressed once.

Metadata is the second. kraken-biom's own --metadata goes through pandas, which
its bioconda package does not depend on, and indexes the frame with .loc on the
sample list, so a sample missing from the file dies in a pandas traceback. The
biom API needs neither, and a missing sample is worth a sentence rather than a
KeyError.

The third is the format. BIOM 1.0 (JSON) is the default here because it is what
phyloseq's import_biom reads without rhdf5, and because h5py is not guaranteed
to be in the container. --format hdf5 is available where it is, and says so
plainly when it is not.
"""

import argparse
import os
import sys

import biom

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from filter_abundance import read_table, taxid_column

GENERATED_BY = "reanaTax"


def keep_taxids(path):
    """The taxids still present in a filtered combined table."""
    _comments, _header_line, header, rows = read_table(path)
    column = taxid_column(header)
    if column is None:
        raise SystemExit(
            f"biom_export: {path} has no taxonomy_id/taxid column, so there is "
            "nothing to match BIOM observation ids against."
        )
    return {row[column].strip() for row in rows if column < len(row) and row[column].strip()}


def read_metadata(path, samples):
    """{sample: {column: value}} from a TSV whose first column is the sample id."""
    with open(path, encoding="utf-8") as handle:
        lines = [line.rstrip("\r\n") for line in handle if line.strip()]
    if not lines:
        raise SystemExit(f"biom_export: {path} is empty")
    header = lines[0].lstrip("#").split("\t")
    table = {}
    for line in lines[1:]:
        fields = line.split("\t")
        table[fields[0].strip()] = {
            name.strip(): (fields[i].strip() if i < len(fields) else "")
            for i, name in enumerate(header)
        }
    missing = [name for name in samples if name not in table]
    if missing:
        raise SystemExit(
            f"biom_export: {len(missing)} sample(s) classified in this run are absent "
            f"from {path}, and BIOM sample metadata must cover every sample: "
            + ", ".join(missing[:10])
            + (" ..." if len(missing) > 10 else "")
        )
    return table


def write(table, path, fmt):
    if fmt == "hdf5":
        try:
            from biom.util import biom_open
        except ImportError as error:  # pragma: no cover - depends on the container
            raise SystemExit(f"biom_export: --format hdf5 is unavailable here ({error}).")
        try:
            with biom_open(path, "w") as handle:
                table.to_hdf5(handle, GENERATED_BY)
        except RuntimeError as error:
            raise SystemExit(
                f"biom_export: writing HDF5 failed ({error}). This build of biom-format "
                "has no h5py; use the default --format json, or `biom convert` afterwards."
            )
        return
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(table.to_json(GENERATED_BY))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--biom", required=True, help="the table kraken-biom wrote")
    parser.add_argument("--output", required=True, help="where the finished table goes")
    parser.add_argument("--keep-table", help="filtered combined table; its taxids are the ones kept")
    parser.add_argument("--metadata", help="TSV, sample id in column 1, attached as sample metadata")
    parser.add_argument("--format", choices=("json", "hdf5"), default="json")
    parser.add_argument("--unfiltered", help="also write the pre-filter table here, if filtering changed anything")
    args = parser.parse_args()

    table = biom.load_table(args.biom)
    samples = list(table.ids(axis="sample"))
    before = table.shape[0]

    columns = 0
    if args.metadata:
        metadata = read_metadata(args.metadata, samples)
        columns = max((len(row) for row in metadata.values()), default=0)
        table.add_metadata({name: metadata[name] for name in samples}, axis="sample")

    if args.keep_table:
        keep = keep_taxids(args.keep_table)
        if args.unfiltered and before and len(keep & set(table.ids(axis="observation"))) < before:
            write(table, args.unfiltered, args.format)
        table.filter(lambda _values, taxid, _md: taxid in keep, axis="observation", inplace=True)
        table.remove_empty(axis="observation", inplace=True)

    after = table.shape[0]
    if after == 0:
        raise SystemExit(
            "biom_export: every observation was filtered out. The filtered combined table "
            "and the per-sample reports share no taxid, which usually means the wrong "
            "table was handed to --keep-table."
        )

    write(table, args.output, args.format)
    print(
        f"[biom_export] {after} taxa x {len(samples)} samples"
        + (f"; dropped {before - after} not in {os.path.basename(args.keep_table)}" if args.keep_table else "")
        + (f"; {columns} metadata column(s)" if args.metadata else ""),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Normalise a ``fastq-dl`` run-info table into a slim, parser-safe run sheet.

``fastq-dl --only-download-metadata`` writes *every* field the archive knows
about (the ENA portal API is queried with ``fields=all``) through
``csv.DictWriter``.  Free-text columns such as ``sample_title`` or
``description`` are therefore quoted and may legitimately contain tabs, double
quotes and even embedded newlines.  Nextflow's ``splitCsv`` is line-based and
would mis-parse those rows, so this script re-emits only the handful of columns
the pipeline actually needs, with every whitespace run collapsed to a single
space and the separator switched to a comma.

One output row per sequencing run.  Missing columns are emitted as empty
strings rather than failing, because the SRA fallback provider in ``fastq-dl``
returns a slightly smaller field set than ENA does.
"""

import argparse
import csv
import re
import sys

# Columns copied straight through from the run-info table (in output order).
COLUMNS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "library_layout",
    "library_strategy",
    "library_source",
    "instrument_platform",
    "instrument_model",
    "scientific_name",
    "sample_title",
    "read_count",
    "base_count",
]

WHITESPACE = re.compile(r"\s+")

# csv.field_size_limit defaults to 128 kB; some ENA free-text fields exceed it.
csv.field_size_limit(min(sys.maxsize, 2**31 - 1))


def clean(value):
    """Collapse whitespace and strip characters that break downstream parsing."""
    if value is None:
        return ""
    return WHITESPACE.sub(" ", str(value)).strip()


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_info", help="fastq-dl '*-run-info.tsv' file")
    parser.add_argument("output", help="normalised run sheet to write (CSV)")
    parser.add_argument(
        "--query",
        default="",
        help="the accession that was queried, recorded in the 'query' column",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    with open(args.run_info, newline="", encoding="utf-8") as fin:
        rows = list(csv.DictReader(fin, delimiter="\t"))

    if not rows:
        sys.exit(f"ERROR: no runs found in '{args.run_info}' for query '{args.query}'")

    seen = set()
    with open(args.output, "w", newline="", encoding="utf-8") as fout:
        writer = csv.DictWriter(fout, fieldnames=["query"] + COLUMNS, delimiter=",")
        writer.writeheader()
        for row in rows:
            run = clean(row.get("run_accession"))
            if not run:
                continue
            # A run can be returned twice when several accessions in one query
            # resolve to overlapping sets of runs.
            if run in seen:
                continue
            seen.add(run)
            record = {column: clean(row.get(column)) for column in COLUMNS}
            record["query"] = clean(args.query)
            # Fall back to the run accession so that grouping never collapses
            # unrelated runs under an empty key.
            for key in ("experiment_accession", "sample_accession"):
                if not record[key]:
                    record[key] = run
            writer.writerow(record)

    if not seen:
        sys.exit(f"ERROR: no usable run accessions in '{args.run_info}'")


if __name__ == "__main__":
    main()

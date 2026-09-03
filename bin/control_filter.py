#!/usr/bin/env python3
"""
control_filter.py -- what a negative control can and cannot tell you.

Every other filter in this pipeline reaches a verdict about a TAXON: minimizer
breadth, host k-mer share, shuffled-read chance, decontam's prevalence score.
All of them are global - the taxon is real, or it is not - and that is exactly
the shape of answer this dataset cannot use.

The failure that motivated this was CSI-Microbes (Robinson et al. 2024): a
plexWell plate where 124 wells were deliberately infected with Fusobacterium
nucleatum and 34 were not. Fusobacterium appears in 13 of the 34, at 1-78 reads.
Those are false positives, and they are false positives made of PERFECTLY GOOD
Fusobacterium reads - index hopping, well-to-well carryover, ambient template on
a shared plate. There is nothing wrong with the reads. There is only something
wrong with the well they are in.

That rules out two whole families of filter at once:

  * Anything judging read quality or evidence breadth. A hopped read has the
    same minimizers as the read it hopped from. Measured on this plate, the
    distinct-minimizers-per-read ratio that Kraken2's manual and the
    exploreMetaTax app both recommend has AUC 0.415 - it ranks the false
    positives ABOVE the true positives, because the ratio is inversely tied to
    read count and every false positive sits at 1-78 reads.

  * Anything reaching a single verdict per taxon, decontam's prevalence method
    included. Fusobacterium is in 6 of the 12 blanks AND is the organism the
    experiment is about. A per-taxon call has to either delete it (recall 62.9%
    -> 0%) or spare it (specificity unchanged at 61.8%). Both are wrong, and no
    threshold escapes it, because the taxon is genuinely present in both places.

What works is a PER-SAMPLE test with an external reference level. A control
library measures how much of each taxon arrives without a sample: reagent
contamination, ambient template, and - on a plate - the hopping rate. A taxon in
a real library is kept only where it clears that level by a stated margin, and
zeroed where it does not. The same taxon can be signal in one library and
carryover in the next, which is the thing a global drop list cannot express and
the whole reason this exists.

Measured on that plate, against the mean control level with a 50x margin:
specificity 61.8% -> 100%, recall 62.9% -> 48.4%, balanced accuracy 62.3% ->
74.2%, and 87% of background reads gone. That recall cost is real and it is not
an artefact of the threshold: the wells it gives up are wells whose entire
Fusobacterium evidence is a handful of reads, indistinguishable by construction
from the carryover in the well next to it. The filter does not resolve that
ambiguity. It declines to call it.

Levels are compared as counts per million CLASSIFIED reads, never as raw reads.
Depth normalisation is the single largest improvement available on that plate -
reads-per-million reaches 73.5% specificity at full recall where raw reads
reaches 61.8% - and without it the threshold means a different thing in every
library.

--prevalence answers a different and much cruder question, for the cohorts that
have no controls at all: is this taxon in EVERY library? A reagent contaminant
introduced at extraction is; a biological signal usually is not. On the same
plate, dropping taxa present in all 180 libraries removes 87 taxa and 77.5% of
all microbial reads without touching Fusobacterium.

That is a large effect from a crude rule, and the rule is only sound where a
universally present taxon CANNOT be real. It is wrong for a mono-culture, wrong
for a dominant gut commensal, wrong for anything where one organism is expected
everywhere. It is off by default and it says what it removed.

Neither test replaces --decontam. decontam is given an external measurement of
the KIT and answers "is this taxon reagent"; this is given control LIBRARIES and
answers "is there more of this taxon here than arrives on its own". Where both
are available, run both.
"""
import argparse
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from filter_abundance import column_totals, read_table, taxid_column, to_float

# Per-sample column suffixes across the two combined layouts this pipeline
# produces. `_num`/`_all` carry the counts; `_frac`/`_lvl` travel with them and
# have to be zeroed in step or the table would contradict itself.
COUNT_SUFFIXES = ("_num", "_all")
COMPANION_SUFFIXES = ("_frac", "_lvl")
ANNOTATION = {
    "name", "taxid", "taxonomy_id", "taxonomy_lvl", "lvl_type",
    "perc", "tot_all", "tot_frac", "tot_lvl",
}


def sample_columns(header):
    """{sample: {'count': i, 'companion': [j, ...]}} keyed by column suffix.

    A sample owns more than one column - Bracken publishes `_num` beside
    `_frac`, combine_kreports `_all` beside `_lvl` - and a cell is only zeroed
    coherently if all of them go together.
    """
    samples = {}
    for index, field in enumerate(header):
        name = field.strip()
        if name.lower() in ANNOTATION:
            continue
        for suffix in COUNT_SUFFIXES:
            if name.endswith(suffix):
                samples.setdefault(name[: -len(suffix)], {})["count"] = index
                break
        else:
            for suffix in COMPANION_SUFFIXES:
                if name.endswith(suffix):
                    entry = samples.setdefault(name[: -len(suffix)], {})
                    entry.setdefault("companion", []).append(index)
                    break
    return {
        name: {"count": entry["count"], "companion": entry.get("companion", [])}
        for name, entry in samples.items()
        if "count" in entry
    }


def read_controls(value, known):
    """Control sample IDs, from a comma-separated list or a file of one per line.

    Every ID has to match a column, and an ID that matches nothing is fatal
    rather than ignored: a typo'd control name would silently turn the filter
    into a no-op on a cohort that has fewer controls than the user believes, and
    that is worse than not running it.
    """
    if os.path.exists(value):
        with open(value, encoding="utf-8") as handle:
            wanted = [line.strip() for line in handle if line.strip() and not line.startswith("#")]
    else:
        wanted = [entry.strip() for entry in value.split(",") if entry.strip()]

    missing = [entry for entry in wanted if entry not in known]
    if missing:
        raise SystemExit(
            "control_filter: no column for control sample(s) "
            + ", ".join(missing)
            + ". The table holds: "
            + ", ".join(sorted(known))
        )
    if not wanted:
        raise SystemExit("control_filter: --controls named no samples.")
    return wanted


LEVELS = {"mean": statistics.fmean, "median": statistics.median, "max": max}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", required=True, help="combined Bracken or Kraken2 table")
    parser.add_argument("--output", required=True, help="the same table with failing cells zeroed")
    parser.add_argument("--evidence", required=True, help="per-taxon evidence table to write")
    parser.add_argument("--drop-list", required=True,
                        help="taxids condemned in EVERY real sample, for the abundance filter")
    parser.add_argument("--drop-cells",
                        help="taxid<TAB>sample for every cell zeroed, so the same verdict can be "
                             "applied to objects this script never sees - the cell-by-taxon "
                             "matrix above all, which is built from the per-read assignments "
                             "rather than from this table and would otherwise keep every cell "
                             "the cohort tables just lost")
    parser.add_argument("--mqc", help="MultiQC custom-content section to write")
    parser.add_argument("--controls",
                        help="negative-control sample IDs: comma-separated, or a file of one "
                             "per line. Without this only --prevalence can run.")
    parser.add_argument("--ratio", type=float, default=10.0,
                        help="a taxon is kept in a sample only where it reaches this multiple "
                             "of its control level (counts per million classified reads)")
    parser.add_argument("--statistic", choices=sorted(LEVELS), default="mean",
                        help="which control level to compare against: mean (default), median, "
                             "or max. 'max' is the strict reading - clear the worst blank, not "
                             "the average one - and needs a smaller --ratio to do the same work.")
    parser.add_argument("--floor-reads", type=float, default=1,
                        help="the control level never falls below what this many reads in a "
                             "typical control library would be. A taxon absent from every "
                             "control has NOT been shown to be at zero - only below the "
                             "detection limit - and treating the two as the same lets any "
                             "taxon the blanks happened to miss clear any ratio on two reads. "
                             "0 disables, and restores that hole.")
    parser.add_argument("--min-reads", type=float, default=2,
                        help="a cell also needs this many reads before the ratio is consulted. "
                             "A single read cannot clear any margin meaningfully, and a taxon "
                             "absent from every control has a level of zero, which the ratio "
                             "test alone would wave through on one read.")
    parser.add_argument("--prevalence", type=float, default=0.0,
                        help="drop any taxon present in at least this FRACTION of all libraries "
                             "(1.0 = present in every one). 0 disables. Read the module "
                             "docstring before using it: it is only sound where a universally "
                             "present taxon cannot be real.")
    parser.add_argument("--prevalence-min-reads", type=float, default=1,
                        help="reads in a library for --prevalence to count it as present")
    args = parser.parse_args()

    if not args.controls and args.prevalence <= 0:
        raise SystemExit(
            "control_filter: nothing to do - neither --controls nor --prevalence was given."
        )

    comments, header_line, header, rows = read_table(args.input)
    samples = sample_columns(header)
    if not samples:
        raise SystemExit(f"control_filter: no per-sample columns found in {args.input}")

    controls = read_controls(args.controls, samples) if args.controls else []
    real = [name for name in samples if name not in set(controls)]
    if args.controls and not real:
        raise SystemExit("control_filter: every sample was named as a control.")

    # Counts per million CLASSIFIED reads. column_totals() already knows the two
    # denominators apart: a Bracken table is flat and its column sums to the
    # sample total, a combine_kreports table is a hierarchy whose column sums to
    # several times it, so there the root row's clade count is used instead.
    count_cols = [entry["count"] for entry in samples.values()]
    totals = column_totals(header, rows, count_cols)

    def cpm(row, column):
        if column >= len(row):
            return 0.0
        return 1e6 * to_float(row[column]) / (totals[column] or 1.0)

    def reads(row, column):
        return to_float(row[column]) if column < len(row) else 0.0

    # One read in a control library of typical depth, in the same counts-per-
    # million the levels are measured in.
    control_depth = (
        statistics.fmean([totals[samples[name]["count"]] for name in controls]) if controls else 0.0
    )
    floor_cpm = (
        1e6 * args.floor_reads / control_depth if (controls and control_depth and args.floor_reads > 0) else 0.0
    )

    taxid_col = taxid_column(header)

    def taxid_of(row):
        return row[taxid_col].strip() if taxid_col is not None and taxid_col < len(row) else ""

    name_col = next((i for i, h in enumerate(header) if h.strip().lower() == "name"), None)
    n_libraries = len(samples)
    prevalence_floor = args.prevalence * n_libraries if args.prevalence > 0 else None

    condemned, evidence, zeroed_cells, zeroed_reads, total_reads = [], [], 0, 0.0, 0.0
    zeroed_pairs = []
    counts = {}
    for row in rows:
        seen = [name for name in samples
                if reads(row, samples[name]["count"]) >= args.prevalence_min_reads]
        row_reads = sum(reads(row, entry["count"]) for entry in samples.values())
        total_reads += row_reads

        # Bound before the branch, because a row condemned by --prevalence never
        # reaches the control test and the evidence row is written for every
        # row regardless of which test spoke.
        level, kept_in, zeroed_in, judged = 0.0, [], [], False

        # --prevalence first, and it condemns the whole row: it is a statement
        # about the taxon across the cohort, not about any one library, so
        # zeroing it sample by sample would misrepresent the reason.
        if prevalence_floor is not None and len(seen) >= prevalence_floor:
            verdict = "ubiquitous"
        elif controls:
            level_values = [cpm(row, samples[name]["count"]) for name in controls]
            level = LEVELS[args.statistic](level_values) if level_values else 0.0
            # The detection limit, as a level in its own right. Absence from the
            # blanks is a bound, not a measurement: all it says is "under one
            # read in each of them". Without this floor the ratio test degrades
            # to `cpm > 0` for every taxon the controls happened to miss, which
            # on a species-level table is most of them - measured on the
            # CSI-Microbes plate, it left the specificity flat at 86.4% however
            # high --ratio went, because the species that leaked into the
            # uninfected wells was in one blank out of twelve.
            level = max(level, floor_cpm)
            judged = True
            for name in real:
                column = samples[name]["count"]
                if reads(row, column) <= 0:
                    continue
                passes = reads(row, column) >= args.min_reads and cpm(row, column) > args.ratio * level
                (kept_in if passes else zeroed_in).append(name)
            for name in zeroed_in:
                entry = samples[name]
                zeroed_pairs.append((taxid_of(row), name))
                zeroed_reads += reads(row, entry["count"])
                zeroed_cells += 1
                for column in [entry["count"]] + entry["companion"]:
                    if column < len(row):
                        row[column] = "0"
            if zeroed_in and not kept_in:
                verdict = "control_level_everywhere"
            elif zeroed_in:
                verdict = "control_level_somewhere"
            else:
                verdict = "kept"
        else:
            verdict = "kept"

        counts[verdict] = counts.get(verdict, 0) + 1
        taxid = taxid_of(row)
        if verdict in ("ubiquitous", "control_level_everywhere") and taxid:
            condemned.append(taxid)
        evidence.append([
            taxid,
            row[name_col].strip() if name_col is not None and name_col < len(row) else "",
            f"{row_reads:.0f}",
            len(seen),
            f"{len(seen) / n_libraries:.4f}",
            len(controls),
            # NA rather than 0 where the control test never ran: a row
            # condemned by --prevalence short-circuits before the level is
            # computed, and a printed 0.000 would read as a measurement.
            f"{level:.3f}" if judged else "NA",
            f"{args.ratio * level:.3f}" if judged else "NA",
            len(kept_in) if judged else "NA",
            len(zeroed_in) if judged else "NA",
            verdict,
        ])

    with open(args.output, "w", encoding="utf-8") as handle:
        for line in comments:
            handle.write(line + "\n")
        handle.write(header_line + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")

    with open(args.evidence, "w", encoding="utf-8") as handle:
        handle.write("\t".join([
            "taxid", "name", "reads", "libraries_seen", "prevalence", "controls",
            "control_cpm", "threshold_cpm", "samples_kept", "samples_zeroed", "verdict",
        ]) + "\n")
        for entry in sorted(evidence, key=lambda e: -to_float(e[2])):
            handle.write("\t".join(str(field) for field in entry) + "\n")

    with open(args.drop_list, "w", encoding="utf-8") as handle:
        for taxid in condemned:
            handle.write(f"{taxid}\n")

    if args.drop_cells:
        with open(args.drop_cells, "w", encoding="utf-8") as handle:
            handle.write("taxid\tsample\n")
            for taxid, sample in zeroed_pairs:
                if taxid:
                    handle.write(f"{taxid}\t{sample}\n")

    if args.mqc:
        with open(args.mqc, "w", encoding="utf-8") as handle:
            handle.write(
                "# id: 'reanatax_control_filter'\n"
                "# section_name: 'Negative-control filter'\n"
                "# description: 'Taxa judged against the level that arrives WITHOUT a sample.\n"
                "#     The verdict is per library, not per taxon: the same organism can be signal\n"
                "#     in one well and carryover in the next, which is what a control library\n"
                "#     measures and what no per-read evidence can see. control_level_somewhere\n"
                "#     means the taxon survived in at least one library and was zeroed in others.'\n"
                "# plot_type: 'bargraph'\n"
                "# pconfig:\n"
                "#     id: 'reanatax_control_filter_plot'\n"
                "#     title: 'reanaTax: negative-control filter'\n"
                "#     ylab: 'Taxa'\n"
            )
            keys = sorted(counts)
            handle.write("Sample\t" + "\t".join(keys) + "\n")
            handle.write("all taxa\t" + "\t".join(str(counts[key]) for key in keys) + "\n")

    summary = ", ".join(f"{count} {call}" for call, count in sorted(counts.items()))
    print(
        f"[control_filter] {len(rows)} taxa: {summary}. "
        + (
            f"{len(controls)} control(s) vs {len(real)} sample(s); zeroed {zeroed_cells} cell(s) "
            f"holding {zeroed_reads:.0f} of {total_reads:.0f} reads "
            f"({100 * zeroed_reads / total_reads if total_reads else 0:.1f}%) "
            f"below {args.ratio:g}x the control {args.statistic} "
            f"(floored at {floor_cpm:.2f} CPM = {args.floor_reads:g} read(s) in a "
            f"{control_depth:.0f}-read control). "
            if controls
            else ""
        )
        + (
            f"--prevalence dropped taxa present in >= {prevalence_floor:.0f}/{n_libraries} libraries. "
            if prevalence_floor is not None
            else ""
        )
        + f"{len(condemned)} taxid(s) listed for removal.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())

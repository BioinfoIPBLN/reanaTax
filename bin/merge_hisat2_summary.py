#!/usr/bin/env python3
"""Sum a set of HISAT2 summary logs into one, as if a single alignment had run.

When a library is aligned in chunks (--hisat2_chunk_size) HISAT2 writes one
summary per chunk. MultiQC would show each of them as a separate sample and
bin/read_accounting.py would read the chunk id as the sample name, so the
chunks have to be folded back into a single file in HISAT2's own format before
anything downstream sees them.

Every count is additive across chunks because a chunk is a disjoint slice of the
same library: `-s`/`-u` partition the reads, they do not resample them. The
percentages are not additive and are recomputed from the summed counts.
"""

import argparse
import re
import sys

# The summary is Bowtie2's, and HISAT2 emits one of two shapes depending on
# whether the input was paired. Each entry is (key, regex); a key that no chunk
# reports stays absent rather than becoming a zero, so a single-end summary is
# never emitted with paired-only lines.
PAIRED_FIELDS = [
    ("total", r"^\s*(\d+)\s+reads;\s+of these:"),
    ("paired", r"^\s*(\d+)\s+\([\d.]+%\)\s+were paired;"),
    ("conc_0", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly 0 times"),
    ("conc_1", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly exactly 1 time"),
    ("conc_multi", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned concordantly >1 times"),
    ("disc_1", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned discordantly 1 time"),
    ("mate_0", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned 0 times\s*$"),
    ("mate_1", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned exactly 1 time\s*$"),
    ("mate_multi", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned >1 times\s*$"),
]

SINGLE_FIELDS = [
    ("total", r"^\s*(\d+)\s+reads;\s+of these:"),
    ("unpaired", r"^\s*(\d+)\s+\([\d.]+%\)\s+were unpaired;"),
    ("un_0", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned 0 times\s*$"),
    ("un_1", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned exactly 1 time\s*$"),
    ("un_multi", r"^\s*(\d+)\s+\([\d.]+%\)\s+aligned >1 times\s*$"),
]


def parse(path):
    with open(path, encoding="utf-8") as handle:
        text = handle.read()

    # "were paired" and "were unpaired" are mutually exclusive, and the
    # mate-level lines of a paired summary use the same wording as the only
    # alignment lines of a single-end one - so the layout has to be decided
    # first, or a paired summary parses as single-end with the wrong numbers.
    paired = bool(re.search(r"were paired;", text))
    fields = PAIRED_FIELDS if paired else SINGLE_FIELDS

    counts = {}
    for key, pattern in fields:
        match = re.search(pattern, text, re.M)
        if match:
            counts[key] = int(match.group(1))
    return paired, counts


def pct(part, whole):
    """HISAT2 prints 0.00% rather than nothing when the denominator is zero,
    which is what an empty chunk produces - a chunk whose offset is past the
    end of the file is legitimate, not an error."""
    return 0.0 if not whole else 100.0 * part / whole


def render_paired(c):
    total = c.get("total", 0)
    paired = c.get("paired", 0)
    conc_0 = c.get("conc_0", 0)
    disc_1 = c.get("disc_1", 0)
    unaligned_pairs = conc_0 - disc_1
    unaligned_mates = unaligned_pairs * 2
    mate_0 = c.get("mate_0", 0)

    # Bowtie2's "overall alignment rate" for paired input is over MATES, not
    # pairs: every mate of every pair, minus the ones that aligned nowhere.
    rate = pct(2 * total - mate_0, 2 * total)

    return "\n".join([
        f"{total} reads; of these:",
        f"  {paired} ({pct(paired, total):.2f}%) were paired; of these:",
        f"    {conc_0} ({pct(conc_0, paired):.2f}%) aligned concordantly 0 times",
        f"    {c.get('conc_1', 0)} ({pct(c.get('conc_1', 0), paired):.2f}%) aligned concordantly exactly 1 time",
        f"    {c.get('conc_multi', 0)} ({pct(c.get('conc_multi', 0), paired):.2f}%) aligned concordantly >1 times",
        "    ----",
        f"    {conc_0} pairs aligned concordantly 0 times; of these:",
        f"      {disc_1} ({pct(disc_1, conc_0):.2f}%) aligned discordantly 1 time",
        "    ----",
        f"    {unaligned_pairs} pairs aligned 0 times concordantly or discordantly; of these:",
        f"      {unaligned_mates} mates make up the pairs; of these:",
        f"        {mate_0} ({pct(mate_0, unaligned_mates):.2f}%) aligned 0 times",
        f"        {c.get('mate_1', 0)} ({pct(c.get('mate_1', 0), unaligned_mates):.2f}%) aligned exactly 1 time",
        f"        {c.get('mate_multi', 0)} ({pct(c.get('mate_multi', 0), unaligned_mates):.2f}%) aligned >1 times",
        f"{rate:.2f}% overall alignment rate",
        "",
    ])


def render_single(c):
    total = c.get("total", 0)
    unpaired = c.get("unpaired", 0)
    un_0 = c.get("un_0", 0)
    rate = pct(total - un_0, total)

    return "\n".join([
        f"{total} reads; of these:",
        f"  {unpaired} ({pct(unpaired, total):.2f}%) were unpaired; of these:",
        f"    {un_0} ({pct(un_0, unpaired):.2f}%) aligned 0 times",
        f"    {c.get('un_1', 0)} ({pct(c.get('un_1', 0), unpaired):.2f}%) aligned exactly 1 time",
        f"    {c.get('un_multi', 0)} ({pct(c.get('un_multi', 0), unpaired):.2f}%) aligned >1 times",
        f"{rate:.2f}% overall alignment rate",
        "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("summaries", nargs="+", help="per-chunk HISAT2 summary logs")
    parser.add_argument("--output", required=True, help="merged summary to write")
    args = parser.parse_args()

    totals = {}
    layouts = set()
    for path in args.summaries:
        paired, counts = parse(path)
        if not counts:
            print(f"[merge_hisat2_summary] {path} holds no recognisable HISAT2 summary",
                  file=sys.stderr)
            continue
        layouts.add(paired)
        for key, value in counts.items():
            totals[key] = totals.get(key, 0) + value

    if not totals:
        sys.exit("[merge_hisat2_summary] none of the given files held a HISAT2 summary; "
                 "refusing to write an empty one rather than reporting a library with no reads")

    # A sample cannot be half paired: if the chunks disagree, one of them was
    # aligned against a different input and summing them would invent a library
    # that never existed.
    if len(layouts) > 1:
        sys.exit("[merge_hisat2_summary] chunk summaries disagree on whether the input "
                 "was paired; they cannot be from the same library")

    text = render_paired(totals) if layouts.pop() else render_single(totals)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(text)


if __name__ == "__main__":
    main()

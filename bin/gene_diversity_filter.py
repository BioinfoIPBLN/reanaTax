#!/usr/bin/env python3
"""
gene_diversity_filter.py -- does a taxon's evidence spread over its genome?

The idea is PRISM's. Every filter upstream of this one counts reads, or counts
the k-mers behind reads; none of them can see WHERE on the organism those reads
landed. That distinction is the whole difference between a species that is
present and a species whose name has been attached to one conserved stretch of
sequence. A thousand reads on a 16S gene, on a ribosomal protein, on a
transposase shared across half a phylum, is a thousand reads - and every
abundance threshold in this pipeline passes it.

A genuine organism transcribes hundreds of genes. Its reads fall on many
distinct gene families and no single family carries the bulk of them. A
spurious call has one locus and nothing else, because there was only ever one
region of homology.

This does not need a new reference or a new alignment. HUMAnN already produces
exactly the required table: gene-family abundances STRATIFIED BY SPECIES, so
every `UniRef90_xxx|g__Genus.s__Species` row says how much of that gene family
was attributed to that organism. Two quantities per taxon fall straight out:

  gene_families   how many distinct families the taxon's reads reached at all.
  top_fraction    what share of the taxon's total abundance its single largest
                  family carries. A taxon at 0.95 is one locus wearing a
                  species name, however many reads it has.

Run the same table regrouped onto KEGG orthologs (--humann_regroup uniref90_ko)
and the same numbers describe PRODUCTS rather than gene families - a taxon can
reach many UniRef90 families that all encode the same thing, and the product
view is the one that catches it.

PRISM's own definitions, from sjdlabgroup/PRISM functions.R, are reproduced
where the input allows, and the columns carry PRISM's names so the two can be
lined up:

    prism_genbank_stats = function(prod){
      prod %>% group_by(staxids) %>%
        summarize(fprod = n(), fugene = length(unique(gene)),
                  fuprod = length(unique(product)),
                  prod_div = product %>% table() %>% vegan::diversity() %>% mean(),
                  gene_div = gene %>% table() %>% vegan::diversity() %>% mean(),
                  .groups = 'drop') %>%
        mutate(fprod = fprod/sum(fprod), fugene = fugene/sum(fugene),
               fuprod = fuprod/sum(fuprod))
    }

So `gene_div` is RAW SHANNON ENTROPY - `vegan::diversity()`'s default index -
not an evenness on 0..1, and `fprod`/`fugene`/`fuprod` are shares of the
across-taxa total rather than absolute counts. Both are matched here.

Three things differ, and they are differences of INPUT, not of formula:

  provenance   PRISM builds its `prod` table from BLAST hits intersected with
               GenBank annotation, so one row is one READ on one annotated
               feature. Here a row is one HUMAnN gene family with its abundance
               for that species. Same question, different instrument: PRISM
               pays for a BLAST pass against nt, this reads a table HUMAnN has
               already written.
  units        PRISM counts reads; HUMAnN reports RPK-derived abundances. The
               entropy of an abundance vector is the natural analogue, but it
               is not the entropy of a read count and the two will not agree to
               the decimal.
  scope        PRISM normalises within one sample. These are pooled over the
               cohort, like every other evidence filter here, because whether a
               taxon's reads have breadth is a property of the taxon and a
               single shallow library is a poor place to settle it.

What is NOT reproduced is PRISM's score. That is an XGBoost model
(`prismxg.RDS`) over FORTY features, of which these five are five; the other
thirty-five are BLAST multi-mapping ratios, k-mer taxonomy proportions at seven
ranks, and Kraken lineage/misclassification statistics, none of which exist
without PRISM's own BLAST and Kraken passes. The `verdict` column below is
therefore this pipeline's own threshold rule and is labelled as such - it is
NOT a PRISM call. Run PRISM itself (--run_prism) for that.

Like every other evidence filter here, nothing is removed. The taxids go to the
abundance filter so one step owns every removal.
"""
import argparse
import csv
import gzip
import math
import re
import sys

SPECIES = re.compile(r"\|.*s__(?P<species>[^|]+)$")
GENUS_ONLY = re.compile(r"\|g__(?P<genus>[^.|]+)$")
UNWANTED = ("UNMAPPED", "UNINTEGRATED", "UNGROUPED")


def open_maybe_gzip(path):
    return gzip.open(path, "rt", errors="replace") if str(path).endswith(".gz") else open(
        path, encoding="utf-8", errors="replace"
    )


def humann_name_to_taxon(label):
    """`g__Escherichia.s__Escherichia_coli` -> `Escherichia coli`, or None."""
    match = SPECIES.search(label)
    if match:
        return match.group("species").replace("_", " ").strip()
    match = GENUS_ONLY.search(label)
    if match:
        return match.group("genus").replace("_", " ").strip()
    return None


def read_stratified(paths):
    """{taxon_name: {gene_family: abundance}} pooled over every sample given."""
    per_taxon = {}
    for path in paths:
        with open_maybe_gzip(path) as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, None)
            if header is None:
                continue
            for row in reader:
                if not row or not row[0]:
                    continue
                label = row[0].strip()
                if label.startswith("#") or label.split("|")[0] in UNWANTED:
                    continue
                if "|" not in label:
                    # The unstratified total for a gene family. Its per-species
                    # rows follow, so counting it too would double everything.
                    continue
                taxon = humann_name_to_taxon(label)
                if not taxon:
                    continue
                family = label.split("|")[0]
                total = 0.0
                for value in row[1:]:
                    try:
                        total += float(value)
                    except ValueError:
                        continue
                if total <= 0:
                    continue
                bucket = per_taxon.setdefault(taxon, {})
                bucket[family] = bucket.get(family, 0.0) + total
    return per_taxon


def read_taxid_map(path):
    """{lowercased taxon name: taxid} from the combined Bracken table."""
    mapping = {}
    if not path:
        return mapping
    with open_maybe_gzip(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            return mapping
        lookup = {name.strip().lower(): name for name in reader.fieldnames}
        name_key = lookup.get("name")
        taxid_key = lookup.get("taxonomy_id") or lookup.get("taxid")
        if not name_key or not taxid_key:
            sys.exit(
                f"[gene_diversity] '{path}' has no name/taxonomy_id columns, so the taxa this "
                "finds could not be turned into a drop list. Pass the combined Bracken table."
            )
        for row in reader:
            name = (row.get(name_key) or "").strip()
            taxid = (row.get(taxid_key) or "").strip()
            if name and taxid:
                mapping[name.lower()] = taxid
    return mapping


def shannon(abundances):
    """Raw Shannon entropy, natural log - what `vegan::diversity()` returns.

    Reported unscaled because that is the quantity PRISM's model was trained
    on. It is bounded above by log(number of features), so it rises with the
    number of features as well as with how evenly they are used; `evenness()`
    below separates the two.
    """
    total = sum(abundances)
    if total <= 0:
        return 0.0
    entropy = 0.0
    for value in abundances:
        if value <= 0:
            continue
        share = value / total
        entropy -= share * math.log(share)
    return entropy


def evenness(abundances):
    """Shannon evenness on 0..1. One family gives 0; a flat spread gives 1.

    Not a PRISM quantity, and reported alongside rather than instead of the
    entropy: a taxon on two features used equally and a taxon on two hundred
    used equally are both perfectly even, and only the entropy tells them
    apart.
    """
    if len(abundances) < 2:
        return 0.0
    return shannon(abundances) / math.log(len(abundances))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genefamilies", nargs="+", required=True,
                        help="HUMAnN gene-family tables, stratified by species")
    parser.add_argument("--bracken", default=None,
                        help="combined Bracken table, for the name -> taxid map")
    parser.add_argument("--min-genes", type=int, default=10,
                        help="a taxon must reach this many distinct families")
    parser.add_argument("--max-top-fraction", type=float, default=0.5,
                        help="...and its largest family must carry no more than this share")
    parser.add_argument("--min-abundance", type=float, default=0.0,
                        help="below this pooled abundance a taxon is not judged")
    parser.add_argument("--label", default="gene", help="'gene' or 'product'; names the columns")
    parser.add_argument("--prefix", default="reanatax")
    args = parser.parse_args()

    per_taxon = read_stratified(args.genefamilies)
    if not per_taxon:
        sys.exit(
            "[gene_diversity] no species-stratified row in any input. HUMAnN writes gene families "
            "as `UniRef90_xxx|g__Genus.s__Species`; a table with only unstratified rows carries no "
            "information about which organism a gene family came from, and this cannot run on it."
        )

    taxids = read_taxid_map(args.bracken)

    # PRISM's `fprod`, `fugene` and `fuprod` are shares of the across-taxa
    # total, not absolute counts, so the denominators are needed before any row
    # can be written.
    total_abundance = sum(sum(families.values()) for families in per_taxon.values())
    total_features = sum(len(families) for families in per_taxon.values())

    # PRISM names the unique-feature column for the KIND of feature, so the
    # column a reader looks for depends on which table was handed in.
    unique_column = "fugene" if args.label == "gene" else "fuprod"
    div_column = "gene_div" if args.label == "gene" else "prod_div"

    rows = []
    for taxon, families in sorted(per_taxon.items(), key=lambda item: -sum(item[1].values())):
        abundances = sorted(families.values(), reverse=True)
        total = sum(abundances)
        if total < args.min_abundance:
            verdict = "untested_low_abundance"
        elif len(abundances) < args.min_genes or abundances[0] / total > args.max_top_fraction:
            verdict = "single_locus"
        else:
            verdict = "clean"
        rows.append({
            "taxid": taxids.get(taxon.lower(), ""),
            "name": taxon,
            "features": len(abundances),
            "abundance": round(total, 6),
            "fprod": round(total / total_abundance, 6) if total_abundance else 0.0,
            unique_column: round(len(abundances) / total_features, 6) if total_features else 0.0,
            div_column: round(shannon(abundances), 4),
            "evenness": round(evenness(abundances), 4),
            "top_feature_fraction": round(abundances[0] / total, 4) if total else 0.0,
            "verdict": verdict,
        })

    columns = ["taxid", "name", "features", "abundance", "fprod", unique_column,
               div_column, "evenness", "top_feature_fraction", "verdict"]
    with open(f"{args.prefix}.{args.label}_diversity.tsv", "w", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in columns) + "\n")

    flagged = [row for row in rows if row["verdict"] == "single_locus"]
    with_taxid = [row for row in flagged if row["taxid"]]
    with open(f"{args.prefix}.{args.label}_diversity_drop.txt", "w", encoding="utf-8") as handle:
        for row in with_taxid:
            handle.write(f"{row['taxid']}\n")
    if len(flagged) > len(with_taxid):
        # Said out loud rather than dropped silently: HUMAnN names a species by
        # its MetaPhlAn clade name, and a taxon Bracken never called has no
        # taxid to put on the list.
        print(
            f"[gene_diversity] {len(flagged) - len(with_taxid)} flagged taxa have no taxid in the "
            "Bracken table and could not be put on the drop list: "
            + ", ".join(row["name"] for row in flagged if not row["taxid"])[:400],
            file=sys.stderr,
        )

    untested = sum(1 for row in rows if row["verdict"].startswith("untested"))
    with open(f"{args.prefix}_{args.label}_diversity_mqc.tsv", "w", encoding="utf-8") as handle:
        handle.write("\n".join([
            f"# id: 'reanatax_{args.label}_diversity'",
            f"# section_name: 'Breadth of evidence across {args.label} families'",
            "# description: 'Whether a taxon's reads spread over its genome or pile onto one",
            "#     locus, read from HUMAnN's species-stratified gene-family table. A thousand",
            "#     reads on a single conserved gene is a thousand reads, and every abundance",
            "#     threshold passes it; only breadth distinguishes an organism that is present",
            "#     from a name attached to one stretch of homologous sequence. Taxa are flagged",
            f"#     below {args.min_genes} distinct families, or when one family carries more than",
            f"#     {args.max_top_fraction:g} of the taxon''s abundance. The evidence table also",
            f"#     carries PRISM''s own {div_column}/fprod/{unique_column} columns, but the",
            "#     verdict here is a threshold rule, not PRISM''s model score.'",
            "# plot_type: 'bargraph'",
            "# pconfig:",
            f"#     id: 'reanatax_{args.label}_diversity_plot'",
            f"#     title: 'reanaTax: {args.label}-family breadth'",
            "#     ylab: 'Taxa'",
            "Sample\tSingle locus\tClean\tUntested",
            f"all taxa\t{len(flagged)}\t{len(rows) - len(flagged) - untested}\t{untested}",
            "",
        ]))

    print(
        f"[gene_diversity] {len(rows)} taxa; {len(flagged)} carried by fewer than "
        f"{args.min_genes} {args.label} families or concentrated above "
        f"{args.max_top_fraction:g} on one; {len(with_taxid)} taxid(s) listed for removal.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

# BioinfoIPBLN/reanatax: Output

## Introduction

This document describes the output produced by the pipeline. All paths are relative to `--outdir`.

## Pipeline overview

- [Download](#download) — raw reads and archive metadata fetched from ENA/SRA
- [Reference](#reference) — host genome and HISAT2 index
- [FastQC](#fastqc) — read quality before and after trimming
- [fastp](#fastp) — adapter and quality trimming
- [HISAT2](#hisat2) — host alignment, host BAMs and non-host FASTQs
- [Kraken2](#kraken2) — taxonomic classification
- [Bracken](#bracken) — abundance re-estimation
- [Krona](#krona) — interactive taxonomy charts
- [Assembly](#assembly) — contigs from the non-host fraction, and what they say about each taxon
- [BIOM](#biom) — the same counts in BIOM format, for tools outside this pipeline
- [MultiQC](#multiqc) — aggregate report
- [Pipeline information](#pipeline-information) — run metadata, versions and reports

### Download

Only produced when `--input_accessions` is used.

<details markdown="1">
<summary>Output files</summary>

- `download/metadata/`
  - `<accession>-run-info.tsv`: every field ENA/SRA holds for each run under the queried accession.
  - `<accession>.runsheet.csv`: the slim, normalised run table the pipeline actually reads — one row per run with its experiment, sample and study accession, library layout and title.
- `download/fastq/`
  - `<run>-run-info.tsv`: per-run download record, including which archive served the file.

</details>

Disable with `--save_download_metadata false`. The FASTQ files themselves are not published by default — they are intermediates, and the run tables plus the accessions are enough to reproduce them exactly.

### Reference

<details markdown="1">
<summary>Output files</summary>

- `reference/genome/`
  - `*_genomic.fna.gz`: the genome as downloaded from NCBI (only when `--host_accession`/`--host_taxid` was used).
  - `*.fna` / `*.fasta`: the decompressed genome that was indexed.
- `reference/hisat2/`
  - `*.ht2`: the HISAT2 index.

</details>

Published by default so later runs can skip the build with `--hisat2_index <outdir>/reference/hisat2`. Disable with `--save_reference false`.

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/raw/`
  - `<sample>_raw*_fastqc.html`: report for the reads as they arrived.
- `fastqc/trimmed/`
  - `<sample>_trimmed*_fastqc.html`: report for the reads after fastp.

</details>

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives per-base quality, adapter content, duplication and GC. The raw/trimmed pair is what tells you whether trimming did what you expected. The `.zip` archives are not published — their contents are in the MultiQC report.

### fastp

<details markdown="1">
<summary>Output files</summary>

- `fastp/`
  - `<sample>.fastp.html`, `<sample>.fastp.json`: trimming report.
  - `<sample>.fastp.log`: tool log.
- `fastp/trimmed/` (with `--save_trimmed`)
  - `<sample>*.fastp.fastq.gz`: the trimmed reads.
- `fastp/failed/` (with `--save_trimmed_fail`)
  - `<sample>*.fail.fastq.gz`: reads fastp discarded.

</details>

[fastp](https://github.com/OpenGene/fastp) does adapter detection and quality trimming in one pass. The JSON feeds MultiQC.

### HISAT2

<details markdown="1">
<summary>Output files</summary>

- `hisat2/`
  - `<sample>.host.sorted.bam`, `.bai`: reads that aligned to the host, sorted and indexed.
- `hisat2/log/`
  - `<sample>.host.hisat2.summary.log`: alignment summary, including the overall alignment rate.
- `hisat2/samtools_stats/`
  - `*.stats`, `*.flagstat`, `*.idxstats`: alignment statistics.
- `unaligned/`
  - `<sample>.host.unmapped_1.fastq.gz`, `<sample>.host.unmapped_2.fastq.gz` (paired-end)
  - `<sample>.host.unmapped.fastq.gz` (single-end)

</details>

The **unaligned FASTQs are the pipeline's main intermediate product** — they are the non-host fraction that Kraken2 classifies, and they are also the right input for any downstream assembly or targeted analysis you want to run yourself.

For paired-end data these are HISAT2's `--un-conc-gz` output: pairs that did not align *concordantly*. A pair where only one mate hit the host is therefore kept, which is the conservative choice for depletion.

The overall alignment rate in the summary log is the number to look at first — it is the host fraction of the library.

Turn either output off with `--save_host_bam false` / `--save_unaligned false`.

### Qualimap

<details markdown="1">
<summary>Output files</summary>

- `qualimap/<sample>/`
  - `qualimapReport.html`: the BamQC report for that sample's host BAM.
  - `genome_results.txt`: the same numbers as plain text.
  - `raw_data_qualimapReport/`: the table behind each plot.

</details>

[Qualimap](http://qualimap.conesalab.org/) BamQC covers what the samtools statistics above cannot: coverage depth and its uniformity across the reference, duplication rate, GC content of the mapped reads, and the mapping-quality and insert-size distributions. Its headline numbers also appear in MultiQC.

Read it with the host-depletion context in mind. A low genome fraction and shallow, patchy coverage are the *expected* result for a library that is not dominated by host — they are only a warning sign when paired with poor mapping quality, which usually means the host reference is the wrong one.

`--skip_qualimap` turns the step off; `--qualimap_gff` adds feature-level statistics.

### Host gene counts

<details markdown="1">
<summary>Output files</summary>

- `host_counts/`
  - `<sample>.host.featureCounts.tsv`: per-gene counts for the host alignment.
  - `<sample>.host.featureCounts.tsv.summary`: assignment summary, also in MultiQC.

</details>

Only written with `--quantify_host` and `--gtf`. These come from the same library as the microbial profile, which is what makes correlating the two defensible — see [usage.md](usage.md).

### Kraken2

<details markdown="1">
<summary>Output files</summary>

- `kraken2/`
  - `<sample>.kraken2.report.txt`: per-sample classification report.
  - `kraken2_combined_report.txt`: all samples in one table.
- `kraken2/reads/` (with `--kraken2_save_reads`)
  - `<sample>.classified*.fastq.gz`, `<sample>.unclassified*.fastq.gz`
- `kraken2/read_assignments/` (with `--kraken2_save_readclassifications`)
  - `<sample>.kraken2.classifiedreads.txt`: one line per read. These files are large.

</details>

[Kraken2](https://ccb.jhu.edu/software/kraken2/) assigns each read to a taxon by exact k-mer matching. The report columns are: percentage of reads in the clade, reads in the clade, reads assigned directly to the taxon, rank code, NCBI taxonomy ID, and name.

Read the percentages as *read* abundance, not organism abundance — that is what Bracken corrects.

A sample in which Kraken2 classified nothing (its report holds only the `unclassified` row) is still published and still reaches MultiQC, but is excluded from Bracken and from the combined tables, with a warning naming the sample. Both `combine_kreports.py` and Bracken fail outright on such reports, and one empty sample should not take down the run.

### Evidence filters

<details markdown="1">
<summary>Output files</summary>

- `minimizer_filter/` (with `--minimizer_filter`)
  - `*.minimizer_evidence.tsv`: every taxon with its reads, distinct minimizers, duplication, coverage and verdict. With `--minimizer_correlation` it also carries the three Spearman coefficients and the largest of their BH-adjusted p-values.
  - `*.minimizer_drop.txt`: the taxids that failed.
- `host_kmer_filter/` (with `--host_kmer_filter`)
  - `*.host_kmer_evidence.tsv`: per taxon, its reads, how many of them carried host k-mers, and the resulting fraction.
  - `*.host_kmer_drop.txt`: the taxids that failed.

</details>

Neither filter removes anything itself. Both write a taxid list that the abundance filter consumes, so one step owns every removal from the combined tables and one `.removed.tsv` records them all — and a taxon condemned by either is removed, since surviving one test is no argument against the other.

- `cleanup/` (with `--cleanup_intermediates`)
  - `<sample>.cleanup.log`: one line per intermediate, saying whether it was removed and how large it was, or kept because it lay outside the work directory. The kept lines are the audit trail for a samplesheet's own FASTQs, which are never touched.

- `decontam/` (with `--decontam`)
  - `reanatax.decontam_evidence.tsv`: every taxon with its decontam score, the prevalence/frequency sub-scores, and whether it was called a contaminant.
  - `reanatax.decontam_drop.txt`: the taxids that were.

`--decontam` is the only one of these given an external measurement of the kit, and so the only one that can separate a reagent contaminant from a genuinely rare organism. With `--decontam_batch_column` the evidence table also carries `n_batches_flagged` and `batches_flagged`: under the default `minimum` rule a taxon condemned in one batch of six and one condemned in all six get the same verdict, and they are not the same claim.

- `control_filter/` (with `--negative_controls` or `--prevalence_filter`)
  - `kraken2_combined_report.control_evidence.tsv` and `bracken_combined_<level>.control_evidence.tsv`: per taxon, its total reads, how many libraries carry it and at what prevalence, the control level and the threshold that follows from it (both counts per million classified reads), how many libraries it survived in and how many it was zeroed in, and a verdict — `kept`, `control_level_somewhere`, `control_level_everywhere` or `ubiquitous`.
  - `*.control_drop.txt`: the taxids condemned outright, being those that failed in *every* library plus those `--prevalence_filter` removed.
  - `*.control_cells.tsv`: `taxid<TAB>sample` for every cell that was zeroed. This is the per-library half of the verdict, written out because objects built from the per-read assignments rather than from the combined tables — the cell-by-taxon matrix above all — cannot see the zeroing otherwise. `--sc_apply_drop_list` feeds it to the single-cell filter, so the cohort profile and the single-cell profile cannot disagree about a cell one of them removed. The kreport's list is the one used there: Bracken's table is species-only, and the matrix carries whatever rank Kraken2 assigned.

This is the only step whose verdict is **per library rather than per taxon**, and the only one that addresses index hopping and well-to-well carryover at all: those put genuine reads of a genuine organism into the wrong library, so no evidence filter can find fault with the reads, and the same taxon is then signal in one library and carryover in the next. It therefore rewrites the table — the failing cells are zeroed — rather than only naming taxids, and it runs *before* the abundance filter so the surviving fractions are recomputed against what is left. `control_level_somewhere` is the verdict a per-taxon filter cannot produce.

The rewritten tables themselves are not published: they are intermediates that the abundance filter reads, and two tables in one directory differing only in which cells are zero would be a trap. `kraken2/` and `bracken/` hold the published before-and-after as usual, and the evidence table says which cells went.

Exact on the Bracken table, which is flat and holds each read once. On the `combine_kreports` hierarchy a zeroed clade leaves its ancestors' clade counts stale — the same limitation `--drop_host_taxon` carries, with the same answer: take Bracken downstream.

The control libraries are never filtered — a control cannot be judged against itself — and stay in the table as columns. Exclude them downstream, or with `--da_samples_to_drop`.

- `shuffle_control/` (with `--shuffle_control`)
  - `reanatax.shuffle_evidence.tsv`: per taxon, its real reads, the reads it collected from the shuffled copy, that count scaled to the real library size, the ratio, and a verdict — `clean`, `composition_only`, or one of the two `untested_*` states.
  - `reanatax.shuffle_drop.txt`: the taxids reported from shuffled sequence.
  - `kraken2/`: the shuffled classification itself, named after the same samples but kept out of MultiQC so it is not read as a second cohort.
  - `reads/<sample>.shuffle_stats.tsv`: how many fragments were shuffled, by which method, at what GC. The shuffled FASTQs themselves are not published — they are the size of the library and are a means, not a result.

- `gene_diversity/` (with `--gene_diversity_filter`)
  - `reanatax.gene_diversity.tsv` and, with `--humann_regroup`, `reanatax.product_diversity.tsv`: per taxon, how many distinct gene families (or products) its reads reached, what share of its abundance the largest one carries, and the Shannon evenness across them.
  - `reanatax.*_diversity_drop.txt`: the taxids with no breadth.
  - `bracken_combined_<level>_genediv.tsv`: the Bracken table with those taxa removed. It lives here rather than in `bracken/` because the judgement was HUMAnN's, made downstream of everything in that directory.

Each filter asks a different question of the same reports. `--minimizer_filter` asks whether a taxon's reads cover enough of its reference; `--host_kmer_filter` asks whether they are host sequence the aligner missed; `--shuffle_control` asks how much of the signal the database would have produced from composition alone; `--gene_diversity_filter` asks whether the reads spread over the genome or pile onto one locus; `--negative_controls` asks whether there is more of it here than arrives without a sample. Surviving one is no argument against the others. See [usage](usage.md#filtering-on-evidence-not-abundance---minimizer_filter).

The first four are all blind to carryover and index hopping, and measurably so: on the CSI-Microbes plate the distinct-minimizers-per-read ratio ranks the carryover calls *above* the true positives (AUC 0.415). Hopped reads are the same reads that the neighbouring well produced correctly, so there is nothing in them for an evidence filter to find. Only `--negative_controls` speaks to that.

`--shuffle_control` is blind to host carry-over — carry-over reads are real sequence, and shuffling removes them — so a carried-over taxon passes it perfectly. Read it beside `host_kmer_filter/`, never instead of it.

### Single cell

<details markdown="1">
<summary>Output files</summary>

- `starsolo/` (with `--single_cell`)
  - `<sample>_Solo.out/`: STARsolo's cell-by-gene matrices — raw and filtered, plus its barcode/UMI statistics. This is the host half of the experiment.
  - `<sample>.barcode_stats.tsv`: how many non-host reads carried a corrected barcode, and how many distinct barcodes they came from.
  - `log/<sample>.Log.final.out`: STAR's alignment summary, also in MultiQC.
- `cell_taxa/`
  - `<sample>.cell_taxa.tsv`: the cell-by-taxon matrix, long format — `sample, barcode, taxid, rank, name, count`.
  - `<sample>.cell_taxa_summary.tsv`: how many reads were counted and why each dropped read was dropped (no barcode, host k-mer, below min-frac, low complexity, wrong rank).
  - `<sample>.sc_kmer_evidence.tsv`, `<sample>.sc_kmer_drop.txt` (with `--sc_kmer_denoise`): per taxon, the Spearman correlation between its total and distinct k-mers across barcodes, and the verdict. `saturated_distinct` means the distinct count never moved while the read count did.

</details>

  - `reanatax.cell_type_enrichment.tsv` (with `--sc_cell_metadata`): per taxon and cell type — observed and expected infected cells, log2FC, Stouffer-combined p and BH-adjusted q.
  - `reanatax.cooccurrence.tsv`, `reanatax.doublet_check.tsv`: taxon pairs found together in single cells more often than chance, and the doublet control that says whether to believe it.
  - `<sample>.cell_taxa_sweep.tsv`: the same matrix at `--sc_min_umis` 1 through 5. Check your conclusion survives the sweep before reporting it.

Join `cell_taxa.tsv` to the Solo matrix **on the barcode** — that is the whole point of the branch, and the two are guaranteed to use the same barcode vocabulary because both come from STARsolo's corrected `CB`.

The barcode-level verdicts are **per sample and are not applied** to the combined tables — barcodes only mean anything within their own library. They are published for you to apply.

- `singlecell/ambient/` (with `--sc_ambient`)
  - `<sample>.sc_ambient.tsv`: per taxon, how many cells and how many empty droplets carry it, the UMIs in each pool, the prevalence and rate ratios, and a verdict. `ambient_undecided` is **not** a contamination call — it means the taxon is as common in the medium as in the cells, which a reagent contaminant and a genuine extracellular organism both look like. Cross-reference `decontam/` to separate them.
  - `<sample>.sc_ambient_drop.txt` (with `--sc_ambient_drop`): the ambient and cell-depleted taxa. Only meaningful if the question is about intracellular microbes.
- `singlecell/host_de/` (with `--sc_host_de`)
  - `<sample>.sc_host_de.tsv`: per taxon, cell type and gene — detection rates and mean expression on both sides, log2 fold change, Wilcoxon p and the q adjusted within that taxon x cell type.
  - `<sample>.sc_host_de_groups.tsv`: which taxon x cell-type groups were testable at all, with the number of infected and bystander cells behind each. Read this first: a group missing from it had too few cells on one side, which is a different statement from finding nothing.
- `singlecell/reanatax.cell_taxa.tsv` (with `--sc_plate_based`): the same long-format matrix, built from the per-cell Kraken2 reports instead of from barcodes. `barcode` holds the cell id and `sample` holds the patient.

A `cell_type` of `ALL_POOLED` in `sc_host_de.tsv` means `--sc_host_de_force_pooled` was used and the comparison is **not** within cell type — the difference between the cell types that carry the taxon is inseparable from the response to carrying it.

Read `cell_taxa_summary.tsv` before the matrix. A large `no_barcode` count means the chemistry or whitelist is wrong; a large `host_kmer` count means real carry-over; a large `below_min_frac` count means reads are scattering across the taxonomy rather than landing coherently.

### PRISM

<details markdown="1">
<summary>Output files</summary>

- `prism/<sample>/prism_out/` (with `--run_prism`)
  - `<sample>-counts.csv`: PRISM's per-species summary — reads confirmed and the model's score.
  - `<sample>-results.csv`: the per-read table after all filtering and scoring.
  - `<sample>_1.fa`, `<sample>_2.fa`: the retained microbial reads.
  - `data/<sample>-xgmat.csv`: the full forty-feature matrix the model was given. This is the file to read when a score surprises you.
- `prism/`
  - `reanatax.prism_reads.tsv`, `reanatax.prism_score.tsv`: taxa by sample. A taxon absent from a sample is `0` in the reads matrix and **`NA`** in the score matrix — a zero score would read as "PRISM was confident this is a contaminant", which is the opposite of "PRISM never saw it".
  - `reanatax.prism_evidence.tsv`: per taxon, how many samples saw it, its read totals, and the median/min/max score behind the verdict.
  - `reanatax.prism_drop.txt`, `bracken_combined_<level>_prism.tsv` (with `--prism_filter`): the taxa called contaminants, and the Bracken table without them.

</details>

`verdict` has four values, and `untested_low_depth` is not a pass: below `--prism_min_reads` there is nothing to confirm either way, which is PRISM's own position.

### PathSeq

<details markdown="1">
<summary>Output files</summary>

- `pathseq/` (with `--pathseq_microbe_bwa_image`)
  - `<sample>.pathseq.scores.txt`: PathSeq's own per-sample table.
  - `<sample>.filter_metrics.txt`, `<sample>.score_metrics.txt`: how many reads survived each stage.
  - `pathseq_reads.tsv`, `pathseq_unambiguous.tsv`, `pathseq_score.tsv`, `pathseq_score_normalized.tsv`: taxa by sample, one matrix per quantity.
  - `pathseq_lineage.tsv`: the taxonomy string for each taxon.
  - `bam/` (with `--pathseq_save_bam`): every read PathSeq aligned, tagged with its call.

</details>

`unambiguous` is the column to set beside a Kraken2 species count: it holds the reads that aligned to that taxon and nowhere else. `reads` includes reads shared with other taxa, because PathSeq divides an ambiguous read between the taxa it hits rather than pushing it to their common ancestor — so the `reads` column does **not** sum to the library. `score` is length-normalised and is comparable between taxa within a sample in a way `reads` is not.

Nothing here is converted into a Kraken report. Where PathSeq and Kraken2 disagree about a taxon, that disagreement is the interesting part and is left visible.

### Co-occurrence

<details markdown="1">
<summary>Output files</summary>

- `sparcc/` (with `--run_sparcc`)
  - `reanatax.sparcc_otu.tsv`, `reanatax.sparcc_taxa.tsv`: what went into the network, after the prevalence filter, and the id-to-name map.
  - `reanatax.sparcc_correlation.tsv`, `reanatax.sparcc_covariance.tsv`, `reanatax.sparcc_pvalues.tsv`: FastSpar's square matrices.
  - `reanatax.sparcc_edges_all.tsv`: every pair as a row, with rho, the bootstrap p, and the BH q adjusted across all pairs at once.
  - `reanatax.sparcc_edges_significant.tsv`: the subset clearing both `--sparcc_p_threshold` and `--sparcc_min_correlation`.

</details>

Use the adjusted `q`, not `p`. A 200-taxon network is 19,900 pairs, so at an unadjusted 0.05 roughly a thousand edges are expected from nothing at all — which is how co-occurrence networks acquired their reputation.

### Diversity

<details markdown="1">
<summary>Output files</summary>

- `diversity/` (with `--run_diversity`)
  - `reanatax.alpha_diversity.tsv`: per sample — reads, observed taxa, Shannon, Simpson, inverse Simpson, Pielou's evenness.
  - `reanatax.beta_variance.tsv`: per metadata variable — adjusted R², F, permutation p, BH-adjusted p, PERMANOVA R² and p, and the cumulative non-redundant R² from forward selection.
  - `reanatax.aitchison_distance.tsv`, `reanatax.ordination.tsv`: the sample-by-sample distance matrix and PCoA coordinates.
  - `reanatax.alpha_diversity.png`, `reanatax.ordination.png`.

</details>

Read `beta_variance.tsv` **before** any per-taxon test. A variable that explains more of the community variation than the one you care about is either a confounder or the actual story. The gap between a variable's `r2_adj` and its `cumulative_r2_adj` is how much of its apparent effect is shared with variables already in the model. `r2_adj` is adjusted, so it is comparable between variables with different numbers of levels and goes negative when a variable explains less than chance.

See [usage](usage.md#diversity-and-what-explains-it---run_diversity) for the caveats on small cohorts and on the CLR zero replacement.

### Bracken

<details markdown="1">
<summary>Output files</summary>

- `bracken/`
  - `<sample>.bracken_<level>.tsv`: re-estimated abundances for that sample.
  - `<sample>.bracken_<level>.kraken2.report_bracken.txt`: the same estimate in Kraken2 report format.
  - `bracken_combined_<level>.txt`: all samples in one table — this is the file to take into R or Python.

</details>

[Bracken](https://github.com/jenniferlu717/Bracken) redistributes reads that Kraken2 could only place at a higher rank down to the requested level (`--bracken_level`, default species). **Use the Bracken table, not the raw Kraken2 report, for abundance comparisons between samples.**

### Krona

<details markdown="1">
<summary>Output files</summary>

- `krona/`
  - `<sample>.krona.html`: self-contained interactive chart.

</details>

[Krona](https://github.com/marbl/Krona) renders the composition as a zoomable hierarchy. Built from the Bracken-corrected report when Bracken ran, otherwise from the Kraken2 report. Open the HTML directly in a browser — no server needed.

### Host clade

Only produced with `--drop_host_clade`.

<details markdown="1">
<summary>Output files</summary>

- `host_clade/`
  - `reanatax.host_clade_evidence.tsv`: every taxon removed as the host under another name — taxid, name, rank, reads, and the clade anchor.
  - `reanatax.host_clade_drop.txt`: the taxids, as handed to the abundance filter.

</details>

A host read that misses the host genome is assigned to the nearest relative the database holds, so `--drop_host_taxon` alone leaves most of the leakage behind. Read the evidence table before trusting a wide rank: it is the list of what was deleted, and at `class` or `phylum` it is long.

### Assembly

Only produced with `--assembly`.

<details markdown="1">
<summary>Output files</summary>

- `assembly/`
  - `pool_<key>.contigs.fa.gz`: the contigs, one file per pool (`pool_all`, `pool_<group>` or `pool_<sample>`).
  - `pool_<key>.contigs.kraken2.report.txt`: Kraken2's classification of those contigs.
  - `pool_<key>.megahit.log` / `.spades.log`: the assembler's own log.
  - `reanatax.contig_evidence.tsv`: one row per taxon — reads, pools, contig count, total bp, longest contig, N50, verdict.
  - `reanatax.contig_drop.txt`: the taxids condemned, written whether or not `--contig_filter` is set, and read by the abundance filter only when it is.

</details>

The one line of evidence the other filters are blind to. They all count things about a taxon's reads; a contig is a longer sequence, and reads placed by index hopping or piled on a single conserved locus do not assemble into one.

Contigs roll up the taxonomy, so a contig assigned to a species also supports every rank above it — the lineage comes from the contig report's own indentation. Read `verdict` with the read gate in mind: `not_judged` means the taxon never had enough reads to assemble anything, and is a statement about depth rather than about the taxon.

Nothing is removed unless `--contig_filter` is set; without it the drop list is published as a record and the combined tables are untouched.

### BIOM

Only produced with `--export_biom`.

<details markdown="1">
<summary>Output files</summary>

- `biom/`
  - `bracken_biom.biom` (or `kraken2_biom.biom` with `--skip_bracken`): the classification as a BIOM 1.0 (JSON) table, observation ids are NCBI taxids, observation metadata is the full `k__…;s__…` lineage.
  - `*.unfiltered.biom`: the same table before the filters were applied, written only when they removed something.

</details>

Nothing in this pipeline reads it. The combined tables are already a taxa-by-sample count matrix — which is why the differential-abundance step skips the `kraken-biom` round trip `EzAppDiffShot` does — so this exists to hand a run to QIIME 2, phyloseq, `microbiome` or anything else that speaks BIOM.

Built by [kraken-biom](https://github.com/smdabdoub/kraken-biom) from the same per-sample reports Krona uses, then **restricted to the taxa that survived the filters**, so it agrees with `*_filtered.tsv` rather than contradicting it. The unfiltered table is published beside it when the two differ, because a BIOM handed to another tool is exactly where an unexplained discrepancy is hardest to chase down.

Every rank from domain to species is present as its own observation (`--max D --min S`), not just species. Counts stay additive: a species row carries its clade's reads, and a row above it carries only the reads that stopped at that rank. This matters more here than it looks — on the CSI-Microbes plate every false positive the negative-control filter removed sat at *genus*, and kraken-biom's own default of `--max O` would have started the table below the ranks where that signal lives.

BIOM 1.0 (JSON) is the default because phyloseq's `import_biom` reads it without `rhdf5`. For BIOM 2.1 (HDF5), either pass `ext.args2 = '--format hdf5'` to `KRAKENBIOM` or run `biom convert -i table.biom -o table.hdf5.biom --to-hdf5` afterwards.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: the aggregate report.
  - `multiqc_data/`: parsed numbers behind every plot.
  - `multiqc_plots/`: static images of the plots.

</details>

[MultiQC](http://multiqc.info) collects FastQC (raw and trimmed, shown as separate sections), fastp, HISAT2, samtools, Qualimap and Kraken2 into a single page. The General Statistics table carries the host alignment rate per sample, which is usually the fastest way to spot a sample that behaved differently from the rest.

### Functional profiling

<details markdown="1">
<summary>Output files</summary>

- `humann/`
  - `<sample>_genefamilies.tsv.gz`, `<sample>_pathabundance.tsv.gz`: HUMAnN output.
  - `taxonomic_profile/<sample>.mpa.txt`: Bracken's kreport (or Kraken2's under `--skip_bracken`) translated into MetaPhlAn format, which is what HUMAnN used to choose pangenomes. The leading `#mpa_v30_CHOCOPhlAn_201901` line is the database declaration HUMAnN 3.6 refuses to run without; it also states which ChocoPhlAn the profile can be matched against — see [usage](usage.md#functional-profiling-optional).
- `humann/regrouped/`, `humann/normalised/`: gene families regrouped (default: KEGG orthologs) and renormalised (default: CPM).

</details>

Only written with `--run_humann`. These are the files the bundled [exploreMetaTax](../apps/exploreMetaTax/) app reads in its Gene families and Stratified tabs.

### Read accounting

<details markdown="1">
<summary>Output files</summary>

- `metaphlan/` (with `--run_metaphlan`)
  - `<sample>_metaphlan.txt`: per-sample marker-gene profile.
  - `metaphlan_combined.tsv`: all samples in one table.
  - `log/<sample>.metaphlan.log`: the run log, including the marker-mapping rate.

- `sylph/` (with `--run_sylph`)
  - `<sample>.sylph.tsv`: sylph's output for one sample, one row per genome called, with its adjusted ANI, taxonomic and sequence abundance, and effective coverage. `Sample_file` holds the sample id rather than the read file name.
  - `sylph_combined.genomes.tsv`: every sample's rows in one table.
  - `taxprof/<sample>.sylphmpa` (with `--sylph_taxonomy`): the sylph-tax taxonomic profile, MetaPhlAn-style, with ANI and coverage on the genome (`t__`) rows.
  - `sylph_combined.relative_abundance.tsv`, `sylph_combined.sequence_abundance.tsv` (with `--sylph_taxonomy`): clades by samples. Relative abundance is coverage-normalised, as MetaPhlAn reports; sequence abundance is the share of reads assigned, as a Kraken2 percentage is. A sample with nothing detected is a column of zeros, not a missing column.

- `differential_abundance/` (with `--da_metadata`)
  - `<comparison>.<method>.results.tsv`: one row per taxon tested - `lfc`, `pvalue`, `qvalue`, `significant`.
  - `<comparison>.<method>.volcano.png`: log2 fold change against -log10(p), significant taxa labelled.

- `read_accounting/`
  - `reanatax.read_accounting.tsv`: one row per sample, raw reads through to classified reads.
  - `reanatax_polya_mqc.tsv`: poly(A) carry-over diagnostic, with `--run_polya_check`.

</details>

The first table to open when a sample looks odd. `host_carryover_pct` and `nonhost_pct_of_raw` also appear in MultiQC's general statistics; a high carry-over means host depletion is leaking rather than that the sample is unusual.

### AI annotations

<details markdown="1">
<summary>Output files</summary>

- `ai/`
  - `taxonomy.ai_insight.md`: the LLM's narrative of the combined taxonomic profile.
- `multiqc_ai/`
  - `multiqc_report.html`: the MultiQC report with a summary under each section, and with the LLM endpoint scrubbed out.
  - `multiqc_data/section_mappings.txt`: which exported data table each section's summary was built from — the first thing to check if a summary looks like it is describing the wrong plot.
  - `multiqc_data/section_prompts.txt`: which sections were sent, and the model used.
- `qualimap_ai/<sample>/`
  - `qualimapReport.html`: the Qualimap report with an AI summary box under its Summary heading.

</details>

These directories only exist when `--llm_endpoint` was given. `qualimap_ai/` is an annotated copy alongside the pristine `qualimap/`; `multiqc_ai/` **replaces** `multiqc/`, because the unannotated report still contains the run's command line and therefore the endpoint. See [usage.md](usage.md) for the full behaviour.

Every box carries the model name, the timestamp, the generation rate and the table it was built from, under a "MUST always verify against the data" label. They are a reading aid for a report you should still read.

A sample whose summary reads `AI summary timed out. Please try again` was attempted and did not get a response within `--llm_timeout`; one with no box at all was either not selected in `--ai_insights` or had no data table to send.

## Exploring the results interactively

The taxonomic tables are meant to be explored, not just read. **exploreMetaTax**
is a Shiny app for exactly that — upload the `kraken2/` and `bracken/` tables
from your `--outdir`, add a metadata sheet, and get filtering, alpha and beta
diversity, rarefaction, PCA, a LifemapR taxonomy tree and publication-ready
exports.

Two ways to run it:

- **Hosted (recommended):** <https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax>
- **Self-contained:** [`apps/exploreMetaTax/`](../apps/exploreMetaTax/) ships a
  [shinylive](https://posit-dev.github.io/r-shinylive/) build that runs entirely
  in your browser — no server, and your data never leaves your machine. Serve it
  with any static web server; see that folder's README.

> [!TIP]
> The hosted app is likely to be **more up to date and considerably faster**.
> Prefer it, and use the self-contained build when the hosted instance is
> unreachable or your data must not leave your machine. If the self-contained
> version feels slow or misbehaves, try the hosted one before reporting a
> problem: R compiled to WebAssembly is several times slower than native R, and
> large datasets feel it.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`, `pipeline_dag_*.html`: Nextflow run reports.
  - `reanatax_software_mqc_versions.yml`: version of every tool that ran.
  - `params_*.json`: the parameters the run used.

</details>

The trace file records CPU, memory and walltime actually used per task — the basis for tuning resource requests on your own data.

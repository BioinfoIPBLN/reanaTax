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

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: the aggregate report.
  - `multiqc_data/`: parsed numbers behind every plot.
  - `multiqc_plots/`: static images of the plots.

</details>

[MultiQC](http://multiqc.info) collects FastQC (raw and trimmed, shown as separate sections), fastp, HISAT2, samtools, Qualimap and Kraken2 into a single page. The General Statistics table carries the host alignment rate per sample, which is usually the fastest way to spot a sample that behaved differently from the rest.

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

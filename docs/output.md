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

[MultiQC](http://multiqc.info) collects FastQC (raw and trimmed, shown as separate sections), fastp, HISAT2, samtools and Kraken2 into a single page. The General Statistics table carries the host alignment rate per sample, which is usually the fastest way to spot a sample that behaved differently from the rest.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report_*.html`, `execution_timeline_*.html`, `execution_trace_*.txt`, `pipeline_dag_*.html`: Nextflow run reports.
  - `reanatax_software_mqc_versions.yml`: version of every tool that ran.
  - `params_*.json`: the parameters the run used.

</details>

The trace file records CPU, memory and walltime actually used per task — the basis for tuning resource requests on your own data.

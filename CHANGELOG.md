# BioinfoIPBLN/reanatax: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

Initial release of BioinfoIPBLN/reanatax, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Three input routes: a samplesheet (`--input`), a folder of FASTQ files with automatic mate detection (`--input_dir`), or ENA/SRA accessions downloaded with fastq-dl (`--input_accessions`).
- Accessions are resolved to their run lists before downloading, so an umbrella accession such as a BioProject fans out into parallel per-run downloads. Runs are then grouped back into samples by `--group_runs_by`.
- Host genome handling: use a prebuilt HISAT2 index, a local FASTA, or let the pipeline fetch the assembly from NCBI by accession or taxonomy ID.
- Host depletion with `hisat2 --very-sensitive`, keeping both the aligned reads (sorted, indexed BAM plus samtools statistics) and the unaligned reads (FASTQ).
- Taxonomic classification of the non-host fraction with Kraken2, abundance re-estimation with Bracken, per-sample Krona charts and combined cross-sample tables.
- Alignment QC of the host BAM: samtools stats/flagstat/idxstats plus Qualimap BamQC, both fed into MultiQC.
- Optional AI annotations against any OpenAI-compatible endpoint, ported from [reanalyzerGSE](https://github.com/BioinfoIPBLN/reanalyzerGSE): a narrative of the combined taxonomic profile at the top of the MultiQC report, per-section summaries throughout it, and a summary box in each Qualimap report. Entirely opt-in via `--llm_endpoint`, text-only, strictly serialised, and the endpoint/API key are scrubbed from the published reports.
- `local` and `slurm` execution profiles, with per-process resource tuning and capped, retried downloads.

### `Fixed`

### `Dependencies`

### `Deprecated`

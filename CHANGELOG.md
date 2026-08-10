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
- A self-contained [shinylive](https://posit-dev.github.io/r-shinylive/) build of the exploreMetaTax Shiny app under `apps/exploreMetaTax/`, for exploring the Kraken2/Bracken tables in a browser with no server. The hosted app at <https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax> remains the faster and more current option.
- Depletion against up to two host references in turn (e.g. GRCh38 then T2T-CHM13): `--fasta`, `--hisat2_index`, `--host_accession` and `--host_taxid` accept comma-separated lists and can be mixed.
- Per-sample read accounting (`read_accounting/`) with a MultiQC read-fate bargraph, plus a host carry-over metric in the general statistics that audits host depletion using the host genome already present in the Kraken2 database.
- `--bracken_read_length auto`: the trimmed read length is measured per sample from fastp and snapped to a k-mer distribution the Bracken database actually ships, instead of being assumed.
- A sparse-taxon filter on the combined tables (`--min_rel_abundance`, `--min_samples`), with the dropped taxa written out rather than silently discarded.
- Optional HUMAnN 3 functional profiling of the non-host fraction (`--run_humann`), with the taxonomic profile derived from the Kraken2 report via KrakenTools' `kreport2mpa.py` rather than a second classifier, plus regrouping to KEGG orthologs and renormalisation.
- Optional host gene quantification with featureCounts (`--quantify_host`), so host expression and microbial composition come from the same library.
- Optional poly(A) carry-over diagnostic (`--run_polya_check`) comparing internal poly-A/T runs against a length- and composition-matched permutation null.
- `local` and `slurm` execution profiles, with per-process resource tuning and capped, retried downloads.

### `Fixed`

- Corrected the documented effect of `--very-sensitive`: `-k 50` costs alignment time, not disk. The HISAT2 module pipes through `samtools view -F 256`, so secondary alignments never reach the BAM — which is also why the BAM can be counted directly with featureCounts.
- `hasClassifiedReads` no longer uses a `while` loop, which `nextflow lint` rejects.

### `Dependencies`

### `Deprecated`

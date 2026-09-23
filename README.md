# BioinfoIPBLN/reanatax

[![GitHub Actions CI Status](https://github.com/BioinfoIPBLN/reanatax/actions/workflows/nf-test.yml/badge.svg)](https://github.com/BioinfoIPBLN/reanatax/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/BioinfoIPBLN/reanatax/actions/workflows/linting.yml/badge.svg)](https://github.com/BioinfoIPBLN/reanatax/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/BioinfoIPBLN/reanatax)

## Introduction

**BioinfoIPBLN/reanatax** takes raw sequencing data - either fetched straight from the public archives or already on disk - strips out the host, and taxonomically profiles what is left. It is built for the reanalysis scenario: you have a BioProject accession from a paper and want to know which non-host organisms are in that data.

Given a BioProject/SRA/ENA accession, a folder of FASTQ files, or a samplesheet, the pipeline downloads and QCs the reads, trims them, aligns them against a host genome that it can fetch from NCBI for you, keeps **both** halves of that split (the host BAM and the non-host FASTQ), and classifies the non-host fraction with Kraken2 and Bracken. Everything lands in a single MultiQC report plus per-sample and combined abundance tables.

1. Fetch reads from ENA/SRA ([`fastq-dl`](https://github.com/rpetit3/fastq-dl)) - each accession is first resolved to its runs so they download in parallel
2. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
3. Adapter and quality trimming ([`fastp`](https://github.com/OpenGene/fastp)), followed by a second `FastQC`
4. Host genome retrieval ([`ncbi-genome-download`](https://github.com/kblin/ncbi-genome-download)) and indexing ([`HISAT2`](https://daehwankimlab.github.io/hisat2/))
5. Host depletion (`HISAT2 --very-sensitive`), saving the aligned reads as sorted, indexed BAM ([`SAMtools`](http://www.htslib.org/)) and the unaligned reads as FASTQ
6. Alignment QC of the host BAM ([`SAMtools`](http://www.htslib.org/) stats/flagstat/idxstats and [`Qualimap`](http://qualimap.conesalab.org/) BamQC)
7. Taxonomic classification of the non-host fraction ([`Kraken2`](https://ccb.jhu.edu/software/kraken2/))
8. Abundance re-estimation ([`Bracken`](https://github.com/jenniferlu717/Bracken)) and interactive charts ([`Krona`](https://github.com/marbl/Krona))
9. Aggregate report ([`MultiQC`](http://multiqc.info/))
10. Optional AI summaries written into those reports by any OpenAI-compatible LLM endpoint (off unless `--llm_endpoint` is given)

Results can be explored interactively with **exploreMetaTax**, either at
<https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax> (recommended — more up to
date and faster) or with the browser-only [shinylive build](apps/exploreMetaTax/)
shipped in this repository, which needs no server and keeps your data local.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

Pick exactly one of the three input routes.

**From public accessions** - a BioProject, study, sample, experiment or run:

```bash
nextflow run BioinfoIPBLN/reanatax \
   -profile local,singularity \
   --input_accessions PRJNA682076 \
   --host_accession GCF_000001405.40 --ncbi_group vertebrate_mammalian \
   --kraken2_db /data/kraken2/Standard \
   --outdir ./results
```

**From a folder of FASTQ files** - mates are paired from their file names:

```bash
nextflow run BioinfoIPBLN/reanatax \
   -profile local,singularity \
   --input_dir /data/my_reads \
   --fasta /data/genomes/host.fa.gz \
   --kraken2_db /data/kraken2/Standard \
   --outdir ./results
```

**From a samplesheet** - use this when you need explicit control over sample names or have several runs per sample:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
```

```bash
nextflow run BioinfoIPBLN/reanatax \
   -profile slurm,singularity --slurm_queue <partition> \
   --input samplesheet.csv \
   --hisat2_index /data/genomes/host_hisat2 \
   --kraken2_db /data/kraken2/Standard \
   --outdir ./results
```

Add `-profile local` to run on the current machine (CPU and RAM ceilings are detected automatically) or `-profile slurm` to submit every task to the scheduler.

See [docs/usage.md](docs/usage.md) for the full parameter reference - in particular how to size the Kraken2 memory request, and why `--hisat2_max_alignments 1` is usually worth adding.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

BioinfoIPBLN/reanatax was originally written by [Jose L. Ruiz](https://github.com/josruirod) (Functional Genomics Center Zurich, University of Zurich / ETH Zurich).

## AI-assisted development

AI coding agents were used to assist with code, comments, and documentation. The design and implementation of the software remain the responsibility of the authors, and all AI-assisted changes were reviewed and tested by the authors.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- If you use BioinfoIPBLN/reanatax for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).

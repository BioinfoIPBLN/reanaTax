# BioinfoIPBLN/reanatax: Usage

## Introduction

reanatax takes raw sequencing data — either downloaded from the public archives or already on disk — removes the host fraction, and taxonomically profiles what is left. It is built for the reanalysis scenario: you have a BioProject accession from a paper and want to know what non-host organisms are in it.

The pipeline runs, in order:

1. **Data retrieval** — [`fastq-dl`](https://github.com/rpetit3/fastq-dl) (only when accessions are given)
2. **Read QC** — [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. **Trimming** — [`fastp`](https://github.com/OpenGene/fastp), followed by a second `FastQC`
4. **Host depletion** — [`HISAT2`](https://daehwankimlab.github.io/hisat2/), keeping both the aligned BAM and the unaligned FASTQ
5. **Classification** — [`Kraken2`](https://ccb.jhu.edu/software/kraken2/) → [`Bracken`](https://github.com/jenniferlu717/Bracken) → [`Krona`](https://github.com/marbl/Krona)
6. **Reporting** — [`MultiQC`](https://multiqc.info/)

## Choosing an input route

Exactly one of the three input options must be given. The pipeline stops with an explicit error if you give zero or more than one.

### 1. Public accessions (`--input_accessions`)

```bash
--input_accessions PRJNA682076
--input_accessions 'SRX9626017,ERX1234253,SRR13191702'
--input_accessions accessions.txt
```

`accessions.txt` holds one accession per line; blank lines and `#` comments are ignored.

Supported accession types are exactly those `fastq-dl` accepts:

| Type       | Prefixes            | Example      |
| ---------- | ------------------- | ------------ |
| BioProject | PRJEB, PRJNA, PRJDB | `PRJNA480016` |
| Study      | ERP, DRP, SRP       | `SRP158268`   |
| BioSample  | SAMD, SAME, SAMN    | `SAMN06479985` |
| Sample     | ERS, DRS, SRS       | `SRS2024210`  |
| Experiment | ERX, DRX, SRX       | `SRX4563689`  |
| Run        | ERR, DRR, SRR       | `SRR7706354`  |

> [!NOTE]
> GEO accessions (`GSE*`/`GSM*`) are **not** supported by fastq-dl. Open the GEO page and use the linked SRA study (`SRP…`) or BioProject (`PRJNA…`) instead. The pipeline detects GEO accessions and tells you this rather than failing later.

Each accession is first resolved to its list of runs, and every run is then downloaded as its own task. A 200-run BioProject therefore downloads with up to `--max_download_forks` transfers in flight instead of one long serial job.

**Grouping runs into samples.** Submissions routinely split one library across several runs. `--group_runs_by` decides what becomes a sample:

- `experiment` (default) — runs of the same experiment are concatenated. This is what most submissions mean by "one sample".
- `sample` — everything sequenced from the same BioSample is concatenated.
- `run` — every run stays its own sample. Use this when a study mixes single- and paired-end runs under one experiment.

**Being a good archive client.** `--max_download_forks` (default `4`) caps concurrent transfers. ENA and SRA throttle aggressive clients, so raise it gradually. On clusters whose compute nodes have no outbound internet, add `--download_local` to run just the download tasks on the submission host.

### 2. A folder of FASTQ files (`--input_dir`)

```bash
--input_dir /data/my_reads
```

Files are found with `--fastq_pattern` (default `**.{fastq,fq}{,.gz}`, which matches both the folder itself and any subdirectory) and mates are paired from the file names. A trailing `_1`/`_2`, `_R1`/`_R2` or `_R1_001`/`_R2_001` (with `.` or `_` as the separator) marks a mate; anything else is treated as single-end. Sample names are whatever precedes that suffix.

If your sample names legitimately end in `_1`, or the folder is single-end only, add `--single_end` to switch mate detection off entirely.

Folders with more than two files per sample are rejected — use a samplesheet for those.

### 3. A samplesheet (`--input`)

```bash
--input samplesheet.csv
```

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
```

| Column    | Description                                                                                                       |
| --------- | ----------------------------------------------------------------------------------------------------------------- |
| `sample`  | Sample name. Repeat it across rows to declare multiple runs of the same sample; they are concatenated.              |
| `fastq_1` | Full path to the FASTQ file for read 1. Must end `.fastq.gz` or `.fq.gz`.                                            |
| `fastq_2` | Full path to the FASTQ file for read 2. Leave empty for single-end data.                                             |

All runs of one sample must have the same endedness.

## The host genome

Give exactly one of these, in decreasing order of precedence:

| Option             | Use when                                                                        |
| ------------------ | ------------------------------------------------------------------------------- |
| `--hisat2_index`   | You already have an index. Fastest — nothing is built. Accepts a directory or a `.tar.gz`. |
| `--fasta`          | You have the genome on disk (plain or gzipped). The index is built once and reused within the run. |
| `--host_accession` | You know the assembly, e.g. `GCF_000001405.40`.                                  |
| `--host_taxid`     | You only know the organism, e.g. `9606`. Resolves to that taxon's RefSeq reference assembly. |

`GCF_` accessions come from RefSeq and `GCA_` from GenBank; the right section is chosen for you. Set `--ncbi_group` to the organism's group (`vertebrate_mammalian`, `bacteria`, `fungi`, `plant`, `viral`, …) — the default `all` works but downloads every group's assembly summary first and is noticeably slower.

`--save_reference` (on by default) publishes the genome and the index under `<outdir>/reference/`, so subsequent runs can skip the build with `--hisat2_index <outdir>/reference/hisat2`.

**Splice-aware indexing.** Pass `--gtf` if the input is RNA-seq and you want splice-aware host capture. Building such an index for a vertebrate genome needs roughly 200 GB of RAM (`--hisat2_build_memory`); below that threshold hisat2-build silently falls back to a genome-only index.

To skip host depletion entirely and classify the trimmed reads directly, use `--skip_host_removal`.

### A note on `--very-sensitive`

The default `--hisat2_args '--very-sensitive'` expands to `--bowtie2-dp 2 -k 50 --score-min L,0,-1`. The `-k 50` part means HISAT2 reports **up to 50 alignments per read**, which makes the host BAM several times larger and the alignment slower — without changing which reads are classified as host.

If the host BAM is only there for QC or for counting host reads, add:

```bash
--hisat2_max_alignments 1
```

This keeps the sensitivity of the search (`--bowtie2-dp 2`, `--score-min L,0,-1`) but writes a single best alignment per read.

The non-host FASTQs are HISAT2's `--un-conc-gz` output, i.e. pairs that did not align **concordantly**. That is deliberately the conservative choice for depletion: a pair where only one mate hit the host still goes forward to classification.

## QC of the host BAM

When `--save_host_bam` is on (the default) the sorted BAM is put through two complementary QC steps, both of which land in the MultiQC report:

- **samtools** `stats`, `flagstat` and `idxstats` — read-level counts: how many reads mapped, how many pairs are proper, and the per-reference breakdown.
- **Qualimap BamQC** — the things counts cannot tell you: coverage depth and how evenly it is spread over the reference, duplication rate, GC of the mapped reads, and the mapping-quality and insert-size distributions.

That second set is what distinguishes *"this library barely contains host"* from *"the host reference is wrong"*. Both produce a low alignment rate; only the latter also shows patchy coverage concentrated in a few repetitive regions at low mapping quality.

Qualimap is the slowest part of the BAM QC, so `--skip_qualimap` turns it off for very large BAMs. Feature-level statistics need an annotation, which is opt-in and deliberately separate from `--gtf`:

```bash
--qualimap_gff /data/genomes/host.gtf
```

`--gtf` only controls whether the HISAT2 index is splice-aware; handing a full vertebrate annotation to BamQC makes it dramatically slower, so it is not inherited.

## The Kraken2 database

`--kraken2_db` is required (unless `--skip_kraken2`). It accepts a database directory or a `.tar.gz` of one. Prebuilt databases are published at <https://benlangmead.github.io/aws-indexes/k2>.

**Memory is the thing to get right.** Kraken2 loads the whole database into RAM. The default request is 72 GB, which suits `Standard-8` or `PlusPF` but not `core_nt` (~700 GB). Either raise the request:

```bash
--kraken2_memory '700.GB'
```

or read the database from disk instead:

```bash
--kraken2_memory_mapping
```

Memory mapping drops the request to 16 GB. It is slower in the worst case, but when the database already sits in the page cache — a fast shared filesystem on a node with plenty of free RAM — it costs almost nothing and is what makes many-sample runs against a large database practical.

**Bracken** re-estimates abundances from the Kraken2 report. It needs a `databaseNmers.kmer_distrib` file in the database matching `--bracken_read_length` (default `100`); check the read length FastQC reports and pick the closest value the database provides. Point `--bracken_db` elsewhere if the distributions live outside the Kraken2 database directory. `--skip_bracken` turns the step off.

`--kraken2_report_minimizer_data` adds distinct-minimizer columns that are useful for filtering false positives, but neither Bracken nor MultiQC can read the resulting report — the pipeline requires `--skip_bracken` alongside it.

## AI annotations (optional)

The pipeline can have a large language model write short summaries into the HTML reports. This is a port of the AI features in [reanalyzerGSE](https://github.com/BioinfoIPBLN/reanalyzerGSE) and shares their design and their `bin/llm_common.py` plumbing.

**Nothing happens unless you ask for it.** With `--llm_endpoint` unset — the default — no AI process runs and no network call is made.

```bash
nextflow run BioinfoIPBLN/reanatax \
   -profile local,singularity \
   --input_dir /data/my_reads --fasta host.fa --kraken2_db /data/kraken2/Standard \
   --llm_endpoint http://your-llm-host:8000/v1/chat/completions \
   --llm_model your-model-name \
   --outdir ./results
```

Any OpenAI-compatible `/v1/chat/completions` endpoint works: a self-hosted vLLM, Ollama or llama.cpp server, or a commercial API. `--llm_api_key` defaults to `dummy`, which is what most local servers expect.

### What gets annotated

`--ai_insights` selects which annotations to generate. It takes `all` (default), `no`, or a comma-separated subset:

| Value      | What it does                                                                                                                          |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `taxonomy` | Reads the combined Kraken2 and Bracken tables and writes a few sentences on the dominant taxa, how consistent samples are, and what looks like contamination. Appears as the **first section of the MultiQC report** and as `ai/taxonomy.ai_insight.md`. |
| `multiqc`  | One summary under every section of the MultiQC report, generated from that section's own exported data table. Works for every MultiQC module, not just the ones this pipeline runs. |
| `qualimap` | A summary box at the top of each per-sample Qualimap report.                                                                            |

`--multiqc_ai_builtin` additionally turns on MultiQC's own AI feature, which writes one summary of the whole report. It is complementary to `--ai_insights multiqc`, which annotates the individual sections.

`--ai_qualimap_sections` adds a box to every individual Qualimap plot as well. That is one extra request per plot per sample, issued one at a time, so for a real cohort it means hours of waiting — hence off by default.

### How it behaves

- **Text only.** Only tables the pipeline already produced are sent — never a figure, never a read, never a BAM. `bin/llm_common.py` rejects multimodal message content outright, so no caller can change that by accident.
- **Strictly sequential.** At most one request is in flight at any moment, because a self-hosted model is usually one server that copes badly with concurrency. This is enforced by the shape of the workflow (`LLM_INSIGHT` → `MULTIQC` → `MULTIQC_AI` → `QUALIMAP_AI` is a chain) plus `maxForks = 1`, with an in-process lock and a `flock` as a backstop. Raising `maxForks` will break it.
- **Never fatal.** An unreachable endpoint, an HTTP error or an empty answer costs you the summary and nothing else; the reports are written either way. Transient failures are retried three times (10 s, 30 s backoff). On a timeout — `--llm_timeout`, 300 s by default — the box says so rather than silently disappearing, so you can tell a failed summary from one that was never requested.
- **Implausible answers are dropped.** A response longer than `LLM_MAX_ANSWER_CHARS` (4000) is treated as model degradation and not shown.
- **Always verify.** Every box is stamped with the model, the date and the table it was built from, under a "MUST always verify against the data" label. Treat the output as a reading aid, not a result.

### Handling of the endpoint and key

The endpoint and API key are treated as deployment secrets; the **model name is not**, and is deliberately kept in every box as provenance.

- They are kept out of the parameter summary that is embedded in the MultiQC report (`validation.summary.hideParams`).
- `MULTIQC_AI` scrubs them from the report and from `multiqc_data/` before publishing — which is why, with `--llm_endpoint` set, the raw MultiQC report is not published at all and `multiqc_ai/` holds the only copy. The pipeline's own methods section quotes the full command line, so the raw report would otherwise contain whatever you typed.
- The parameter dump in `pipeline_info/params_*.json` is scrubbed at the end of the run.

Two places are **not** covered, because nothing in the pipeline can reach them:

- Nextflow's work directory keeps every task's `.command.sh`, as it does for any parameter.
- `pipeline_info/execution_report_*.html` quotes the run command line and each task script, and Nextflow writes it after the pipeline's last hook. Passing the credentials through `-params-file secrets.json` or a private `-c` config keeps them off the command line, which removes most of this.

## Running the pipeline

Typical invocation:

```bash
nextflow run BioinfoIPBLN/reanatax \
    -profile local,singularity \
    --input_accessions PRJNA682076 \
    --host_accession GCF_000001405.40 \
    --ncbi_group vertebrate_mammalian \
    --kraken2_db /data/kraken2/Standard \
    --outdir ./results
```

This creates in your working directory:

```console
work/                # Nextflow scratch, safe to delete once finished
results/             # everything from --outdir
.nextflow.log        # log file from Nextflow
```

Use `-resume` to restart a failed or modified run from where it left off, and `-params-file params.yaml` to keep the parameters under version control instead of on the command line:

```yaml title="params.yaml"
input_accessions: "PRJNA682076"
host_accession: "GCF_000001405.40"
kraken2_db: "/data/kraken2/Standard"
outdir: "./results"
```

### Local execution

```bash
-profile local,singularity
```

CPU and memory ceilings are detected from the machine, and Nextflow then schedules tasks so that the running tasks' combined requests never exceed them — which is what keeps a 72 GB Kraken2 task from being started four times at once. Override the detection with `--resource_limit_cpus`, `--resource_limit_memory` and `--resource_limit_time`.

### SLURM execution

```bash
-profile slurm,singularity --slurm_queue <partition>
```

Every process becomes one `sbatch` job requesting exactly the CPU/memory/time of its resource label, so the scheduler decides what runs concurrently. Useful knobs:

| Option                    | Purpose                                                                 |
| ------------------------- | ----------------------------------------------------------------------- |
| `--slurm_queue`           | Partition to submit to.                                                  |
| `--slurm_options`         | Appended to every `sbatch`, e.g. `'--account=myproject'`.                |
| `--slurm_queue_size`      | Jobs Nextflow keeps queued at once (default `100`).                     |
| `--resource_limit_*`      | Cap per-job requests when the partition's nodes are smaller than the defaults. |
| `--download_local`        | Keep download tasks on the submission host when compute nodes are firewalled. |

Submissions are rate-limited to 20 jobs per minute and failed jobs are retried up to three times, because cluster jobs die for reasons — preemption, node failure, a filesystem hiccup — that have nothing to do with the pipeline.

### Resource tuning

Defaults live in [`conf/base.config`](../conf/base.config), with per-process overrides for the steps whose needs are unusual: downloads (network-bound, 2 CPUs, capped concurrency, retried), `hisat2-build` (memory scales with whether a GTF was given), `HISAT2_ALIGN` and `KRAKEN2_KRAKEN2`.

To change anything else, write a config and pass it with `-c`:

```groovy title="custom.config"
process {
    withName: 'HISAT2_ALIGN' {
        cpus   = 24
        memory = 64.GB
    }
}
```

`-c` never overrides parameters — use `--param` or `-params-file` for those.

## Reproducibility

Pin the pipeline version with `-r`:

```bash
nextflow run BioinfoIPBLN/reanatax -r 1.0.0 ...
```

and prefer `-profile singularity`, `docker`, `apptainer`, `podman` or `conda` over installing the tools yourself. Tool versions are recorded in `<outdir>/pipeline_info/reanatax_software_mqc_versions.yml` and in the MultiQC report.

If you use Singularity and want to reuse images across runs, set:

```bash
export NXF_SINGULARITY_CACHEDIR=/path/to/cache
```

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen.

### `-profile`

Configuration presets, comma-separated. Order matters: later profiles override earlier ones.

Execution: `local`, `slurm`
Containers/environments: `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `conda`, `mamba`, `wave`
Testing: `test`, `test_accession`, `test_full`

> [!TIP]
> We highly recommend using Docker or Singularity containers for full pipeline reproducibility. `conda` is supported as a fallback where containers are not possible.

### `-resume`

Restart from the last successful step. Processes whose inputs and code are unchanged reuse their cached results. You can also supply a specific run name or session ID (`nextflow log` lists them).

### `-c`

Supply an additional config file, for example to change resource requests or add an institutional profile.

## Custom configuration

See the [nf-core configuration docs](https://nf-co.re/docs/usage/getting_started/configuration) for the full picture, including how to point Nextflow at an existing institutional profile with `-profile <institute>`.

## Running in the background

```bash
nextflow run BioinfoIPBLN/reanatax -bg ...
```

Or use `screen`/`tmux`, or submit the Nextflow process itself as a job.

## Nextflow memory requirements

The Nextflow Java process itself can claim excessive memory. Cap it in `~/.bashrc`:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

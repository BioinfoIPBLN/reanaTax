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

## Removing rRNA (optional)

For a total-RNA / RiboZero library, ribosomal RNA is typically 94-99% of the
reads. Kraken2 will happily classify it - rRNA is conserved, so it misassigns
across taxa - and it swamps HUMAnN. `--remove_rrna` runs SortMeRNA to strip it.

```bash
nextflow run . -profile local,singularity \
  --input_accessions PRJNA379542 \
  --remove_rrna --sortmerna_db /data/rrna/smr_v4.3_sensitive_db.fasta \
  --host GCF_000001405.40 \
  --kraken2_db /data/k2_core_nt \
  --outdir results
```

Leave it off for poly(A)-selected RNA or for DNA libraries: there is no rRNA to
remove and the step is pure cost.

The reference FASTAs are not bundled with SortMeRNA. They ship as a separate
`database.tar.gz` download from the SortMeRNA releases page (Silva 138 SSURef
NR99 + RFAM 14.1), giving `smr_*_default_db.fasta`, `*_fast_db.fasta` and
`*_sensitive_db.fasta`. Pass one or several, comma-separated.

### Where it runs, and why

Between trimming and host depletion - **not** on the non-host fraction, which is
where Sola-Leyva et al. (Hum Reprod 2021) put it. Two reasons:

- **Cost.** Depleting rRNA first shrinks HISAT2's input by up to 20x, instead of
  making HISAT2 align a library that is almost entirely rRNA.
- **Accounting.** `classified_pct` is derived from the last HISAT2 pass's
  unaligned count, which is exactly Kraken2's input. A filter placed after
  HISAT2 makes that denominator silently too large; placed before, every
  downstream number stays true by construction.

The trade is that host rRNA is removed before host depletion can count it, so
`host_removed_reads` gets smaller. That is the honest number - an rRNA read is
neither host transcriptome nor microbial signal.

Pairs are kept in step with `--paired_in`, which flags **both** mates as rRNA
when **either** aligns. The surviving output therefore holds only pairs in which
neither mate is rRNA - the same rule as `--require_both_mates_unmapped` in host
depletion, and it means no `repair.sh` re-pairing pass is needed.

`rrna_removed_reads` and `rrna_pct` are added to the read-accounting table, and
the SortMeRNA log reaches MultiQC as a `% rRNA` general-statistics column. The
accounting figure is computed as trimmed reads minus what the first host pass
received, not read off the log: with `--paired_in` the log's per-read E-value
counts overstate what actually left the step.

### Index reuse and versions

The index is built once per run and shared by every sample. Point
`--sortmerna_index` at a prebuilt one to skip that.

The pipeline pins SortMeRNA **4.3.7**, which is the newest build on Bioconda.
Versions 5, 6 and 7 are published on conda-forge only. If you want one of those:

```groovy
process {
    withName: 'SORTMERNA_INDEX|SORTMERNA_READS' {
        container = '<your sortmerna 7 image>'
    }
}
```

The command-line options used here are unchanged across those releases, but
**SortMeRNA 6.0 changed the index format** (CMPH to BBHash), so an index built
by 4.x cannot be read by 6.x/7.x. Rebuild it, or drop `--sortmerna_index`.

## The host genome

`--host` takes the reference in whatever form you have it, and works out which it is from the value:

| You pass                                | Detected as | Notes                                                              |
| --------------------------------------- | ----------- | ------------------------------------------------------------------ |
| `GCF_000001405.40` / `GCA_009914755.4`  | Accession   | Downloaded from NCBI.                                               |
| `9606`                                  | Taxonomy ID | Resolves to that taxon's RefSeq reference assembly.                 |
| a directory, or a `.tar.gz`/`.tgz`      | HISAT2 index | Fastest — nothing is built. The directory must hold `*.ht2` files. |
| any other path or URL                   | Genome FASTA | Plain or gzipped. The index is built once and reused in the run.   |

```bash
--host GCF_000001405.40                    # download it
--host /data/genomes/GRCh38.fa             # build from a local genome
--host results_smoke/reference/hisat2      # reuse an index built earlier
```

`--fasta`, `--hisat2_index`, `--host_accession` and `--host_taxid` still work, but they are grouped by kind rather than kept in the order you wrote them, so a mixed set is always ordered index, FASTA, accession, taxid. `--host` preserves its own order, which is what makes the primary reference yours to choose.

`GCF_` accessions come from RefSeq and `GCA_` from GenBank; the right section is chosen for you. Set `--ncbi_group` to the organism's group (`vertebrate_mammalian`, `bacteria`, `fungi`, `plant`, `viral`, …) — the default `all` works but downloads every group's assembly summary first and is noticeably slower.

`--save_reference` (on by default) publishes the genome and the index under `<outdir>/reference/`, so subsequent runs can skip the build with `--host <outdir>/reference/hisat2`.

**Splice-aware indexing.** Pass `--gtf` if the input is RNA-seq and you want splice-aware host capture. Building such an index for a vertebrate genome needs roughly 200 GB of RAM (`--hisat2_build_memory`); below that threshold hisat2-build silently falls back to a genome-only index.

To skip host depletion entirely and classify the trimmed reads directly, use `--skip_host_removal`.

### A note on `--very-sensitive`

The default `--hisat2_args '--very-sensitive'` expands to `--bowtie2-dp 2 -k 50 --score-min L,0,-1`. The `-k 50` part makes HISAT2 look for **up to 50 alignments per read**, which costs alignment time without changing which reads are classified as host.

It does **not** bloat the BAM: the module pipes HISAT2 through `samtools view -F 256`, so secondary alignments are discarded before anything is written and the BAM holds one primary alignment per aligned read. That is also what makes the BAM directly usable for `--gtf`-driven host quantification.

**`--hisat2_max_alignments` does not work with a preset.** It appends `-k N` *after* `--very-sensitive`, and HISAT2 ignores it — measured, not assumed: a run with `--very-sensitive -k 1` produced identical alignment rates (0.44% / 28.77% on Swab 26, to two decimals) and an identical NH-tag distribution to a `-k 50` run. The preset is applied last and wins. The pipeline warns when you combine them.

To actually cap the alignments, replace the preset with its expansion:

```bash
--hisat2_args '--bowtie2-dp 2 --score-min L,0,-1' --hisat2_max_alignments 1
```

That keeps the sensitivity of the search and only stops HISAT2 hunting for alignments that are thrown away anyway. Untested at the time of writing.

The non-host FASTQs are HISAT2's `--un-conc-gz` output, i.e. pairs that did not align **concordantly**. A pair where only one mate hit the host therefore still goes forward to classification, and is then usually classified as host.

### `--require_both_mates_unmapped`

This switches to a variant module that drops `--no-mixed --no-discordant` and splits the alignment stream on explicit SAM flags instead: `-F 4 -F 8` is the host BAM, `-f 12` (both mates unmapped) is the non-host FASTQ. Pairing is never broken — HISAT2 emits both mates adjacently and `samtools view` preserves that order, so `samtools fastq` stays in sync without a collate or re-pairing pass.

**On by default.** Measured across four vulvar-swab RNA libraries (PRJNA1392516):

| | carry-over before | after | non-host reads lost |
| --- | --- | --- | --- |
| Swab 7 | 0.52% | 0.00% | 2.85% |
| Swab 21 | 0.57% | 0.00% | 3.58% |
| Swab 22 | 2.93% | 0.00% | 3.59% |
| Swab 26 | 13.44% | 0.00% | 2.70% |

The two costs are not equivalent, which is why this is the default rather than an opt-in. The loss is **uniform** — ~3% of reads from every sample, composition unchanged, so it cannot bias comparisons between samples. Carry-over is **sample-specific**, it ranged 0.5–13.5% across four libraries from one study, and you cannot know it is low until after you have measured it. On the matching DNA libraries it reached 62%.

Set `--require_both_mates_unmapped false` to get the `--un-conc-gz` behaviour back. Note it is slower: HISAT2 now looks for discordant and single-mate alignments it previously skipped.

### Multiple host references

`--host` accepts up to **four** references, comma-separated and freely mixed in kind. Reads are depleted against them in the order given — a read has to fail against **every one** of them to be called non-host — and the first is the primary one that `--gtf`, Qualimap and `--quantify_host` describe:

```bash
--host /data/genomes/GRCh38.fa,GCA_009914755.4        # local genome, then download T2T-CHM13
--host results_smoke/reference/hisat2,GCF_009914755.1 # reuse an index, then download T2T-CHM13
--host GCF_000005575.2,GCF_000001405.40,/refs/UniVec.fa  # mosquito, then its blood meal, then vector
```

A fifth reference is refused up front rather than dropped, because a host reference that is silently ignored leaves host sequence in the non-host fraction and nothing downstream can tell that from a real microbial signal. The limit is the number of depletion passes the workflow declares — Nextflow needs one alias per invocation of a subworkflow, so the passes are named at compile time. Raising it means adding a `HOST_DEPLETION_HISAT2` alias and a `PREPARE_HOST_REFERENCE` alias in `workflows/reanatax.nf`, its two naming rules in `conf/modules.config`, and a higher `maxHostPasses()`.

**Two assemblies of one species** is the common case and worth the second pass on its own. GRCh38 is missing most centromeric and satellite sequence and a good deal of structurally variant genome; reads from those regions do not align to it, survive depletion, and then get classified as whatever microbe shares a k-mer with them. T2T-CHM13 contains that sequence, so the second pass removes them. Monteleone et al. (*Microbiome* 2026) use exactly this pairing for the same reason.

**Further references** are for sequence the library genuinely contains but the primary assembly does not:

- a **second organism** — the blood meal in an engorged insect, a xenograft's host, a co-cultured feeder line, the diet in a gut sample;
- a **synthetic reference** — UniVec, PhiX, or the cloning vector a construct came from. These are not hosts in any biological sense, but they are exactly the same operation: align, keep what fails to align. UniVec in particular is small enough that its index build is free next to a genome's.

The **first** reference given is the primary one: it is the one `--gtf` describes, the one Qualimap reports on, and the one `--quantify_host` counts against. Later passes get no GTF — an annotation belongs to one assembly, and handing it to another would build a nonsense splice index — and no Qualimap report, since one reference's leftovers aligned to a different reference answer no question.

Each pass is named apart in the outputs (`<id>.host1`, `<id>.host2`, `<id>.host3`, and `<id>.host` for the last one), so `read_accounting.tsv` reports what each reference removed separately. Whether you needed a given reference is answerable after the fact from that table, and from the **host carry-over** column described below.

### Synthetic sequence (`--univec`)

UniVec is NCBI's set of vector, adapter, linker and primer sequences. Reads
from them are common in real libraries, and they classify to whatever organism
the vector was derived from — which for most cloning vectors is *E. coli* or a
coliphage.

```bash
--univec
```

It is a **separate pass**, not another `--host` entry, and for one reason: a
vector database needs different alignment stringency from a genome.
`--very-sensitive` expands to `--score-min L,0,-1`, which lets a 100 bp read
carry roughly twenty-five mismatches and still align. Against a genome that is
what catches diverged, repetitive and structurally variant host sequence.
Against a 3 kb plasmid backbone it deletes real reads. So this pass does **not**
inherit `--hisat2_args`; it runs end-to-end at HISAT2's own default score
threshold, and `--univec_hisat2_args` replaces the arguments wholesale rather
than adding to them.

| Parameter | Default | Meaning |
|---|---|---|
| `--univec_source` | NCBI `UniVec_Core` | the database; a local FASTA works offline |
| `--univec_hisat2_args` | `--no-spliced-alignment --score-min L,0,-0.2` | replaces `--hisat2_args` for this pass |
| `--univec_save_bam` | false | keep the vector alignments |

**Read this before turning it on.** UniVec's plasmid and phage backbones are
built from real *E. coli* and coliphage sequence — pBR322, the pUC series, M13,
lambda. Depleting against UniVec therefore removes genuine *E. coli* and
coliphage reads along with the vector artefacts. On a gut, stool or clinical
isolate study where those organisms are part of the question, this deletes
signal, not noise. On a tissue or swab study where they are not expected, it is
a cheap win. `--univec_source` accepts the full `UniVec` set as well as
`UniVec_Core`; Core is the reduced-redundancy subset and is the safer default.

It runs **last**, after every `--host` pass, so `unaligned/` holds its
leftovers and the host pass stops publishing there — exactly one process
publishes non-host reads, whatever the configuration. `read_accounting.tsv`
gains a `univec_*` block and its `nonhost_reads` is taken from this pass. It
also works with `--skip_host_removal`, which screens vectors with no host
genome at all.

### Disk (`--compression_level`, `--alignment_output_format`, `--cleanup_intermediates`)

A reanalysis run is large on disk in a way the classification step gets blamed
for and is not responsible for. One CSI-Microbes plate cohort — 180 wells,
110 Gbp — left **1.2 TB** in `work/`: roughly 496 GB of alignments and 413 GB of
FASTQ. Three knobs address that, and they are independent.

| Parameter | Default | What it does |
|---|---|---|
| `--compression_level` | 6 | level for every BGZF/gzip stream the pipeline writes |
| `--alignment_output_format` | `bam` | `cram` for the sorted host alignments |
| `--cleanup_intermediates` | `false` | delete intermediates once nothing reads them |

**`--compression_level`** reaches `samtools sort -l`, `samtools view -l`,
`samtools fastq -c`, `fastp -z` and the `pigz`/`gzip` calls in the local
modules. 6 is the zlib and samtools default; 9 buys roughly 5–10% on BAM for
2–3× the CPU, which is rarely the right trade on a large cohort, and 1 is worth
considering for a run whose intermediates are deleted anyway. fastp refuses 0,
so it is clamped to 1 there. The vendored nf-core `HISAT2_ALIGN` writes its
intermediate at its own default and is the one gap — its output is re-encoded by
the sort, which does honour the setting. The **split** aligner, which is what a
chunked run uses, honours it throughout.

**`--alignment_output_format cram`** typically halves the host alignments. It
comes with three hard requirements, all checked at startup rather than left to
fail per sample at the end of a run:

- **`--host` is required.** CRAM stores differences from a reference, so the
  file is unreadable without one.
- **Qualimap is refused.** `qualimap bamqc` takes `-bam` and cannot open a CRAM.
  Add `--skip_qualimap`; samtools stats, flagstat and idxstats all read CRAM and
  still run.
- **`--quantify_host` is refused.** featureCounts (Subread) reads SAM and BAM
  only.

The reference is indexed once per depletion pass rather than rebuilt inside
every task, which is what htslib would otherwise do — several gigabytes, per
sample. Targeted alignments stay BAM whatever this is set to: they hold one
taxon's reads, so the saving would be negligible, and they are encoded against
the target genome rather than against `--host`.

**`--cleanup_intermediates`** removes two things, each as soon as the pipeline
can prove nothing else reads it:

- **The per-chunk HISAT2 BAMs**, from inside the merge task once `samtools cat`
  has returned. Done there rather than in a cleanup process of its own, because
  "the merge succeeded" is then a property of the script rather than of a
  channel join — `set -e` is in effect, so a failed merge never reaches those
  lines. Only with `--hisat2_chunk_size`, since otherwise there are no chunks.
- **The downloaded FASTQs**, once trimming and raw FastQC have both read them.
  Ordering is enforced by data: the cleanup task takes those tasks' *outputs* as
  inputs, so it cannot be scheduled before them.

**Nothing outside the work directory is ever deleted.** A staged filename is a
symlink, so it is the link target that would have to go, and a target is removed
only when it lies under `workDir`. A FASTQ named in a samplesheet is your file,
wherever you put it; it is skipped, recorded in the cleanup log, and left alone.
Under `--skip_trimming` the downloads are the working read set and are likewise
left alone — the pipeline says so rather than silently doing half the job.

**Off by default, because a deleted intermediate cannot be resumed.** Nextflow
will re-run whatever produced it, and for the downloads that means fetching them
from the archive again. Turn it on for a run you are confident about, not for
one you are still debugging.

### Chunked alignment (`--hisat2_chunk_size`)

HISAT2's memory grows with how much a single process has aligned, and only a fresh process resets it. With `--very-sensitive` on libraries of a few hundred million pairs this is not a small effect — a single alignment has been measured at 68 GB peak RSS, and larger libraries have been reported in the hundreds of GB. `--hisat2_chunk_size` aligns the library in slices, one HISAT2 process each, so the ceiling is set by the chunk rather than by the library:

```bash
--hisat2_chunk_size 5000000                              # 5 M read pairs per task
--hisat2_chunk_size 5000000 --hisat2_chunk_memory 32.GB  # ...and what to request for one
```

**The slices are offsets, not files.** HISAT2's `-s/--skip` and `-u/--upto` count reads *or pairs*, so each chunk task opens the same FASTQ and aligns its own window. Nothing is split, recompressed, or written to scratch — which matters, because physically splitting a 70 Gbp library costs a full decompress/recompress pass and doubles the disk while it runs. The price instead is that each task decompresses past its own offset, which is minutes against hours of `--very-sensitive` alignment.

**The bigger win is usually parallelism.** The chunks are independent tasks, so the longest step in the pipeline stops being one very long task and becomes as many short ones as the executor will schedule. A library that took 16 hours in one pass runs in roughly one, given the cores.

The per-chunk BAM, summary log and non-host FASTQs are joined back into exactly the files an unchunked run produces, under the same names:

- **BAM** — `samtools cat`, which concatenates without decompressing or re-sorting; the sort happens downstream as it always did. Every chunk is given the same `--rg-id`, because `samtools cat` takes the header of the first file for the whole output and a per-chunk read group would leave later records pointing at an `@RG` that is not in the merged header.
- **Summary log** — summed by `bin/merge_hisat2_summary.py` and re-emitted in HISAT2's own format with the percentages recomputed. Every count is additive because the chunks are disjoint slices of one library. MultiQC and `read_accounting.tsv` therefore see one alignment per sample, not fifty.
- **Non-host FASTQs** — concatenated. `cat a.gz b.gz` is itself a valid gzip stream, so there is no decompression pass.

Chunking applies to every depletion pass, and composes with multiple `--host` references.

> **Chunk size is a memory knob, not a speed knob.** Smaller chunks bound memory harder and parallelise further, but each task pays the index load and the decompression to reach its offset, so very small chunks waste both. A few million pairs is a reasonable starting point. `--hisat2_chunk_memory` is separate because the peak a chunk reaches scales with the chunk size, which the pipeline cannot infer.

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

## Host gene counts

`--quantify_host` (with `--gtf`) runs featureCounts on the host BAM and writes a per-sample count matrix to `host_counts/`.

The point is not the counts on their own — it is that they come from the *same library* as the microbial profile. Host expression and microbiome composition measured from one set of molecules share their technical batch effects, so correlating them is defensible in a way that correlating a separate RNA-seq run against a separate 16S run is not. That is the core methodological argument of Monteleone et al. (*Microbiome* 2026), and `--quantify_host` is what puts both halves in your `--outdir`.

No special alignment settings are needed: the HISAT2 module already pipes through `samtools view -F 256`, so each read is counted once.

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

**Bracken** re-estimates abundances from the Kraken2 report. It needs a `databaseNmers.kmer_distrib` file in the database matching `--bracken_read_length` (default `auto`, which measures the trimmed length from fastp and snaps it to the nearest distribution the database ships); pass a number to pin it instead. Point `--bracken_db` elsewhere if the distributions live outside the Kraken2 database directory. `--skip_bracken` turns the step off.

**Stringency.** `--kraken2_min_hit_groups` (Kraken2's `--minimum-hit-groups`) defaults to `2`, which is Kraken2's own default. For host-dominated or low-biomass libraries — RNA-seq mined for microbial signal is both — raising it is the standard tightening, and Monteleone et al. use `3`:

```bash
--kraken2_min_hit_groups 3
```

Each hit group is a run of consecutive k-mers matching the same taxon, so requiring more of them demands that a classification rest on more than one lucky k-mer. It costs sensitivity for rare taxa; that is the trade you are making.

**Sparse taxa.** The combined tables are filtered to taxa exceeding `--min_rel_abundance` (default `0.001`, i.e. 0.1%) in at least `--min_samples` samples, following the same paper. Per-sample reports are left untouched — prevalence is a cross-sample property. This is a sparsity filter, **not** a contamination caller: with no negative controls it cannot tell a reagent contaminant from a rare real organism, it only removes taxa too sparse for any statistic to speak about. What it dropped is written to a `.removed.tsv` beside each filtered table.

`--kraken2_report_minimizer_data` adds distinct-minimizer columns that are useful for filtering false positives, but neither Bracken nor MultiQC can read the resulting report — the pipeline requires `--skip_bracken` alongside it.

### Filtering on evidence, not abundance (`--minimizer_filter`)

`--min_rel_abundance` asks how many reads landed on a taxon. That is precisely
the question a spurious call passes: reads piled on one conserved gene or one
repeat look abundant. What separates them from a real organism is **how much of
the reference those reads cover** — and both classifiers report it.

```bash
# Kraken2 route: the columns need the flag that writes them
--kraken2_report_minimizer_data --minimizer_filter

# KrakenUniq route: its reports carry kmers/dup/cov already
--krakenuniq_db /data/krakenuniq_std --skip_kraken2 --minimizer_filter
```

Three quantities are derived, whichever classifier produced the report:

| | meaning | parameter |
| --- | --- | --- |
| `distinct` | distinct minimizers (Kraken2) or k-mers (KrakenUniq) the reads cover | `--minimizer_min_distinct` (10) |
| `duplication` | observed ÷ distinct — how often the same ones are re-hit | `--minimizer_max_duplication` (100) |
| `coverage` | distinct observed ÷ distinct held in the database for that clade | `--minimizer_min_coverage` (0, off) |

A taxon holding fewer than `--minimizer_min_reads` (50) reads is **never judged**:
too few reads to demonstrate breadth either way, and condemning them would just
re-derive the abundance filter that already ran. The unclassified and root rows
are skipped too — bookkeeping, not organisms.

Coverage is off by default because it is confounded with depth: a genuinely rare
organism also covers little of its reference. Turn it on when your libraries are
deep enough that a real taxon should be well covered.

#### Where the numbers come from

Kraken2's `--report-minimizer-data` inserts two columns — minimizers observed for
the clade, and how many of those were distinct. Its own manual's example makes
the point: 182 reads assigned to an influenza strain, but only **18 distinct
minimizers** behind them. KrakenUniq reports the same idea on k-mers rather than
minimizers (roughly five times the resolution at Kraken2's default k=35, l=31),
and adds `cov` outright.

For the Kraken2 route the coverage denominator is read from `inspect.txt` in the
database directory, whose second column is the distinct minimizers held per
clade; the coverage test is skipped if that file is not there.

#### Output, and how to calibrate

`minimizer_filter/*.minimizer_evidence.tsv` lists every taxon with its reads,
distinct count, duplication, coverage and verdict. Read that **before** trusting
the defaults — they are deliberately conservative, chosen to catch egregious
cases without deleting real low-abundance taxa, not tuned to any particular
library. Sort by `duplication` descending and see where the real organisms stop.

The taxids that fail are handed to the abundance filter rather than removed on
the spot, so one step owns every removal from the combined tables and one
`.removed.tsv` records them all.

#### Does the evidence scale? (`--minimizer_correlation`)

The thresholds above ask whether a taxon's evidence is broad **enough**.
`--minimizer_correlation` adds a second, independent question, taken from
SAHMI's sample-level denoising: does that evidence **scale**? Across samples, a
real organism that is more abundant in one library contributes proportionally
more k-mers and more distinct k-mers there too, so all three of

```
reads ~ k-mers        reads ~ distinct        k-mers ~ distinct
```

rise together. A reagent contaminant sitting at a flat low level in every
library breaks that proportionality — its read count moves and its distinct
count does not. All three Spearman correlations must clear
`--minimizer_correlation_p` (0.05, BH-adjusted), at genus or species rank, for a
taxon seen in at least three samples.

The two tests catch different things and neither subsumes the other, so they are
combined by **union**: a taxon is dropped if either condemns it. A constant
contaminant with plenty of distinct minimizers passes the thresholds and fails
the correlation; a genuine bloom in a single library is judged on its thresholds
alone.

> **At least five samples.** With n samples the smallest attainable two-sided
> Spearman p-value is 2/n!, so at n=4 it is 0.083 and no taxon can clear 0.05
> however clean its evidence. The run is refused rather than silently dropping
> everything it tested.

The evidence table gains `rho_reads_observed`, `rho_reads_distinct`,
`rho_observed_distinct` and `cor_q_max`; taxa the test could not judge — wrong
rank, too few samples, or a count that never varies — read `NA` and keep
whatever the thresholds said.

### Host leakage at k-mer level (`--host_kmer_filter`)

Host depletion is an alignment problem, and alignment is not exhaustive.
Whatever HISAT2 fails to place stays in the non-host fraction, gets classified,
and is reported as a microbe. Neither the abundance nor the minimizer filter can
see this: a carried-over host read is a genuine read with genuine k-mers.

SAHMI attacks it one level down, and `--host_kmer_filter` applies the same test
to bulk data. Kraken2's read-level output records which taxon each *run of
k-mers* in a read was assigned to, so a read carrying host k-mers is
recognisable even when the read as a whole was called a bacterium — which is
what happens to chimeric and repeat-derived reads. A taxon whose reads are
mostly host-tainted is host leakage wearing a species name.

```bash
--kraken2_save_readclassifications --host_kmer_filter
```

| Parameter | Default | Meaning |
|---|---|---|
| `--host_kmer_taxid` | `--host_carryover_taxid` | taxid(s) counted as host, comma-separated |
| `--host_kmer_max_fraction` | 0.5 | drop a taxon once this share of its reads are tainted |
| `--host_kmer_min_reads` | 10 | below this, leave the taxon to the abundance filter |

Two things to know:

- **The host must be in the Kraken2 database**, or no k-mer can ever be assigned
  to it and every count is zero. `core_nt` contains *Anopheles gambiae* (taxid
  7165); PlusPF does not. The filter **fails the run** rather than reporting a
  clean result when it sees no host k-mer at all.
- **`--host_kmer_taxid` takes a list.** The host species alone by default; adding
  its genus or family catches k-mers Kraken2 could only place higher up, at the
  cost of condemning their genuine relatives too.

`host_kmer_filter/*.host_kmer_evidence.tsv` gives, per taxon, its reads, how
many carried host k-mers and the resulting fraction. The MultiQC report gains a
per-library bar of classified reads with and without host k-mers — which is a
direct measurement of carry-over, not the report-level proxy
`--host_carryover_taxid` gives.

### The shuffled-read negative control (`--shuffle_control`)

Every filter above asks whether a taxon's evidence looks strong. This asks a
different question, and it is the only one here with an external null: **how
much evidence would this database have produced from reads that contain no
biological sequence at all?**

A k-mer classifier matches exact substrings, so its false positives are driven
by base composition. An AT-rich read finds AT-rich genomes; a poly-G tail from
two-colour chemistry finds whatever GC-rich organism happens to be in the
database. `--shuffle_control` classifies every library **twice** — once as
sequenced, and once with each read shuffled so that its length, its GC content
and its dinucleotide frequencies survive and none of its 35-mers do. Whatever
the classifier still reports is what it produces from composition alone.

```bash
--shuffle_control --kraken2_use_daemon
```

The result is a ratio, not a p-value: shuffled reads over real reads, per taxon,
pooled across the cohort. Near zero means the taxon needs genuine sequence to be
found. Near one means the database would have reported it from noise, and its
real count is not evidence.

| Parameter | Default | Meaning |
|---|---|---|
| `--shuffle_method` | `dinuc` | `dinuc`/`mono`/`reverse` |
| `--shuffle_seed` | 1 | fixed; a control that moves between runs is not one |
| `--shuffle_reads` | 0 | shuffle only the first N fragments per sample (0 = all) |
| `--shuffle_max_ratio` | 0.1 | flag a taxon at this shuffled/real ratio |
| `--shuffle_min_reads` | 10 | below this many real reads, leave it untested |

The three shuffles differ in how conservative they are.
`dinuc` is an Altschul–Erickson Euler-path shuffle preserving the exact
dinucleotide frequency as well as the base composition, and is the default
because most compositional bias is dinucleotide-level — a shuffle that discards
it understates how many chance matches real data produces. `mono` permutes the
bases and preserves GC content only. `reverse` reads each sequence backwards
(**not** reverse-complemented, which would still map), preserving every
compositional statistic including the trinucleotide spectrum, at the cost of
leaving palindromic and low-complexity stretches intact — which is a reason to
choose it if a homopolymer artefact is what you are chasing.

Three things it deliberately does not do:

- **It does not test significance.** The shuffled library is one draw, and a
  per-taxon Poisson test against it would dress a single observation up as an
  inference.
- **It does not compare per sample.** Chance matches are rare per taxon per
  sample, so a per-sample ratio is mostly 0/0. Counts are pooled over the
  cohort.
- **It does not see host carry-over.** Carry-over reads are real sequence, and
  shuffling removes them, so a carried-over taxon passes this control
  perfectly. That is `--host_kmer_filter`'s job, and the two are not
  substitutes.

`--shuffle_reads` makes the control cheap and raises its detection floor: one
chance read in a subsample stands for several in the full library, so no taxon
below `scale / --shuffle_max_ratio` real reads could be flagged. That floor is
computed, printed, and taxa under it are reported as `untested_below_floor`
rather than `clean` — an untested taxon has not passed anything.

The second Kraken2 pass runs with **exactly** the same arguments as the first.
A control run at a different `--kraken2_confidence` would be measuring a
different classifier.

### Breadth across the genome (`--gene_diversity_filter`)

Every filter to this point counts reads, or counts the k-mers behind reads.
None of them can see **where** those reads landed — and that is the whole
difference between a species that is present and a species whose name has been
attached to one conserved stretch of sequence. A thousand reads on a 16S gene,
on a ribosomal protein, on a transposase shared across half a phylum, is a
thousand reads, and every abundance threshold passes it.

A genuine organism transcribes hundreds of genes. Its reads fall on many
distinct gene families and no single family carries the bulk of them.

This needs no new reference and no new alignment, because HUMAnN already writes
the table: gene-family abundances **stratified by species**. So the filter
requires `--run_humann`.

```bash
--run_humann --gene_diversity_filter
```

| Parameter | Default | Meaning |
|---|---|---|
| `--gene_diversity_min_families` | 10 | distinct families a taxon must reach |
| `--gene_diversity_max_top_fraction` | 0.5 | ...with no one family carrying more than this |
| `--gene_diversity_min_abundance` | 0 | below this pooled abundance, leave it untested |

With `--humann_regroup` set (the default is `uniref90_ko`) it runs a **second**
time on the regrouped table. That is the difference between many gene families
and many *products*: a taxon can reach dozens of UniRef90 families that all
encode the same thing, and only the product view catches it.

Unlike the other evidence filters, this one cannot feed the abundance filter
inside the classification step — HUMAnN runs downstream of it, and a taxid list
cannot travel backwards. Its verdict is applied by a further filter pass whose
output is `gene_diversity/bracken_combined_<level>_genediv.tsv`, named for what
was applied to it.

The idea is PRISM's — that read count without breadth is not evidence, and that
below about ten reads there is nothing to confirm either way. The metric is this
pipeline's own; PRISM's score is not reproduced.

## Single-cell data (`--single_cell`)

scRNA-seq is not a variant of the bulk route, and the reason is worth stating
plainly: **in bulk, host depletion throws the host reads away. In single cell
the host reads are the other half of the experiment** — the cell-by-gene matrix
— and the cell barcode is the only thing that joins them to the microbial reads.
HISAT2 is barcode-blind, so `--single_cell` **replaces** `HOST_DEPLETION_HISAT2`
with STARsolo, which produces both in one pass. SAHMI uses STARsolo for the same
reason.

```bash
--single_cell --host GRCh38 --gtf gencode.v44.gtf --sc_chemistry 10xv3 --sc_whitelist 3M-february-2018.txt --kraken2_db /path/to/k2_core_nt --kraken2_save_readclassifications
```

`fastq_1` must be the **barcode+UMI** read and `fastq_2` the **cDNA** read.
(STARsolo's own `--readFilesIn` takes them in the opposite order; the module
handles that, so you don't.)

### How the two halves stay consistent

The non-host reads come out as an ordinary **single-end FASTQ with the corrected
barcode folded into each read name**. That matters because Kraken2 truncates a
read name at the first whitespace, and `samtools fastq -T` appends its tags after
a tab — left alone, the barcode would be silently dropped on the way into the
classifier.

Because of that, everything downstream runs unchanged: Kraken2, Bracken, Krona,
every evidence filter, the exploreMetaTax bundle, MultiQC. **The pseudobulk
profile and the cell-by-taxon matrix are built from the same classification, so
they cannot disagree**, and any taxon the filters remove is removed from both.

### Filters on the microbial side

| Parameter | Default | From | Meaning |
|---|---|---|---|
| `--sc_min_frac` | 0.5 | SAHMI | share of a read's k-mers inside the assigned taxon's lineage |
| `--sc_max_homopolymer` | 0 (off) | SAHMI `nFilter` | drop reads with a longer single-base run |
| `--sc_host_taxid` | `--host_kmer_taxid` | SAHMI | a read with *any* k-mer of these taxids is dropped |
| `--sc_ranks` | `S,G` | SAHMI | ranks the matrix is tabulated at |
| `--sc_umi_dedup` | `true` | **not** SAHMI | count distinct UMIs, not reads |

Two of these are deliberate departures from SAHMI, and both are corrections:

- **Barcodes come from STARsolo's corrected `CB` tag.** SAHMI's `sckmer.r` reads
  the barcode positionally as `substr(R1, 1, cb_len)`, with no whitelist and no
  error correction — so one sequencing error in a barcode manufactures a new
  "cell". Reads whose barcode STARsolo could not correct carry `CB:NA` and are
  kept in the pseudobulk but excluded from the matrix.
- **UMIs are deduplicated throughout.** SAHMI's `taxa_counts.r` deduplicates too
  — it takes `unique()` over `(barcode, umi, taxid)`, so its matrix is already UMI
  counts. What it does not do is carry that through to `sckmer.r`, whose k-mer
  statistics — and therefore the barcode-level denoising built on them — are
  computed over **undeduplicated reads**. So SAHMI's matrix and its own denoising
  disagree about what a count is. Here the deduplication applies to both. Set
  `--sc_umi_dedup false` to count reads instead.

`--sc_max_homopolymer` is the one filter that has to read the FASTQ — the
sequences are not in Kraken2's output — so it is skipped entirely when set to 0.
SAHMI's own default of 130 is very permissive: on a 150 bp read it only catches
sequences that are ~87% homopolymer.

### Does a taxon's evidence scale across cells? (`--sc_kmer_denoise`)

This is SAHMI's single-cell contribution, and it is the same question
`--minimizer_correlation` asks across samples — at a far finer grain. A study
has a handful of libraries but a library has **thousands of barcodes**, so this
correlation has real power where the sample-level one needs at least five
samples to reach significance at all.

A barcode carrying more reads of a real organism carries proportionally more
**distinct** k-mers too, because those reads sample more of its genome. Ambient
contamination smeared over every droplet, or reads piling on one conserved
locus, saturates: the total climbs and the distinct count does not. Spearman
over barcodes, BH-adjusted, on taxa seen on more than three barcodes.

The distinct k-mers cannot be read out of the Kraken2 report — that is a
per-sample number and the question is per barcode — so they are recovered from
the reads. Kraken2's per-read output run-length encodes which taxon each
consecutive k-mer matched (`taxid:count taxid:count …` in read order), so the
k-mers belonging to a taxon can be pulled straight out of the sequence.

One case the correlation alone cannot express: **a taxon whose distinct count
never moves while its read count does**. Spearman is undefined there (zero
variance), but the meaning is not — that is the contamination signature in its
purest form, and it is failed outright as `saturated_distinct` rather than left
untested.

| Parameter | Default | Meaning |
|---|---|---|
| `--sc_kmer_len` | 35 | Kraken2's k; 35 for the standard databases |
| `--sc_min_barcodes` | 4 | SAHMI tests taxa on more than three barcodes |
| `--sc_correlation_p` | 0.05 | BH-adjusted |
| `--sc_max_barcodes_per_taxon` | 1000 | SAHMI's `nsample`; bounds memory |

> **The verdicts are per sample and are NOT fed into the cohort abundance
> filter.** Barcodes only mean anything inside the library that produced them,
> and a taxon that fails in one sample may be perfectly real in another. The
> evidence table and drop list are published for you to apply; the bulk filters
> (`--minimizer_filter`, `--host_kmer_filter`, `--decontam`) remain the ones
> that act on the combined tables.

### Which cell types carry which microbes (`--sc_cell_metadata`)

The question the branch exists to answer. The statistics here are
[CSI-Microbes'](https://www.science.org/doi/10.1126/sciadv.adj7402) rather than
SAHMI's, for three reasons that paper argues explicitly:

- **Fisher exact, not chi-square.** The chi-square approximation needs expected
  frequencies of at least 5, and that fails for most sample × cell-type × taxon
  combinations. (Lloréns-Rico et al. use a pooled chi-square; on sparse data
  Fisher is the correct exact alternative.)
- **Presence, not abundance.** The matrix is >90% zeros, so the useful question
  is what *fraction of cells* carry a taxon. This is also why ALDEx2/ANCOM-BC2 —
  which model compositional abundance across samples — are the wrong tool here
  and are not used on the cell matrix.
- **Per sample, then combined — never pooled.** Per-sample Fisher p-values are
  combined by Stouffer's Z weighted by expected infected cells.

```bash
--sc_cell_metadata cells.tsv     # barcode, sample, cell_type [, host_umis]
```

> **The `sample` column is mandatory and the script refuses to run without it.**
> Two samples where cell type 2 is clearly enriched *within each* can pool to
> make cell type 1 look enriched, purely from differing cell-type composition —
> Simpson's paradox. reanaTax exists to reanalyse public data, where samples
> have wildly different compositions, so this is a live risk, not a hypothetical
> one. The implementation is validated against the paper's own counterexample.

Effect size is `log2(observed / expected)` where expected = (cells of that type
÷ cells in sample) × infected cells in sample. Across the cohort **observed and
expected are each summed and the ratio taken last** — never an average of
ratios, which would weight a two-cell sample like a two-thousand-cell one.

Cell types are the one input the pipeline cannot produce — STARsolo gives a
count matrix, not clusters — so bring them from Seurat or scanpy. Without
`--sc_cell_metadata` the step is skipped.

**Co-occurrence** (per-sample hypergeometric → Stouffer → BH) is emitted
alongside, and with it a **doublet check**: cells carrying two taxa may simply be
undetected doublets, which would explain the co-occurrence entirely. The check
compares host UMI counts between multi- and singly-infected cells (needs the
optional `host_umis` column). Treat a significant result there as invalidating
any polymicrobial-cell claim.

### Cell association vs the ambient pool (`--sc_ambient`)

**This does not decide what is contamination, and it is an annotation rather
than a filter for that reason.**

An empty droplet contains at least four things a droplet cannot tell apart:
reagent and kit contaminants; genuinely extracellular microbes that were in the
tissue suspension; microbes released from cells lysed during dissociation; and
ambient nucleic-acid soup. A luminal, mucosal or biofilm organism — which for a
gut, skin, oral or vulvar sample is very often the organism the study is *about*
— lands in that pool by construction, not by artefact. Treating everything in
the empty droplets as contamination would delete it.

What the comparison *can* answer is narrower: **is this taxon spatially
associated with cells?** Three outcomes, and only one of them says anything
about contamination:

| `prevalence_ratio` | verdict | reading |
|---|---|---|
| ≫ 1 | `cell_associated` | intracellular or tightly adherent — the strongest evidence droplet data can offer |
| ≈ 1 | `ambient_undecided` | reagent contaminant **or** genuine extracellular organism. Not resolvable here |
| ≪ 1 | `cell_depleted` | reagent, index hopping, or something that does not survive the cells it came with |

The middle row does not resolve and cannot be made to. The only thing that
separates a kit contaminant from a real extracellular organism is an external
measurement of the kit — which is exactly what [`--decontam`](#reagent-contaminants---decontam)'s
blanks provide and what nothing else in this pipeline does. Read the two
together: a taxon that is `ambient_undecided` here **and** commoner in the
blanks than in the samples is a contaminant; `ambient_undecided` and absent from
the blanks is a real organism that lives outside cells.

```bash
--sc_ambient
```

Two ratios are computed, because they can disagree and the disagreement is
informative. `prevalence_ratio` is the share of cells carrying the taxon over
the share of empty droplets carrying it — robust, and blind to how much is
there. `rate_ratio` is the taxon's UMIs per host UMI in cells over the same in
empty droplets — sensitive to abundance, and the one that moves when a taxon is
everywhere but concentrated in cells. Significance is Fisher's exact on the 2×2
presence table, one-sided in whichever direction the ratio points, BH-adjusted
across taxa.

| Parameter | Default | Meaning |
|---|---|---|
| `--sc_ambient_min_empty_umis` | 1 | host UMIs an empty droplet needs to be used |
| `--sc_ambient_max_empty_umis` | 0 | ...and at most this many (0 = no bound) |
| `--sc_ambient_min_droplets` | 5 | droplets a taxon must appear in to be tested |
| `--sc_ambient_enriched_ratio` | 2.0 | ratio for `cell_associated` |
| `--sc_ambient_depleted_ratio` | 0.5 | ratio for `cell_depleted` |
| `--sc_ambient_drop` | false | write a drop list of ambient + depleted taxa |

`--sc_ambient_min_empty_umis` matters more than it looks: a raw 10x matrix has
hundreds of thousands of barcodes with zero or one UMI, and counting those as
empty droplets inflates the denominator until every taxon looks cell-associated.
The upper bound is there for the opposite problem — barcodes just below the knee
carry enough real content to behave like cells.

**`--sc_ambient_drop` is off by default and should usually stay off.** It is
only the right call when the question is specifically about *intracellular*
microbes. If the study is about the community in the tissue, the ambient taxa
are part of the answer, not noise to be removed.

Needs both halves of STARsolo's `Solo.out` tree — `filtered/` to name the cells,
`raw/` to name the empty droplets — so it is refused with `--sc_cell_filter
None`, which makes STARsolo write only one matrix. The raw matrix is streamed
for its column sums only and never held in memory.

### The host half: infected against bystander cells (`--sc_host_de`)

The cell-by-taxon matrix is the **instrument**, not the result. What SAHMI's
paper actually reports is what a cell's own transcriptome does when it is
carrying something, and that comparison is what everything upstream exists to
make trustworthy.

```bash
--sc_host_de --sc_cell_metadata cells.tsv
```

Wilcoxon rank-sum on log-normalised counts — Seurat's `FindMarkers` defaults,
which is what SAHMI calls — between cells carrying a taxon and **bystander cells
of the same type in the same library**. Both halves come out of the same
STARsolo pass keyed on the same corrected barcode, so the two views cannot
disagree about which cell is which.

**Within a cell type, always.** Infection is not distributed at random over cell
types — measuring that is the entire point of the enrichment step above — so a
pooled infected-vs-uninfected test recovers the difference *between* the cell
types that happen to be infected and reports it as a response *to* infection.
Nothing in the output would distinguish the two. `--sc_cell_metadata` is
therefore required; `--sc_host_de_force_pooled` accepts a pooled comparison and
stamps every row `ALL_POOLED` so it cannot be read as a within-type result.

A bystander is a cell of the same type carrying **none of that taxon** — not a
cell carrying no microbe at all. A cell positive for a different taxon is a
valid bystander here, and excluding it would silently restrict the comparison to
the cleanest cells in the library.

| Parameter | Default | Meaning |
|---|---|---|
| `--sc_host_de_min_cells` | 10 | infected *and* bystander cells a group needs |
| `--sc_host_de_min_pct` | 0.1 | a gene must be detected in this share of one side |
| `--sc_host_de_logfc` | 0.25 | minimum \|log2 fold change\| to test |
| `--sc_host_de_top_taxa` | 20 | most-infecting taxa to test (0 = all) |

`--sc_min_umis` decides what "infected" means, applied after whatever denoising
ran upstream. p-values are adjusted **within** each taxon × cell type: each is a
separate question asked of a separate set of cells, and pooling them would let a
group with thousands of tested genes set the threshold for a group with a
hundred.

### Plate-based data (`--sc_plate_based`)

CSI-Microbes runs on two kinds of single-cell data and treats them differently,
because they are different experiments wearing the same name. In droplet data
one library holds thousands of cells and the barcode separates them. In
plate-based data — Smart-seq2 and its relatives — **every well is sequenced as
its own library**, there is no barcode, and the cell *is* the sample.

That layout needs none of the barcode machinery and no STARsolo. Each cell runs
the ordinary bulk route, and the cell-by-taxon matrix is those per-cell reports
stacked; everything downstream of the matrix is shared with the droplet route.
So `--sc_plate_based` and `--single_cell` are mutually exclusive, and giving
both stops the run.

```bash
--sc_plate_based --sc_cell_metadata cells.tsv
```

The samplesheet is one **cell** per row. `--sc_cell_metadata` is required and
must carry the cell id in the first column, a `cell_type` column, and a
`sample` (or `patient`/`donor`/`plate`) column. That last one is not optional
bookkeeping: the enrichment test stratifies by sample and combines across
strata precisely to avoid Simpson's paradox, and without it every cell would be
its own stratum and the design would silently collapse to the pooled test it
exists to avoid.

`--sc_plate_min_reads` (default 10) decides when a cell carries a taxon.
Deliberately **not** `--sc_min_umis`: Smart-seq2 has no UMIs, and one amplified
cDNA molecule can contribute hundreds of reads and a single UMI, so a read
threshold and a UMI threshold are not interchangeable and the read one has to be
stricter.

### Reference

STARsolo needs its own index, so `--single_cell` builds one with
`STAR_GENOMEGENERATE` (or takes `--star_index`). Everything before the index
build — resolving `--host`, downloading from NCBI, gunzipping — is shared with
the bulk route.

`--gtf` is **required**: STARsolo assigns reads to genes at alignment time, so an
index built without an annotation cannot produce a cell-by-gene matrix and no
later flag can rescue it.

Only the first `--host` reference is used; the multi-pass depletion the bulk
route offers is not available here, since a later pass would re-align reads that
no longer have a barcode read beside them.

### The host half: gene expression and host-microbe correlation

`--quantify_host` counts host reads against `--gtf`. Until now those counts were published per sample and read by nothing; they are now merged into a gene x sample matrix and tested between the **same two groups** the microbial side is tested on, so one `--da_metadata` file and one design describe both halves of the library.

```bash
--quantify_host --gtf host.gtf \
--da_metadata meta.tsv --da_grouping tissue \
--da_sample_group salivary_gland --da_ref_group midgut
```

That runs automatically once `--quantify_host` and the `--da_*` contrast are both present. `--host_de_method` selects the engine — `deseq2`, `edger`, or both, exactly like `--da_method` — and `null` skips it.

> **Not ALDEx2/ANCOM-BC2, deliberately.** Those model *compositional* data: a microbial profile carries no absolute scale, only proportions, and recovering something interpretable from that is their whole purpose. Gene counts do have a scale, up to a library-size factor that DESeq2's median-of-ratios and edgeR's TMM estimate directly. Running the microbial methods on gene counts works mechanically, which is the dangerous part — it answers a different question, badly.

| Parameter | Default | Meaning |
|---|---|---|
| `--host_de_method` | `deseq2,edger` | engine(s); `null` to skip |
| `--host_de_min_count` | 10 | a gene needs this many reads... |
| `--host_de_min_samples` | 2 | ...in at least this many samples |

The contrast, covariates, `--da_p_threshold` and `--da_lfc_threshold` are shared with the microbial side. Filtering happens **before** the test: filtering results afterwards would leave the adjusted p-values computed against a gene set that was never tested. Two samples per group is a hard floor — below that the dispersion cannot be estimated and any p-value is an artefact of the model's fallback.

#### Host expression against microbial abundance (`--host_microbe_correlation`)

The third layer Monteleone et al. describe, and the analysis `--quantify_host` exists for: both sides come from the *same* library, so an association between them is not confounded by sample handling the way two separate assays would be. Spearman on log-CPM host expression against relative microbial abundance, BH-adjusted.

```bash
--host_microbe_correlation --host_microbe_top_genes 2000 --host_microbe_top_taxa 50
```

> **This is opt-in because it is easy to run and impossible to interpret on a small cohort.** With n samples the smallest attainable two-sided Spearman p-value is 2/n!, whatever the data. Testing g genes against t taxa is g x t tests, and under BH the best of them must clear alpha/(g x t). On five samples the floor is 2/120 = 0.0167, so **at most three pairs could ever be significant** — against the millions a gene-by-taxon grid contains. Run blind, the analysis returns an empty table that reads like a negative result and is nothing of the kind.
>
> The step computes that arithmetic up front and **refuses**, quoting it, rather than producing the empty table. On a real four-sample cohort a gene with a perfect rho of -1.0 came back at q = 0.449. Cut the grid with `--host_microbe_top_genes`/`--host_microbe_top_taxa`, add samples, or pass the step `--force` through `ext.args` if you want the table regardless.

Output lands in `host_expression/`: the count matrix, gene lengths, one results table per engine, and the correlation table with the feasibility arithmetic in its header.

### One taxon, followed to its own differential expression (`--target_taxid`)

The rest of the pipeline answers *who is there, and how much*. This answers *and what were they doing*, for one organism chosen after the classification has run.

Kraken2's per-read output records which read it put on which taxon, so those reads can be pulled back out with KrakenTools, aligned to that organism's **own** genome, counted against its annotation, and tested between the same two groups as everything else:

```bash
--kraken2_save_readclassifications \
--target_taxid 5858 \
--target_reference GCF_000956335.1 \
--target_gtf refs/pmalariae.gtf \
--da_metadata meta.tsv --da_grouping tissue \
--da_sample_group salivary_gland --da_ref_group midgut
```

| Parameter | Default | Meaning |
|---|---|---|
| `--target_taxid` | — | taxon to extract; comma-separated for several |
| `--target_reference` | — | its genome: accession, taxid, FASTA or HISAT2 index |
| `--target_gtf` | — | its annotation; without it the branch stops at the BAM |
| `--target_include_children` | `true` | also take reads assigned below the taxon |
| `--target_include_parents` | `false` | ...and above it |
| `--target_max_reads` | 100000000 | KrakenTools' own cap, made explicit |
| `--target_de_method` | `deseq2,edger` | engine(s) for the gene-level test |

Almost none of this is new machinery. Extraction is the only new step: the reference preparation, alignment, sorting and gene counting are the **same subworkflow the host pass uses**, and the merge-and-test tail is the same one the host gene counts go through, both aliased rather than copied so there is one implementation of each to keep correct.

> **The alignment rate is the check, and you should look at it.** Reads reach this branch because a k-mer classifier assigned them, and a classifier's false positives arrive looking exactly like its true positives. Reads that genuinely came from the target organism align to its genome; reads Kraken2 mis-assigned do not. `targeted/log/*.target.hisat2.summary.log` is that number, and it also appears in MultiQC alongside the host alignment as a separate `<sample>.target` row.

> **`--target_include_parents` is off for a reason.** A read Kraken2 could only place at the genus is *consistent with* the target species but is not evidence for it. Turning this on quantifies such reads as though it were, and the pipeline warns when you do.

Output lands in `targeted/`: the extracted FASTQs, the sorted BAM and its samtools stats, per-sample gene counts, the merged matrix, and one DE table per engine.

### Pseudoalignment quantification (`--host_transcripts`, `--target_transcripts`)

Where a transcriptome is available, kallisto replaces align-then-count:

```bash
--host_transcripts refs/agam_transcripts.fa --host_tx2gene refs/agam_tx2gene.tsv
--target_transcripts refs/pmalariae_transcripts.fa --target_tx2gene refs/pm_tx2gene.tsv
```

`--host_transcripts` pseudoaligns the **trimmed** reads, so it is independent of host depletion and of `--quantify_host` — you get host expression whether or not the BAM was kept. `--target_transcripts` replaces the entire align-and-count route for the targeted taxon: no genome, no index, no GTF.

Either way the output is the same gene × sample matrix, feeding the same DESeq2/edgeR path, so nothing downstream changes.

| Parameter | Default | Meaning |
|---|---|---|
| `--host_transcripts` | — | host transcriptome FASTA |
| `--host_tx2gene` | — | two-column TSV to sum transcripts to genes |
| `--target_transcripts` | — | transcriptome for `--target_taxid` |
| `--target_tx2gene` | — | as above, for the target |
| `--kallisto_fragment_length` | 200 | single-end only |
| `--kallisto_fragment_sd` | 30 | single-end only |

> **Not offered for host depletion, deliberately.** A transcriptome index can only remove reads that came from an annotated transcript. Intronic, intergenic, repeat and satellite reads have nowhere to pseudoalign, would survive depletion, and are exactly the reads that go on to be classified as spurious microbes — the failure mode pairing GRCh38 with T2T-CHM13 exists to prevent. It is mechanically possible (`--pseudobam` emits unmapped records); it is a worse depletion, so the pipeline does not offer it as an equal alternative.

Three things worth knowing about the numbers:

- **`est_counts` are estimated and fractional.** A read pseudoaligning to several transcripts of one gene is apportioned between them. DESeq2 is handed rounded counts, which is an approximation; the fully correct route is tximport, which passes effective lengths to the model as offsets. Supplying `--host_tx2gene`/`--target_tx2gene` removes most of the difference, because the apportioning that rounding damages is mostly *within* a gene.
- **Single-end input needs `-l` and `-s`.** kallisto measures the fragment-length distribution from paired reads and cannot infer it from one end, so single-end abundances inherit whatever error your estimate carries.
- **The targeted route loses its sanity check.** On the alignment route, the HISAT2 alignment rate is what distinguishes real assignments from Kraken2's false positives. Pseudoalignment gives `p_pseudoaligned` in `run_info.json` instead, which is the nearest equivalent but is against a transcriptome rather than a genome.

### Carryover and index hopping (`--negative_controls`)

Every filter above, `--decontam` included, reaches **one verdict per taxon**:
the organism is real, or it is not. There is a failure mode that shape of answer
cannot describe, and on a multiplexed plate it is the dominant one.

Index hopping, well-to-well carryover and ambient template put **genuine reads
of a genuine organism into the wrong library**. Nothing about those reads is
wrong — they are the same reads, with the same k-mers, that the neighbouring
well produced correctly. And the taxon is then *signal* in one library and
*carryover* in the next, so no single verdict about the taxon can be right for
both.

This was measured, not assumed. On the CSI-Microbes plate (Robinson et al.
2024) — 124 wells deliberately infected with *Fusobacterium nucleatum*, 22
uninfected, 12 empty — the organism turns up in uninfected wells at 1–78 reads,
and:

- **Evidence filters cannot see it.** The distinct-minimizers-per-read ratio
  that Kraken2's manual and the exploreMetaTax app both recommend has **AUC
  0.415** on those calls. Below 0.5: it ranks the false positives *above* the
  true positives, because the ratio is inversely tied to read count and every
  false positive sits at a handful of reads. Duplication reaches 0.656, raw
  distinct minimizers 0.854 — statistically indistinguishable from raw reads
  (0.858), which is to say it carries no information depth does not.
- **`--decontam`'s prevalence method cannot settle it either**, and no threshold
  makes it. *Fusobacterium* is in 6 of the 12 blanks **and** is the organism the
  experiment is about. A per-taxon call has to either delete it (recall 54.8% →
  0) or spare it (specificity unchanged). Both are wrong.

What works is a per-library test against an external reference level. A control
library measures how much of each taxon arrives *without a sample* — reagent
contamination, ambient template, and, on a plate, the hopping rate. A taxon is
kept in a library only where it exceeds that level by a stated margin.

```bash
--negative_controls 'metadata.tsv:condition:Empty Well' --control_ratio 20 --control_min_reads 30
```

| Parameter | Default | Meaning |
|---|---|---|
| `--negative_controls` | — | which libraries are the controls; four forms below |
| `--control_ratio` | 10 | multiple of the control level a taxon must reach |
| `--control_statistic` | `mean` | `mean`/`median`/`max` of the controls |
| `--control_min_reads` | 2 | reads needed before the ratio is consulted |
| `--control_floor_reads` | 1 | the control level never falls below this many reads |
| `--prevalence_filter` | 0 | drop taxa in this fraction of libraries; 0 = off |
| `--prevalence_min_reads` | 1 | reads for `--prevalence_filter` to count a library |

`--negative_controls` takes any of four forms, so a cohort that already declares
its blanks for `--decontam` does not declare them twice:

| Form | Meaning |
|---|---|
| `SRX1,SRX2` | the sample IDs |
| `controls.txt` | a file of one ID per line |
| `meta.tsv:column` | rows whose column is `true`/`yes`/`1`/`control`/`blank`/`negative` |
| `meta.tsv:column:value` | rows whose column equals `value` exactly |

**Levels are counts per million classified reads, never raw reads.** Depth
normalisation is on its own the single largest improvement available on that
plate — reads-per-million reaches 73.5% specificity at full recall where raw
reads reaches 61.8% — and without it the threshold would mean a different thing
in every library.

Measured on the plate, against `benchmark/csi_microbes`'s own scorecard — 124
infected wells, 22 uninfected, 12 empty as the controls:

| Setting | Specificity | Detected | Wells with any taxon | Scorecard |
|---|---|---|---|---|
| off | 30/34 (88.2%) | 62/124 | 180/180 | **1 FAIL**, 3 PASS |
| `--control_ratio 10` | 32/34 (94.1%) | 62/124 | 180/180 | 1 FAIL, 3 PASS |
| **`--control_ratio 20 --control_min_reads 30`** | **33/34 (97.1%)** | 58/124 | 165/180 | **0 FAIL, 4 PASS** |
| `--control_ratio 60` | 34/34 (100%) | 58/124 | 161/180 | 1 FAIL, 4 PASS |

Two things that table is really saying.

**The ratio is not the only lever, and on its own it is the wrong one.**
Pushing `--control_ratio` alone does reach 100% specificity, but it gets there by
applying the same multiple to *every* taxon, and past about 50× it empties whole
libraries of everything — the run then fails a different check. The four false
positives are all small in absolute terms (11–72 reads), so
`--control_min_reads` removes them at a fraction of the collateral damage. Reach
for the read floor first and the ratio second.

**The recall cost is real and it is not a threshold artefact.** Detection falls
62/124 → 58/124. The worst false positive sits at 4,209 reads per million and
only 4 of the 62 detected wells sit below it, so those four are the price. They
are wells whose entire evidence is a handful of reads, indistinguishable by
construction from the carryover in the well beside them. The filter does not
resolve that ambiguity; it declines to call it.

**The verdicts are computed on the combined Kraken2 report, not on Bracken**,
and on a single-cell run that distinction decides whether the filter works at
all. Bracken's table is species-only; the cell-by-taxon matrix carries whatever
rank Kraken2 assigned each read. On this plate every false positive sat at
*genus* Fusobacterium, absent from the Bracken table — so Bracken-derived
verdicts reached none of them and specificity did not move. The whole-taxon drop
list still comes from Bracken, whose flat table counts each read once.

With `--sc_apply_drop_list` (on by default for single-cell runs) the per-library
verdicts reach the cell matrix too, so the cohort profile and the single-cell
profile cannot disagree about a cell one of them zeroed.

**`--control_floor_reads` is what makes the ratio mean anything.** Absence from
the blanks is a *bound*, not a measurement: it says a taxon is under one read in
each of them, not that it is at zero. Without a floor the test collapses to "any
reads at all" for every taxon the controls happened to miss — which on a
species-level table is most of them. Left at 0, the table above stays flat at
86.4% specificity however high `--control_ratio` goes, because the species that
leaked into the uninfected wells appeared in one blank out of twelve.

Note that **only the failing cells are zeroed**; a taxon that fails in *every*
library leaves by the usual route, as a taxid the abundance filter removes. The
control libraries themselves are never filtered — you cannot judge a control
against itself — and stay in the table as columns.

#### When there are no controls (`--prevalence_filter`)

A much cruder question, for the cohorts that have no blanks at all: is this
taxon in **every** library? A reagent contaminant introduced at extraction is; a
biological signal usually is not.

On the same plate, `--prevalence_filter 1.0` removes **six taxa holding 62.6% of
all microbial reads**, and does not touch *Fusobacterium*:

| Reads | Taxon |
|---|---|
| 1,976,110 | synthetic construct |
| 439,008 | *Escherichia coli* |
| 258,164 | *Malassezia restricta* |
| 56,113 | *Microbacterium foliorum* |
| 17,191 | *Dietzia psychralcaliphila* |
| 16,182 | *Microbacterium maritypicum* |

Every one is on Salter et al.'s (BMC Biol 2014) list of reagent contaminants, or
is a cloning-vector artefact. That is a large effect from a crude rule.

**It is only sound where a universally present taxon cannot be real.** It is
wrong for a mono-culture, wrong for a dominant gut commensal, wrong for anything
with an expected core microbiome. It is off by default, it warns below 0.5, and
it names what it removed in `control_filter/`.

**It is a background remover, not a false-positive remover, and the distinction
is worth being clear about.** On this plate it strips 62.6% of the reads and
leaves the specificity failure exactly where it was — 30/34, unchanged at every
threshold from 1.0 down to 0.9. That is not a defect: *Fusobacterium* is in 101
of 180 libraries, so no prevalence rule can touch it without also deleting it
from the wells where it is real. Reach for it to clear the floor, and for
`--negative_controls` to fix a specificity problem.

**Both are complementary to `--decontam`, not replacements.** decontam is given
a measurement of the *kit* and asks whether a taxon is reagent; this is given
control *libraries* and asks whether there is more of a taxon here than turns up
on its own. Where both are available, run both.

### Reagent contaminants (`--decontam`)

Every filter above judges a taxon by its own evidence — how abundant it is, how
much of its reference its reads cover, whether its counts scale, whether they
are host sequence. **None of them can tell a reagent contaminant from a
genuinely rare organism**, because nothing in the data distinguishes the two.
Only a blank can.

[decontam](https://benjjneb.github.io/decontam/) settles it. Its **prevalence**
method compares how often a taxon appears in true samples against how often it
appears in negative controls; commoner in the blanks means contaminant. Its
**frequency** method needs no blanks but does need a post-PCR DNA concentration
per sample, and tests whether a taxon's relative abundance falls as total DNA
rises — the signature of a constant amount of contaminating template diluted by
varying amounts of real sample.

```bash
--decontam --decontam_metadata metadata.tsv --decontam_neg_column sample_type --decontam_neg_value blank
```

| Parameter | Default | Meaning |
|---|---|---|
| `--decontam_metadata` | `--da_metadata` | TSV; first column is the sample id |
| `--decontam_neg_column` | — | column marking negative controls |
| `--decontam_neg_value` | `true` | value marking one (case-insensitive) |
| `--decontam_conc_column` | — | DNA concentration, for `frequency` |
| `--decontam_method` | `prevalence` | `prevalence`/`frequency`/`combined`/`either`/`minimum` |
| `--decontam_threshold` | 0.1 | decontam's own default |
| `--decontam_batch_column` | — | column grouping samples into batches |
| `--decontam_batch_combine` | `minimum` | `minimum`/`product`/`fisher` |

**`--decontam_batch_column` is worth reaching for whenever the cohort spans more
than one run.** Contamination is a property of the *kit*, not of the study. Pool
two sequencing runs, two extraction days or two centres, and a taxon that was
heavy in one batch and absent from the other looks exactly like a taxon that
varies between samples — which is what a real organism looks like. Naming the
column that separates the runs scores each batch against its own blanks and then
combines the verdicts.

The combination rule matters as much as the split. `minimum` (decontam's own
default) takes the smallest p across batches, so a taxon condemned in any one
batch is condemned overall — the right stance when a contaminant was only
present in one run. `product` and `fisher` require agreement across batches and
are correspondingly conservative.

Every batch must contain at least one control and two true samples, or the run
stops and names the batches that do not. That is deliberate: decontam would
otherwise return `NA` for a starved batch, and a taxon absent from every *other*
batch would then never be scored at all.

The evidence table gains two columns the package does not produce —
`n_batches_flagged` and `batches_flagged` — because with `minimum` a taxon
flagged in one batch of six and a taxon flagged in all six carry the same
verdict and are not the same claim. The first is a contaminant of one run; the
second is a contaminant of the kit.

**This is opt-in and stays opt-in.** Most public datasets carry neither blanks
nor concentrations, and a reanalysis pipeline that required them could not run
on the archives it exists to mine. What it does instead is **fail loudly** when
asked for a method whose inputs are absent — a wrong column name, or a
`--decontam_neg_value` that matches no sample, stops the run with a message
naming the columns it did find, rather than reporting that nothing was
contaminated.

Two practical notes:

- **The blanks must have been run through the pipeline** like any other sample,
  so they appear as columns in the combined table. They stay as columns
  afterwards: nothing here drops samples. Exclude them downstream, or with
  `--da_samples_to_drop` for differential abundance.
- **Scored on the combined Bracken table**, not the Kraken2 one. A combined
  kreport is a hierarchy — a clade's count already contains its children's — so
  a prevalence test on it would count the same reads at every rank. Hence
  `--decontam` cannot be combined with `--skip_bracken`.

## Which taxa travel together (`--run_sparcc`)

A co-occurrence network over the cohort, inferred by SparCC (FastSpar's
implementation).

**Why not a correlation matrix.** Taxonomic profiles are compositional: the
counts sum to a library size that has nothing to do with the biology, so when
one taxon rises every other one falls whether or not anything about them is
related. A Pearson or Spearman matrix over such data manufactures negative
correlations everywhere and a spurious positive block among the rare taxa.
SparCC works on log-ratio variances, which are invariant to the total, and
infers the correlations of the underlying basis instead — the same reason
ANCOM-BC2 and ALDEx2 are this pipeline's differential-abundance methods rather
than a t-test on proportions.

```bash
--run_sparcc
```

| Parameter | Default | Meaning |
|---|---|---|
| `--sparcc_min_prevalence` | 0.5 | a taxon must be seen in this share of samples... |
| `--sparcc_min_reads` | 10 | ...with at least this many reads there |
| `--sparcc_top_taxa` | 200 | most abundant survivors kept (0 = all) |
| `--sparcc_iterations` | 50 | SparCC iterations for the point estimate |
| `--sparcc_exclusion_iterations` | 10 | strongly-correlated-pair exclusion rounds |
| `--sparcc_permutations` | 1000 | bootstrap replicates behind the p-values |
| `--sparcc_min_correlation` | 0.3 | \|rho\| an edge needs before it is reported |
| `--sparcc_p_threshold` | 0.05 | on the **adjusted** p |
| `--sparcc_force` | false | run even when the cohort cannot support it |

**The run is refused rather than emptied.** SparCC estimates each pair from the
variance of a log ratio across samples; with a handful of samples those
variances are noise, and bootstrap p-values will be small for some pair purely
because a resample of five samples has very few distinct outcomes. Below ten
samples the run stops. It also stops when the bootstrap p-value floor
(`1 / --sparcc_permutations`) sits above what Benjamini–Hochberg needs for the
number of pairs being tested — with 200 taxa that is 19,900 pairs, so 1,000
replicates cannot call *any* edge significant, and the message says how many
replicates it would take. An empty network from such a run would mean "this
cohort cannot be asked", not "there are no associations".

Taxa are filtered on prevalence before the tool sees them, because a log ratio
is undefined at zero: SparCC adds a pseudocount, and a taxon that is zero in
most samples has its correlations decided by that pseudocount rather than by its
data.

The p-values are adjusted **across every pair at once**, not per taxon and not
per sign. An edge is one test, and splitting the family by sign would let the
positive edges be judged against a smaller family than the negative ones for no
reason but their sign. An unadjusted SparCC network is mostly noise arranged in
a suggestive shape, which is how co-occurrence networks acquired their
reputation.

## Diversity and what explains it (`--run_diversity`)

ALDEx2 and ANCOM-BC2 answer "which taxa differ between two groups". This answers
the question that usually comes first when reanalysing someone else's data:
**how much of the variation does each recorded variable actually explain, and
which variables are redundant with each other?** A dataset where 30% of the
community variation tracks sequencing batch and 2% tracks the disease of
interest is telling you something no per-taxon test will.

The design follows [Lloréns-Rico et al. (Nat Commun 2021)](https://www.nature.com/articles/s41467-021-26500-8),
whose analysis code is public:

- Counts are **CLR-transformed and compared by Euclidean distance** — together
  the Aitchison distance, used rather than Bray-Curtis because sequencing yields
  compositions, where only ratios between taxa carry information.
- Each variable is tested on its own by **distance-based RDA** (`vegan::capscale`),
  giving an adjusted R² and a permutation p-value, BH-corrected across variables.
- A **stepwise model** (`vegan::ordiR2step`) then adds variables in order of the
  variation they explain that the ones already in the model do not. *The gap
  between a variable's univariate R² and its contribution here is the
  confounding* — two variables can each explain 20% and jointly explain 21%.
- **PERMANOVA** (`vegan::adonis2`) on the same distance, since that is what most
  readers expect quoted.

Alpha diversity (observed richness, Shannon, Simpson, inverse Simpson, Pielou's
evenness) comes out alongside, on raw counts.

```bash
--run_diversity --diversity_metadata metadata.tsv
```

**Metadata is optional.** Without it you still get alpha diversity and an
unconstrained Aitchison PCoA; the variance partitioning simply has nothing to
partition by and is skipped rather than faked. It falls back to `--da_metadata`
so one file can serve both.

Unlike differential abundance, this runs on the **filtered** table. The reason
is the mirror image: ALDEx2 and ANCOM-BC2 filter internally, so pre-thinning
would filter twice; diversity indices filter not at all, and observed richness
on an unfiltered Bracken table is very largely a count of database artefacts.

Two things to watch:

- **Small cohorts cannot resolve anything.** With four samples split 2/2 there
  are three distinct permutations, so the smallest attainable p-value is 0.33.
  A negative adjusted R² is not a bug — it means the variable explains less than
  chance would.
- **The zero replacement matters.** CLR needs positive values, so each zero is
  replaced by half the smallest non-zero count in its own sample. A principled
  treatment (`zCompositions`) models the zeros instead of imputing them, and the
  difference grows with sparsity.

### PRISM (`--run_prism`)

[PRISM](https://github.com/sjdlabgroup/PRISM) is not a classifier and does not
replace one. It takes candidate taxa and tries to **confirm** each of them with
evidence no k-mer count can supply: full-length BLAST alignment against `nt`,
the GenBank annotation those alignments land on, host mapping with STAR and
minimap2, and the k-mer composition of the taxon's reads across every rank.
Forty features per taxon go into an XGBoost model and come out as one number —
the probability that the taxon is genuinely present rather than a contaminant.

```bash
--run_prism --prism_path /refs/PRISM \
  --prism_blast_db /refs/blast/nt \
  --prism_star_genome_dir /refs/star_GRCh38 \
  --prism_minimap2_index /refs/GRCh38.mmi
```

**Single-cell libraries need `--prism_barcode_only`.** On 10x chemistry read 1
is barcode plus UMI with no biological sequence, and without the flag PRISM
classifies and BLASTs it as if it were cDNA; the run is refused rather than
allowed to do that. Setting it **drops the barcode read before PRISM starts**,
so Kraken2, STAR, minimap2 and BLAST all see the cDNA read alone. That is
stricter than PRISM's own `--barcode_only`, which reaches only its BLAST step
while its Kraken2 pass still runs `--paired` over both mates. For PRISM's
literal behaviour instead, leave the flag false and pass
`--prism_args '--barcode_only TRUE'` — the module never sets that flag itself,
so there is nothing to collide with.

**It takes the trimmed reads, not the non-host ones.** This is the only tool
here that does. PRISM performs its own host removal, and several of its forty
features are statistics *of* that removal — how much of the library was host,
how much unclassified, how a candidate's k-mers distribute across ranks. Hand
it a fraction that has already had the host taken out and those columns empty,
which changes what the model is scoring without changing anything that would
look wrong in the output. The pipeline wires it to `ch_trimmed_reads` for you.

**Preparing `--prism_path` is the work.** A clone of the repository is not
enough. The directory must also hold:

- `genbank/genbank_indices.RDS` and `genbank/rds/*.RDS` — the processed GenBank
  annotation, distributed separately by the authors rather than in the repo.
- `sorted_accession_map.txt` — the accession→taxid map, built once with
  `blastdbcmd` against your BLAST database.

plus a BLAST `nt` database, a host STAR index and a host minimap2 index. The
run stops at validation if any of the parameters naming these is missing,
rather than failing several hours in.

**There is no public container.** PRISM is an R script that shells out to
Kraken2, minimap2, STAR, BLAST and SeqKit, and no biocontainer carries that
set. Under `-profile conda` the module's `environment.yml` builds it. Under
docker or singularity you must supply `--prism_container` pointing at an image
you built — the run refuses to start otherwise, rather than trying to pull
something that does not exist.

If you build that image yourself, two things about the dependency list are
worth knowing, because both were found by reading the code rather than the
README:

- **`xgboost` is required and is not listed.** `functions.R` calls
  `predict()` on the `xgb.Booster` it `readRDS()`es without ever attaching the
  package, so `predict.xgb.Booster` is unregistered and the scoring step dies
  on an install that follows the README exactly.
- **Pin xgboost to 1.7.** `prismxg.RDS` was written with `saveRDS` on an
  `xgb.Booster`, which xgboost warns against because the model's `handle` is an
  external pointer, and its fields are the 1.x layout. The 2.x/3.x R interface
  changed how those are restored, so the newest package is the riskiest choice.
  `modules/local/prism/run/environment.yml` pins `r-xgboost=1.7.6` for this
  reason.

| Parameter | Default | Meaning |
|---|---|---|
| `--prism_kraken_db` | `--kraken2_db` | PRISM runs its own Kraken2 pass |
| `--prism_min_qcovs` | 80 | PRISM's minimum BLAST query coverage |
| `--prism_min_read_per` | 10000 | PRISM's minimum reads per X to analyse |
| `--prism_min_uniq_frac` | 5 | PRISM's minimum unique k-mer ratio |
| `--prism_max_sample` | 100 | PRISM's maximum reads to sample |
| `--prism_score_threshold` | 0.5 | median score below which a taxon is a contaminant |
| `--prism_min_reads` | 10 | below this, PRISM has nothing to confirm either way |
| `--prism_filter` | false | apply the verdict to the combined Bracken table |

Scores are aggregated across the cohort as the **median over the samples that
saw the taxon**. Not the mean, which one confidently-scored library drags
either way, and not "flagged in any sample", which on a large cohort condemns
everything.

`--prism_filter` is off by default on purpose. PRISM scores the whole trimmed
library; Bracken profiles the non-host fraction. They agree on taxids, but they
were asked about slightly different sets of reads, so look at
`prism/reanatax.prism_evidence.tsv` before letting it remove anything.

#### How this relates to `--gene_diversity_filter`

They share an argument and nothing else. `gene_div` and `prod_div` are two of
PRISM's forty features, and `--gene_diversity_filter` reproduces those two —
plus `fprod`, `fugene` and `fuprod` — from HUMAnN's species-stratified table,
under PRISM's own names and with PRISM's own definitions (raw Shannon entropy
from `vegan::diversity()`, and shares of the across-taxa total rather than raw
counts).

What it cannot reproduce is the **score**. That needs all forty features, and
the other thirty-five are BLAST multi-mapping ratios, k-mer taxonomy at seven
ranks, and Kraken lineage/misclassification statistics — none of which exist
without PRISM's own BLAST and Kraken passes. The `verdict` column in
`gene_diversity/` is therefore this pipeline's threshold rule, not a PRISM
call. It is also far cheaper: no BLAST against `nt`, no reference bundle. Use
it when PRISM's setup is out of reach, and PRISM itself when it is not.

### GATK PathSeq (`--pathseq_microbe_bwa_image`)

A third classification route, and the only one here that is **not k-mer based**.
Kraken2 and KrakenUniq match exact substrings against a database of genomes;
MetaPhlAn matches clade-specific markers. PathSeq *aligns* the surviving reads
to its microbe reference with BWA and scores taxa from those alignments,
distributing an ambiguously-mapping read's weight across the taxa it hits
instead of pushing it up to their common ancestor.

That difference is the reason to run it. A read Kraken2 can only place at a
genus — because its k-mers are shared — contributes nothing at species level;
PathSeq gives it fractional weight at every species it aligns to and reports
`unambiguous` separately, so the two kinds of evidence stay distinguishable. It
is also alignment-based, so it is not fooled by the composition-driven chance
matches `--shuffle_control` exists to measure. **Where PathSeq and Kraken2 agree
on a taxon, two quite different failure modes have both been avoided.**

```bash
--pathseq_microbe_bwa_image microbe.fa.img \
--pathseq_microbe_dict microbe.dict \
--pathseq_taxonomy_db microbe.db
```

Additive, like every other classifier here; `--skip_kraken2` uses it alone. It
runs on the reads that survived host depletion, and its own host filter — if
`--pathseq_host_bwa_image` and `--pathseq_host_kmers` are both given — then acts
as a second, alignment-plus-k-mer pass over what HISAT2 left behind. Give both
or neither.

| Parameter | Default | Meaning |
|---|---|---|
| `--pathseq_host_bwa_image` | — | host BWA image for PathSeq's own filter |
| `--pathseq_host_kmers` | — | the k-mer file that goes with it |
| `--pathseq_save_bam` | false | the per-read BAM is the size of the library |
| `--pathseq_rank` | `species` | PathSeq's `type` kept in the combined tables |
| `--pathseq_args` | — | passed to `PathSeqPipelineSpark` |

**The resource bundle is the price of entry.** This pipeline does not build the
index images: they are large, slow to build and specific to a reference choice
you have to make anyway. Broad distributes prebuilt ones.

The combined output is kept in PathSeq's **own units** rather than converted
into a Kraken report. The two tools disagree about what "a read assigned to a
taxon" means, and a conversion would have to pick one convention and hide the
other; keeping both means a disagreement between the routes stays visible
instead of being averaged away by the format. Four matrices come out — `score`,
`score_normalized`, `reads` and `unambiguous` — and `unambiguous` is the one to
set beside a Kraken2 species count.

### KrakenUniq (`--krakenuniq_db`)

Kraken2 tells you how MANY reads landed on a taxon. KrakenUniq also tells you how
much of that taxon they actually cover: its report carries per-taxon distinct
k-mer counts (`kmers`, `dup`) estimated by HyperLogLog, and `cov`, the fraction
of the reference those k-mers span. A taxon called from reads piled on one
conserved stretch has few distinct k-mers for its read count and near-zero
coverage - which no abundance filter can see.

```bash
nextflow run . -profile local,singularity \
  --input_accessions PRJNA1392516 \
  --host GCF_000001405.40 \
  --kraken2_db /data/k2_core_nt \
  --krakenuniq_db /data/krakenuniq_std \
  --outdir results
```

**It needs its own database.** KrakenUniq cannot read a Kraken2 index, which is
why this is a separate parameter rather than a switch: build one with
`krakenuniq-build`, or fetch a prebuilt one. Nothing about the route runs unless
`--krakenuniq_db` is set.

**It is additive.** Left as above, both classifiers run over the same non-host
reads, which is the only way to ask whether a taxon Kraken2 called is backed by
distinct k-mers or by reads on a single conserved region. Add `--skip_kraken2` to
use KrakenUniq on its own; read accounting then takes its numbers from the
KrakenUniq report instead.

**The database is loaded once, not once per sample.** All libraries of one
layout are classified in a single task after one `krakenuniq --preload`. Set
`--krakenuniq_args '--preload-size 100G'` to classify against a database larger
than the node's memory; that switches the module out of whole-database
preloading.

**Bracken does not run on this route.** Bracken's k-mer distributions are built
alongside a *Kraken2* database, so pairing them with a KrakenUniq report would
estimate one database's abundances from another's statistics. You get the
per-sample reports, the combined table, the sparse-taxon filter and Krona.

#### Kraken2 can approximate this without a second database

`--kraken2_report_minimizer_data` adds two columns to the Kraken2 report: the
minimizers observed for a taxon and, of those, how many were distinct. That is
the same discriminator measured on minimizers rather than k-mers - roughly a
one-in-five subsample at Kraken2's default k=35, l=31. The database-side
denominator for a coverage figure is already in the `inspect.txt` shipped with
the database, whose second column is the distinct minimizers it holds per clade.

So reach for KrakenUniq when you need confident calls at low abundance, where
the minimizer count gets too thin to discriminate; the minimizer columns cover
the rest at no extra cost.

### Keeping the database in memory (`--kraken2_use_daemon`)

The legacy `kraken2` wrapper reads the entire hash into the process on **every**
invocation. Against a 340 GB `core_nt` on NFS that is roughly half an hour of
pure I/O per task, repeated for every sample. `k2 classify --use-daemon` loads
the index once into a background process that outlives the task; every later
sample skips the load entirely.

```bash
nextflow run . -profile local,singularity \
  --input_accessions SRX8592588,SRX8592589 \
  --host GCF_943734735.2 \
  --kraken2_db /srv/GT/databases/kraken2/k2_core_nt_20251015 \
  --kraken2_use_daemon \
  --kraken2_confidence 0.1 \
  --outdir results
```

Four consequences worth understanding before switching it on.

**Samples classify one at a time.** The daemon is addressed through a single
control FIFO with no locking, so concurrent clients would interleave their
handshakes. This costs nothing to prevent - the index is resident, so each
sample starts classifying immediately, and the alternative was paying the load
repeatedly anyway.

Note that `maxForks 1` does not achieve it on its own. Nextflow counts maxForks
per process, and the module is included twice - once for the reads as
sequenced, once aliased for `--shuffle_control` - so the two ran concurrently,
two clients on one FIFO, which wedged the run and left two private copies of
the index resident (2 x 314 GB for a 344 GB database). The serialisation is
therefore an `flock` the process instances share, which the kernel also
releases if a task is killed.

**The database must be an unpacked directory, given as a stable path.** The
daemon keys its resident index on the `--db` argument, so the path is passed
through unstaged and absolute. A staged symlink would resolve differently in
every task and the daemon would load a second copy alongside the first. A
tarballed `--kraken2_db` is rejected up front for the same reason.

**The daemon outlives the pipeline.** Nothing stops it when the run ends, and it
goes on holding the index - 340 GB, in the core_nt case. On a shared node, stop
it when you are done:

```bash
k2 clean --stop-daemon
```

It is also addressed through hardcoded paths in `/tmp` (`classify.pid`,
`classify_stdin`, `classify_stdout`), which means it is shared with anything
else running k2 on that node, not scoped to your run.

**A dead daemon does not report itself.** Nothing in the handshake has a
deadline. If the daemon dies - it can be killed by memory pressure from
something else on a shared node - a client that opens the control FIFO blocks
in `open(2)` and never returns: no error, no CPU, an empty `.k2.log`, and with
`maxForks 1` the rest of the cohort queued behind it. `k2` decides whether a
daemon is already running from `/tmp/classify.pid` and only asks whether that
PID exists, not whether it belongs to a k2 process, so a recycled PID makes the
check succeed forever.

Each classify call is therefore bounded by `--kraken2_daemon_timeout` (default
`3h`, enough for a cold load of a ~340 GB index plus a large library). On
timeout the task fails and is retried, and the retry stops the daemon and
clears `/tmp` so a fresh one starts - which is why the reset happens on retry
and not up front: the daemon is shared with anything else on the node, and a
healthy one serving another run must not be torn down on suspicion. Raise the
timeout for droplet-scale libraries.

If a run does stall, the check is one line:

```bash
cat /tmp/classify.pid            # then: ps -p <pid> -o comm=,cmd=
```

**The daemon is stopped for you.** It runs outside a PID namespace so that it
survives the task that started it - that is what lets the whole cohort share one
resident index - so nothing would otherwise reap it, and a finished run was
measured still holding 316 GB ninety minutes later. The pipeline stops it in
`workflow.onComplete` and `workflow.onError`, and the benchmark launch scripts
also trap `EXIT INT TERM` for the case where Nextflow itself is killed and no
handler runs. Both only ever signal a PID that is currently a `classify` process
**and** owned by the owner of the pid file, so a stale pid file naming a recycled
PID cannot cause the wrong process to be killed.

`--kraken2_use_daemon` and `--kraken2_memory_mapping` are contradictory - one
exists to hold the index in RAM, the other to avoid doing so - and the pipeline
refuses the combination.

### The alternative: tmpfs plus memory-mapping

If you cannot use the daemon, putting the database in `/dev/shm` and mapping it
gets most of the benefit and lets samples run concurrently:

```bash
mkdir -p /dev/shm/k2_core_nt
cp /path/to/k2_core_nt_20251015/*.k2d /dev/shm/k2_core_nt/
cp /path/to/k2_core_nt_20251015/database*mers.kmer_distrib /dev/shm/k2_core_nt/

--kraken2_db /dev/shm/k2_core_nt --kraken2_memory_mapping --kraken2_memory 32.GB
```

The `.kmer_distrib` files come along because `--bracken_db` defaults to
`--kraken2_db` and the `--bracken_read_length auto` logic globs that directory.
Note that memory-mapping is only a good idea over tmpfs or fast local disk - over
NFS every page fault becomes a network round-trip and it is far slower than
simply loading the database.

## Functional profiling (optional)

Taxonomy says who is present; `--run_humann` says what they are doing. HUMAnN 3 produces gene-family and pathway abundances from the non-host reads, optionally regrouped onto KEGG orthologs — the third layer in Monteleone et al.'s host-mRNA / species / pathway analysis.

```bash
--run_humann \
--humann_nucleotide_db /data/humann/chocophlan \
--humann_protein_db    /data/humann/uniref \
--humann_utility_db    /data/humann/utility_mapping
```

It is off by default because it needs those large databases and is comfortably the most expensive step here.

Three things to know about how it is wired:

- **The taxonomic profile comes from Bracken, not MetaPhlAn.** HUMAnN picks which pangenomes to align against from a taxonomic profile and only reads MetaPhlAn's format, so a Kraken-style report is translated with KrakenTools' `kreport2mpa.py`. The alternative — running MetaPhlAn as a second classifier — would mean a second database and a second set of abundances quietly disagreeing with the Bracken tables in the same report. `--run_humann` therefore cannot be combined with `--skip_kraken2`. Which report gets translated depends on `--skip_bracken`: normally it is Bracken's `-w` kreport, because HUMAnN only sees `s__` lines and pushing reads down to species is precisely what Bracken does. With `--skip_bracken` it falls back to Kraken2's own report, and the functional table gets thinner in proportion to how much of the profile sat above species. This is the same ordering the FGCZ HUMAnN SUSHI app applies.
- **Species-level assignments are what count.** HUMAnN keys on `s__` lines. Clades that stay at genus or above contribute nothing to the functional table, so a profile that is mostly higher-rank will produce a thin result. That is a property of the classification, not a failure.
- **Use the species-named ChocoPhlAn, `v201901_v31`.** HUMAnN 3.6 refuses any profile without an mpa database-version comment line, so the pipeline writes `#mpa_v30_CHOCOPhlAn_201901` at the top of the translated profile. That declaration has to be true of the database you point `--humann_nucleotide_db` at: HUMAnN builds the per-sample nucleotide database by matching each `g__`/`s__` pair against ChocoPhlAn's filenames, which are species-named up to `v201901_v31` and SGB-named (`t__SGB…`) from `vJan21` onwards. NCBI species names — all a Kraken2 or Bracken report can offer — match the former and can never match the latter, so an SGB ChocoPhlAn selects nothing and HUMAnN silently degrades to translated search alone.

`--humann_regroup` (default `uniref90_ko`) and `--humann_renorm` (default `cpm`) control the regrouping and normalisation; raw HUMAnN output is in RPK, which is not comparable across libraries of different depth. Both outputs are also what the bundled [exploreMetaTax](../apps/exploreMetaTax/) app reads in its Gene families and Stratified tabs.

## Loading the results into exploreMetaTax

`--exploremetatax_bundle` writes `exploremetatax/reanatax.exploreMetaTax.tar.gz`
for the [exploreMetaTax](https://shiny-public.fgcz.uzh.ch/app/exploreMetaTax)
Shiny app. To use it:

1. Extract the archive.
2. In the app, open **Upload and Filter Data**.
3. Tick **Browse a folder (uploads all matching reports inside)**.
4. Choose the input format on the radio button, Browse to the extracted folder,
   and press **Load Data**.

The same folder serves every format: the app uploads everything in it and keeps
only the files whose names match the chosen format, so switch the radio button
and load again to move between Kraken2, Bracken and the combined table.

### Why this is not just a tar of `--outdir`

The app decides what a file is from its **name**, and reanaTax's published names
do not line up with what it looks for:

| published name | what the app does with it |
| --- | --- |
| `<sample>.bracken_S.tsv` | matches no format - silently skipped, though this is the table the app actually wants |
| `<sample>.bracken_S.kraken2.report_bracken.txt` | matches `bracken`, but it is a kreport, so the parser looks for `new_est_reads` columns that are not there |
| `bracken_combined_S.txt` | misses `combined_bracken`, whose pattern wants "combined" *before* "bracken" |
| `kraken2_combined_report.txt` | false-matches per-sample `kraken2` and loads as an extra sample |

So the bundle renames files on the way in, leaves out the two that can only
mislead, and gunzips the HUMAnN tables. Before the archive is sealed, every name
is checked against every one of the app's patterns and must match exactly the
one format it is meant for - a rename that would reintroduce an ambiguity fails
the run instead of shipping an archive that loads wrongly.

Whatever ran gets included: Kraken2 reports, per-sample Bracken tables, the
combined Bracken table, MetaPhlAn profiles and the merged table, and the HUMAnN
pathabundance / genefamilies / KO tables. `--da_metadata` is added as
`metadata.tsv` if given - load that through the app's separate Metadata box. A
`README.txt` in the archive repeats the steps above.

The HUMAnN `_renorm` table is deliberately excluded: what it contains depends on
`--humann_regroup`, so its name alone cannot say whether it holds KOs or gene
families.

## Read accounting and diagnostics

Every run writes `read_accounting/reanatax.read_accounting.tsv`: one row per sample tracking raw → trimmed → aligned to each host reference → non-host → classified. All of it exists somewhere already; the table is where it finally sits side by side, which is what you need to decide whether a sample has enough non-host reads for its profile to mean anything. Two columns also land in MultiQC's general statistics:

- **Host carry-over** — the share of supposedly non-host reads that Kraken2 still assigns to the host taxon (`--host_carryover_taxid`, default 9606). Kraken2's standard databases contain the host genome deliberately, which makes this a free audit of the aligner. A few percent is normal. Tens of percent means depletion is leaking, and the fix is usually a second host reference.
- **Non-host** — what fraction of the raw library survived to classification.

`--run_polya_check` adds a diagnostic for poly(A)-selected libraries only. Most bacterial transcripts have no poly(A) tail, so microbial reads in such a library need an explanation: oligo-dT internal mispriming, or non-specific carry-over. The check counts internal poly-A/T runs in the non-host reads against a null built by shuffling each read while preserving its length and base composition — the comparison has to be composition-matched, because A/T runs are common in AT-rich genomes for entirely trivial reasons. An enrichment near 1 means carry-over; clearly above 1 means capture, which biases abundances towards A-rich transcripts.

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

## Differential abundance

Which taxa differ between two groups of samples. Modelled on FGCZ's DiffShot apps (`EzAppDiffShot` in ezRun), reduced to what this pipeline already produces: the combined Bracken table is already a taxa-by-sample count matrix, so the `kraken-biom` round-trip DiffShot needs is skipped.

Group membership comes from `--da_metadata`, a TSV whose first column is the sample id and whose remaining columns are factors:

```tsv
SampleID	Condition	Sex
SRX31575908	case	F
SRX31575925	control	F
```

It is a separate file rather than a samplesheet column because the `--input_accessions` and `--input_dir` routes have no samplesheet. The ids must match the ones in the Bracken table — with the default `--group_runs_by experiment` those are `SRX`/`ERX` accessions, not `SRR`/`ERR`.

```bash
--da_metadata metadata.tsv \
--da_grouping Condition --da_sample_group case --da_ref_group control \
--da_method ancombc,aldex2 --da_covariates Sex
```

| Method | What it does |
| --- | --- |
| `ancombc` | ANCOM-BC2: bias-corrected log-ratios, estimates a per-sample sampling fraction, handles structural zeros. The default. |
| `aldex2` | ALDEx2: centred log-ratio space, Monte-Carlo instances from a Dirichlet posterior. `--da_covariates` switches it from the t-test path to `aldex.glm`. |

Each has its own container — no biocontainer carries both — so naming both runs two processes over the same table. They disagree in informative ways; treat a taxon only one of them calls with suspicion.

**Both need integer counts**, which is why they read the combined Bracken table and why MetaPhlAn relative abundances are rejected rather than silently rescaled. They also run on the *unfiltered* combined table on purpose: each applies its own `--da_prevalence_min` / `--da_relabund_min`, and passing a table already thinned by `--min_rel_abundance` would filter twice under two different rules.

Outputs land in `differential_abundance/`: a results TSV per method (`taxon`, `lfc`, `pvalue`, `qvalue`, `significant`), a volcano plot, and a top-25 table that appears as a MultiQC section.

## MetaPhlAn

`--run_metaphlan` profiles the non-host reads with MetaPhlAn as well. The two classifiers answer the same question differently: Kraken2 assigns every read by k-mer and Bracken redistributes what it could only place above species, whereas MetaPhlAn maps against clade-specific marker genes only and never sees most reads. Where they disagree on one library is worth a look.

```bash
--run_metaphlan \
--metaphlan_db /srv/GT/databases/metaphlan_databases \
--metaphlan_index mpa_vJan25_CHOCOPhlAnSGB_202503
```

It is additive rather than exclusive: add `--skip_kraken2` to use MetaPhlAn *instead* of Kraken2, or leave Kraken2 on to run both over one set of reads. The invocation mirrors ezRun's `app-metaphlan.R`, which is what the databases under `/srv/GT/databases/metaphlan_databases` are laid out for.

`--metaphlan_args` defaults to `-t rel_ab_w_read_stats`, which adds the `estimated_number_of_reads_from_the_clade` column. Per-sample profiles and one merged table land in `metaphlan/`.

Differential abundance does not read MetaPhlAn output — see above for why.

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

# PRISM for reanaTax

reanaTax runs [PRISM](https://github.com/sjdlabgroup/PRISM) (Ghaddar et al.,
_Cancer Cell_ 2026) as an external tool: `--run_prism --prism_path <dir>`. PRISM
is not part of reanaTax and is not distributed with it. Its licence (Rutgers'
non-commercial research licence, `LICENSE.PRISM`) limits it to non-commercial
use, and reanaTax is MIT. You obtain PRISM from its authors and accept that
licence yourself. Check with your technology-transfer office before using it in
fee-for-service work for commercial clients.

What lives here is what reanaTax needs PRISM to be:

- `prism_3cf9e43_reanatax.patch`: fixes to PRISM's code, against commit
  `3cf9e43`. It modifies RU-NCRL code and is under that licence, not MIT.
- `prepare_prism.sh`: clones PRISM at that commit, applies the patch, unpacks
  the GenBank annotation, builds the accession map and puts NCBI's
  `taxonomy4blast.sqlite3` beside the BLAST database.
- `LICENSE.PRISM`: PRISM's licence, which the patch carries.

## Preparing a PRISM directory

```bash
assets/prism/prepare_prism.sh /path/to/PRISM /path/to/blast/core_nt /path/to/genbank.zip
```

`genbank.zip` is the authors' processed GenBank annotation (about 2 GB). It is
downloaded by hand from the Google Drive link in PRISM's README. The accession
map is built from the BLAST database you will give `--prism_blast_db`. For
core_nt this takes hours and writes about 3 GB, and it must be rebuilt when that
database changes.

PRISM's tools (Rscript with xgboost 1.7, kraken2, seqkit, STAR, minimap2,
blastn) come from the PATH, `--prism_bin_dir`, `-profile conda` or
`--prism_container`; see the PRISM section of `docs/usage.md`.

## What the patch changes, and why

Each change was found running PRISM inside reanaTax. The first two made PRISM
fail silently; all apply to upstream's code as of `3cf9e43`.

| File                                  | Change                                                                                                         | Why                                                                                                                                                                                                                                                                                                                                                                                  |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `functions.R`, `prism_custom_blastdb` | `join -t $'\t'` → `-t '\t'`                                                                                    | R's `system()` runs `/bin/sh`, which on Debian is dash. Dash does not expand `$'\t'`, so the join produced an empty taxid map, the custom BLAST database had no taxonomy, every hit came back `staxids = 0`, and PRISM reported 0 taxa.                                                                                                                                              |
| `functions.R`, `prism_kraken`         | Mate suffixes stripped before `/1` and `/2` are re-appended, and the ID list made unique                       | IDs that already ended in `/1` became `…/1/1`, and seqkit extracted no reads.                                                                                                                                                                                                                                                                                                        |
| `functions.R`, `prism_kraken`         | The Kraken2 filter checks taxids as a set (awk lookup) instead of one alternation regex of every taxid to drop | With thousands of taxids to drop, the regex was very slow. Same reads kept.                                                                                                                                                                                                                                                                                                          |
| `functions.R`, `prism_subsample`      | `if (length(idx) == 0) next`                                                                                   | A lineage with no Bacteria, Fungi or Viruses anchor (2, 4751, 10239) made `idx:length(lin)` fail, and the step stopped.                                                                                                                                                                                                                                                              |
| `functions.R`, `prism_genbank`        | `res[, .SD[which.min(interval_width)], by = id]` → `res[res[, .I[which.min(interval_width)], by = id]$V1]`     | Same rows. The `.SD` form ran more than 12 h on a mock-community library and the task hit its time limit.                                                                                                                                                                                                                                                                            |
| `functions.R`, `prism_filter`         | BLAST taxids rolled up to their species (new `prism_species_taxids`, from the Kraken2 database's `nodes.dmp`)  | PRISM counts a read as unique when its hits carry one taxid, and keeps a hit only if its taxid is in the sample's Kraken2 report. Many nt entries are labelled with a strain or serovar, so reads hitting two strains of one species were not unique to it, and hits to a strain Kraken2 had not reported were dropped as non-microbial. Strains also came out as rows of their own. |
| `PRISM.R`                             | Finds `nodes.dmp` in `--kraken_db_path` (or its `taxonomy/`) and passes it to the two filtering steps          | Without one, taxids stay as BLAST reports them, as upstream.                                                                                                                                                                                                                                                                                                                         |
| `PRISM.R`                             | A relative `--model_org_taxids` path is made absolute                                                          | `prism_blast` changes directory (`setwd`) and reads the file again on its second call, so a relative path broke the full BLAST.                                                                                                                                                                                                                                                      |
| `kreport2mpa.py`                      | Made executable                                                                                                | `prism_kraken` runs it as `<prism_path>/kreport2mpa.py`, not through python.                                                                                                                                                                                                                                                                                                         |

Not in the patch, because reanaTax does it in its module (`modules/local/prism/run`):

- **Read IDs only in FASTQ headers.** PRISM looks for `/1` anywhere in the
  header line. An ENA header with the mate number in its comment
  (`@ERR2935854.1 HWI-…:2249/1`) made it append suffixes the IDs do not have, so
  no read was extracted and PRISM reported no taxa.
- **Kraken2 with `--memory-mapping`**, through a wrapper, so concurrent PRISM
  tasks share one copy of the database instead of loading it each.
- **Decoys from the run's host** (`--prism_decoy_taxids`) instead of PRISM's
  bundled model organisms.
- **blastn in slices** (`bin/prism_blastn.py`), through a wrapper: the reads
  are cut into contiguous slices, one blastn per four cores searches each, and
  the outputs are joined in order. A single megablast kept about 8.5 of 16
  threads busy, and `-mt_mode 1` was 2.5 times slower than that. On 200,000
  reads, 4 slices of 4 threads took 183 s against 280 s for one blastn of 16,
  with the same best score, species and per-read uniqueness for every read;
  only which of several equally good genomes fill the ten hit slots changes,
  as it already does between thread counts.
- **STAR with `--runThreadN`** (the task's cpus), through a wrapper. PRISM
  calls STAR without it, so on one thread: about 30 min per 10M read pairs.
  On 100,000 simulated human pairs plus 200,000 microbial ones, 1 and 16
  threads gave the same 214,050 alignment records. On the S. enterica isolate
  ERR2935848 (12.2M pairs) STAR took 62 min on one thread and 3 min on 32.
- **`BLASTDB` set to the BLAST database's directory**, so BLAST finds
  `taxonomy4blast.sqlite3` there for nt and for the per-sample database PRISM
  builds in its own directory.

## taxonomy4blast.sqlite3

PRISM filters BLAST by taxid twice: `blastdbcmd -taxidlist` extracts the
per-sample database, and `blastn -taxids` restricts the searches. BLAST includes
a taxid's descendants only when it finds NCBI's `taxonomy4blast.sqlite3`, which
comes in `taxdb.tar.gz` but is not in every copy of the databases, and BLAST
2.16 looks for it on `BLASTDB`, not beside the database it searches. Without it
both match the exact taxid. In core*nt that leaves out 89% of the Salmonella
enterica sequences, which carry a subspecies or serovar taxid, and 20-27% of
those of \_S. cerevisiae* and _C. neoformans_; for most bacteria it is 1-12%.
With the file, the strains come in, and the species roll-up in the patch is
what keeps them from splitting a species' reads.

## The accession map

PRISM's README builds `sorted_accession_map.txt` with `blastdbcmd -outfmt "%T %a"`
sorted on the first field. The join in `prism_custom_blastdb` reads it as
`taxid<TAB>accession` sorted on the accession, which is what
`prepare_prism.sh` writes. Built the README's way, join rejects the file as
unsorted. Note that `blastdbcmd` writes a literal `\t` when given one inside
`-outfmt`; the script passes a real tab.

## Kept as upstream on purpose

`nreads = 10` in `PRISM.R`'s call to `prism_multimapping`: the reads a species
needs, among the sampled ones, to count as uniquely identifiable. A copy at
FGCZ had it lowered to 2 for a one-off smoke test on 2026-10-25, and it stayed
that way until 2026-10-28. PRISM runs from that clone in between used the looser
threshold.

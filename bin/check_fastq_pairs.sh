#!/usr/bin/env bash
# check_fastq_pairs.sh LOG
#
# LOG is the stderr of a paired `samtools fastq -1 -2 -s /dev/null` run on a
# host aligner's both-mates-unmapped stream. Each pair arrives as two adjacent
# records, and samtools pairs them by identical read name, so a sane run sets
# none aside as singletons. When the mates' names differ (a header such as
# `@SRR1.1/1 1 length=100`, where HISAT2 and Bowtie2 keep the /1 because it is
# not at the end of the line) every read is a singleton and the non-host FASTQs
# come out empty with exit 0. This stops the task instead.
set -euo pipefail

log=$1
cat "$log" >&2

singletons=0
reads=0
while IFS= read -r line; do
    if [[ $line =~ discarded\ ([0-9]+)\ singletons ]]; then
        singletons=${BASH_REMATCH[1]}
    elif [[ $line =~ processed\ ([0-9]+)\ reads ]]; then
        reads=${BASH_REMATCH[1]}
    fi
done < "$log"

if (( singletons > 0 && singletons * 100 > reads )); then
    cat >&2 <<EOF
ERROR: samtools fastq set aside ${singletons} of the ${reads} non-host reads as unpaired, so they are missing from the non-host FASTQs.
The two mates of a pair reached it under different names. HISAT2 and Bowtie2 drop a /1 or /2 only at the very end of the header line, so a header like '@SRR1.1/1 1 length=100' keeps it in the read name.
Remove the /1 and /2 from the first word of each header and rerun, e.g.
  zcat in_1.fastq.gz | sed -E '1~4s#^(@[^[:space:]]+)/[12]([[:space:]]|\$)#\1\2#' | gzip > fixed_1.fastq.gz
EOF
    exit 1
fi

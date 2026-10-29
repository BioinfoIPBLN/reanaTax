#!/usr/bin/env bash
# Prepare a directory for reanaTax's --prism_path: sjdlabgroup/PRISM at the
# commit reanaTax's fixes were made against, those fixes applied, the authors'
# processed GenBank annotation unpacked, the accession-to-taxid map built from
# your own BLAST database, and NCBI's taxonomy4blast.sqlite3 put beside that
# database. See README.md beside this script.
#
# Usage: prepare_prism.sh <target_dir> <blast_db_prefix> <genbank.zip>
#   target_dir       where the prepared PRISM goes (cloned if it does not exist)
#   blast_db_prefix  the BLAST nucleotide database --prism_blast_db will name,
#                    e.g. /data/blast/core_nt
#   genbank.zip      the processed GenBank files, downloaded by hand from the
#                    link in PRISM's README (Google Drive, ~2 GB)
# Needs git, unzip, sort, curl, tar, and blastdbcmd on the PATH.
set -euo pipefail
[[ $# -eq 3 ]] || { sed -n '2,14p' "$0"; exit 1; }
target=$1
blast_db=$2
genbank_zip=$3
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
commit=3cf9e43

command -v blastdbcmd >/dev/null || { echo "blastdbcmd is not on the PATH." >&2; exit 1; }
[[ -r ${genbank_zip} ]] || { echo "Cannot read ${genbank_zip}." >&2; exit 1; }

if [[ ! -d ${target}/.git ]]; then
    git clone -q https://github.com/sjdlabgroup/PRISM "${target}"
fi
git -C "${target}" -c advice.detachedHead=false checkout -q "${commit}"
if git -C "${target}" apply --check "${here}/prism_3cf9e43_reanatax.patch" 2>/dev/null; then
    git -C "${target}" apply "${here}/prism_3cf9e43_reanatax.patch"
    echo "fixes applied"
elif git -C "${target}" apply --reverse --check "${here}/prism_3cf9e43_reanatax.patch" 2>/dev/null; then
    echo "fixes already applied"
else
    echo "${target} has local changes the patch does not expect; inspect 'git -C ${target} diff'." >&2
    exit 1
fi

[[ -s ${target}/genbank/genbank_indices.RDS ]] || unzip -q -o "${genbank_zip}" -d "${target}"
[[ -s ${target}/genbank/genbank_indices.RDS ]] || { echo "No genbank/genbank_indices.RDS after unzipping." >&2; exit 1; }

# BLAST matches a taxid's descendants only with this file where it can find it
# (PRISM_RUN points BLASTDB at the database's directory); without it, a species
# taxid misses every sequence labelled with a strain or serovar.
blast_dir=$(dirname "${blast_db}")
if [[ ! -s ${blast_dir}/taxonomy4blast.sqlite3 ]]; then
    curl -fsSL https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz | tar xz -C "${blast_dir}" taxonomy4blast.sqlite3
fi

# taxid<TAB>accession, sorted by accession: the order PRISM's join reads it in.
# PRISM's README (`%T %a`, sorted on the taxid) gives a map that join rejects.
if [[ ! -s ${target}/sorted_accession_map.txt ]]; then
    blastdbcmd -db "${blast_db}" -entry all -outfmt $'%T\t%a' \
        | LC_ALL=C sort -S 8G -k2,2 -u > "${target}/sorted_accession_map.txt.part"
    mv "${target}/sorted_accession_map.txt.part" "${target}/sorted_accession_map.txt"
fi
echo "ready: --prism_path ${target} --prism_blast_db ${blast_db}"

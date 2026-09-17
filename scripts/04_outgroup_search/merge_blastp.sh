#!/bin/bash
# Check that every query chunk has a finished blastp result and merge them.
# Input : ${PROJECT_ROOT}/blast_results/query_chunks/chunk_<NNN>.faa
#         ${PROJECT_ROOT}/blast_results/chunk_results/chunk_<NNN>.tsv
# Output: ${PROJECT_ROOT}/blast_results/blastp_nr.tsv
# Run from the repository root:
#   bash scripts/04_outgroup_search/merge_blastp.sh

set -euo pipefail
set -a; source config/config.env; set +a

DIR=$PROJECT_ROOT/blast_results

missing=()
results=()
for query in "$DIR"/query_chunks/chunk_*.faa; do
    result=$DIR/chunk_results/$(basename "$query" .faa).tsv
    if [[ -f "$result" ]]; then
        results+=("$result")
    else
        missing+=("$(basename "$query" .faa)")
    fi
done

if (( ${#missing[@]} )); then
    echo "ERROR: no finished result for: ${missing[*]} (resubmit those array tasks)" >&2
    exit 1
fi

cat "${results[@]}" > "$DIR/blastp_nr.tsv"
echo "Merged ${#results[@]} chunks: $(wc -l < "$DIR/blastp_nr.tsv") hits," \
     "$(cut -f1 "$DIR/blastp_nr.tsv" | sort -u | wc -l) of $(cat "$DIR"/query_chunks/*.faa | grep -c '^>') queries with hits"
echo "wrote $DIR/blastp_nr.tsv"

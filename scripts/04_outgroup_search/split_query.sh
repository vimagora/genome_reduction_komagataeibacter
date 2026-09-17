#!/bin/bash
# Split the genes of interest into $BLAST_CHUNKS FASTA files of (nearly) equal
# size, one per blastp array task.
# Input : ${PROJECT_ROOT}/region/genes_of_interest.faa
# Output: ${PROJECT_ROOT}/blast_results/query_chunks/chunk_001.faa ...
# Run from the repository root:
#   bash scripts/04_outgroup_search/split_query.sh

set -euo pipefail
set -a; source config/config.env; set +a

QUERY=$PROJECT_ROOT/region/genes_of_interest.faa
OUTDIR=$PROJECT_ROOT/blast_results/query_chunks
[[ -f "$QUERY" ]] || { echo "ERROR: $QUERY not found (run 03_region first)" >&2; exit 1; }

rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

n_seqs=$(grep -c '^>' "$QUERY")
n_chunks=$(( BLAST_CHUNKS < n_seqs ? BLAST_CHUNKS : n_seqs ))
per_chunk=$(( (n_seqs + n_chunks - 1) / n_chunks ))

# Consecutive blocks of sequences, so the merged output keeps the query order
awk -v per="$per_chunk" -v dir="$OUTDIR" '
/^>/ { if (n++ % per == 0) { if (out) close(out); out = sprintf("%s/chunk_%03d.faa", dir, ++c) } }
{ print > out }
' "$QUERY"

echo "Split $n_seqs sequences into $(ls "$OUTDIR" | wc -l) chunks of up to $per_chunk in $OUTDIR"

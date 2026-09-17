#!/bin/bash
# blastp of one query chunk per array task against $BLASTP_DB, restricted to
# $BLASTP_TAXIDS, reporting taxonomy IDs of the hits.
# Input : ${PROJECT_ROOT}/blast_results/query_chunks/chunk_<NNN>.faa (task N = chunk N)
# Output: ${PROJECT_ROOT}/blast_results/chunk_results/chunk_<NNN>.tsv
# Submit from the repository root, one task per chunk:
#   mkdir -p logs/blast
#   N=$(ls "$PROJECT_ROOT"/blast_results/query_chunks/chunk_*.faa | wc -l)
#   sbatch --account="$SLURM_ACCOUNT" --array=1-${N} scripts/04_outgroup_search/blastp_nr.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=blastp
#SBATCH --output=logs/blast/blast_%A_%a.out
#SBATCH --error=logs/blast/blast_%A_%a.err
#SBATCH --partition=small
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail
set -a; source config/config.env; set +a
# Environment modules often reference unset variables, so relax -u while loading them
if [[ -n "$BLAST_SETUP" ]]; then
    set +u; eval "$BLAST_SETUP"; set -u
fi
command -v blastp >/dev/null 2>&1 || { echo "ERROR: blastp not on PATH (check BLAST_SETUP)" >&2; exit 127; }

CHUNK=$(printf "chunk_%03d" "$SLURM_ARRAY_TASK_ID")
QUERY=$PROJECT_ROOT/blast_results/query_chunks/$CHUNK.faa
OUTDIR=$PROJECT_ROOT/blast_results/chunk_results
OUT=$OUTDIR/$CHUNK.tsv
mkdir -p "$OUTDIR"
[[ -f "$QUERY" ]] || { echo "ERROR: $QUERY not found (run split_query.sh first)" >&2; exit 1; }

TAXID_ARGS=()
[[ -n "$BLASTP_TAXIDS" ]] && TAXID_ARGS=(-taxids "$BLASTP_TAXIDS")

echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID -> $CHUNK ($(grep -c '^>' "$QUERY") queries)"
echo "  db     : $BLASTP_DB"
echo "  taxids : ${BLASTP_TAXIDS:-all}"
echo "  cpus   : $SLURM_CPUS_PER_TASK"

# Write to a temporary file so an interrupted task never leaves a partial result
blastp -db "$BLASTP_DB" \
  -query "$QUERY" \
  -out "$OUT.tmp" \
  "${TAXID_ARGS[@]}" \
  -num_threads "$SLURM_CPUS_PER_TASK" \
  -max_target_seqs 500 \
  -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore staxids stitle"
mv "$OUT.tmp" "$OUT"

echo "[$(date '+%F %T')] done."

#!/bin/bash
# OrthoFinder on the included proteomes.
# Input : ${PROJECT_ROOT}/orthofinder/proteins/
# Output: ${PROJECT_ROOT}/orthofinder/run_<timestamp>/Results_<date>/
#         ${PROJECT_ROOT}/orthofinder/results -> link to that Results folder (read by later steps)
# Submit from the repository root:
#   mkdir -p logs/orthofinder
#   sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/orthofinder.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=orthofinder
#SBATCH --output=logs/orthofinder/orthofinder_%j.out
#SBATCH --error=logs/orthofinder/orthofinder_%j.err
#SBATCH --partition=small
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/orthofinder-env/bin:$PATH"

INPUT=$PROJECT_ROOT/orthofinder/proteins
OUTDIR=$PROJECT_ROOT/orthofinder/run_$(date +%Y%m%d_%H%M%S)
TMPDIR_JOB=$TMPDIR/orthofinder_${SLURM_JOB_ID}

mkdir -p "$TMPDIR_JOB"
export TMPDIR="$TMPDIR_JOB"
trap 'rm -rf "$TMPDIR_JOB"' EXIT

# Preflight
[[ -d "$INPUT" ]] || { echo "ERROR: input dir '$INPUT' not found." >&2; exit 127; }
n_faa=$(ls "$INPUT"/*.faa 2>/dev/null | wc -l)
(( n_faa > 0 )) || { echo "ERROR: no .faa files in $INPUT" >&2; exit 127; }
command -v orthofinder >/dev/null 2>&1 || { echo "ERROR: orthofinder not on PATH." >&2; exit 127; }

echo "[$(date '+%F %T')] Starting OrthoFinder"
echo "  input   : $INPUT ($n_faa proteomes)"
echo "  output  : $OUTDIR"
echo "  cpus    : $SLURM_CPUS_PER_TASK"
echo "  tmpdir  : $TMPDIR"

# -t: parallel sequence-search threads
# -a: parallel analysis threads (usually cpus/4)
# -M msa: use MSA for gene trees (more accurate than default dendroblast)
# -S diamond: use DIAMOND for sequence search (fast, accurate, default anyway)
orthofinder \
    -f "$INPUT" \
    -o "$OUTDIR" \
    -t "$SLURM_CPUS_PER_TASK" \
    -a $(( SLURM_CPUS_PER_TASK / 4 )) \
    -S diamond \
    -M msa

RESULTS=$(find "$OUTDIR" -mindepth 1 -maxdepth 1 -type d -name 'Results_*' | sort | tail -1)
[[ -n "$RESULTS" ]] || { echo "ERROR: no Results_* folder in $OUTDIR" >&2; exit 1; }
ln -sfn "$RESULTS" "$PROJECT_ROOT/orthofinder/results"

echo
echo "[$(date '+%F %T')] OrthoFinder done."
echo "  results: $PROJECT_ROOT/orthofinder/results -> $RESULTS"

#!/bin/bash
# GTDB-Tk classify_wf: GTDB taxonomy of the included genomes.
# Input : ${PROJECT_ROOT}/gtdbtk/input/
# Output: ${PROJECT_ROOT}/gtdbtk/results/classify/ (gtdbtk.bac120.summary.tsv is read next)
# Submit from the repository root:
#   mkdir -p logs/gtdbtk
#   sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/gtdbtk_classify.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=gtdbtk_classify
#SBATCH --output=logs/gtdbtk/gtdbtk_classify_%j.out
#SBATCH --error=logs/gtdbtk/gtdbtk_classify_%j.err
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/gtdbtk-env/bin:$PATH"

INPUT=$PROJECT_ROOT/gtdbtk/input
OUTBASE=$PROJECT_ROOT/gtdbtk/results
TMPDIR_JOB=$TMPDIR/gtdbtk_${SLURM_JOB_ID}

mkdir -p "$OUTBASE" "$TMPDIR_JOB"
export TMPDIR="$TMPDIR_JOB"
trap 'rm -rf "$TMPDIR_JOB"' EXIT

# Preflight
[[ -d "$INPUT" ]] || { echo "ERROR: input dir '$INPUT' not found." >&2; exit 127; }
n_fna=$(ls "$INPUT"/*.fna 2>/dev/null | wc -l)
(( n_fna > 0 )) || { echo "ERROR: no .fna files in $INPUT" >&2; exit 127; }
command -v gtdbtk >/dev/null 2>&1 || { echo "ERROR: gtdbtk not on PATH." >&2; exit 127; }

echo "[$(date '+%F %T')] Running classify_wf"
echo "  input   : $INPUT ($n_fna genomes)"
echo "  db      : $GTDBTK_DATA_PATH"
echo "  cpus    : $SLURM_CPUS_PER_TASK"

gtdbtk classify_wf \
    --genome_dir "$INPUT" \
    --extension fna \
    --out_dir "$OUTBASE/classify" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --tmpdir "$TMPDIR"

echo "[$(date '+%F %T')] classify_wf done."

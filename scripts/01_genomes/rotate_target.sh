#!/bin/bash
# Rotate every contig of the target's assembly with dnaapler, so that the
# chromosome starts at dnaA and each plasmid at repA.
# Input : $TARGET_ASSEMBLY (target = the sample in $REGION)
# Output: ${PROJECT_ROOT}/dnaapler/<target>/<target>_reoriented.fasta (read by annotate_target.sh)
# Run from the repository root, in an interactive session or as a job:
#   bash scripts/01_genomes/rotate_target.sh
#   sbatch --account="$SLURM_ACCOUNT" scripts/01_genomes/rotate_target.sh    (mkdir -p logs/dnaapler first)
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=dnaapler
#SBATCH --output=logs/dnaapler/dnaapler_%j.out
#SBATCH --error=logs/dnaapler/dnaapler_%j.err
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/dnaapler-env/bin:$PATH"

TARGET=$(awk -F'\t' 'NR==2 {print $1}' "$REGION")
OUTDIR=$PROJECT_ROOT/dnaapler/$TARGET
THREADS=${SLURM_CPUS_PER_TASK:-4}

[[ -f "$TARGET_ASSEMBLY" ]] || { echo "ERROR: TARGET_ASSEMBLY '$TARGET_ASSEMBLY' not found." >&2; exit 127; }
command -v dnaapler >/dev/null 2>&1 || { echo "ERROR: dnaapler not on PATH." >&2; exit 127; }

echo "[$(date '+%F %T')] Rotating $TARGET"
echo "  input  : $TARGET_ASSEMBLY"
echo "  output : $OUTDIR"

dnaapler all \
    --input "$TARGET_ASSEMBLY" \
    --output "$OUTDIR" \
    --prefix "$TARGET" \
    --threads "$THREADS" \
    --force

ROTATED=$OUTDIR/${TARGET}_reoriented.fasta
[[ -f "$ROTATED" ]] || { echo "ERROR: $ROTATED not written; see $OUTDIR" >&2; exit 1; }

echo "[$(date '+%F %T')] done: $ROTATED"
grep '^>' "$ROTATED"

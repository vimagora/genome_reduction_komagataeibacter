#!/bin/bash
# BUSCO completeness of one genome per array task (task N = data row N of $GENOMES).
# Submit from the repository root, one task per row of $GENOMES:
#   mkdir -p logs/busco
#   sbatch --account="$SLURM_ACCOUNT" --array=1-$(($(wc -l < config/genomes.tsv) - 1))%8 scripts/02_genome_qc/busco_array.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=busco_array
#SBATCH --output=logs/busco/busco_%A_%a.out
#SBATCH --error=logs/busco/busco_%A_%a.err
#SBATCH --partition=small
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/busco-env/bin:$PATH"

INDIR=$PROJECT_ROOT/genomes/fna
OUTBASE=$PROJECT_ROOT/busco
mkdir -p "$OUTBASE" "$BUSCO_DOWNLOADS"

# Row for this array task (skip header, drop blank lines)
LINE=$(tail -n +2 "$GENOMES" | awk -F'\t' 'NF' | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [[ -z "$LINE" ]]; then
    echo "ERROR: no row in $GENOMES for task $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

SAMPLE=$(echo  "$LINE" | cut -f1)
LINEAGE=$(echo "$LINE" | cut -f10)
INFILE="$INDIR/$SAMPLE.fna"

echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID -> $SAMPLE"
echo "  input   : $INFILE"
echo "  output  : $OUTBASE/$SAMPLE"
echo "  lineage : $LINEAGE"

if [[ ! -f "$INFILE" ]]; then
    echo "ERROR: input '$INFILE' not found." >&2
    exit 127
fi
if [[ -z "$LINEAGE" ]]; then
    echo "ERROR: no busco_lineage set for $SAMPLE." >&2
    exit 1
fi
command -v busco >/dev/null 2>&1 || { echo "ERROR: busco not found on PATH." >&2; exit 127; }

cd "$OUTBASE"
busco \
    --in "$INFILE" \
    --out "$SAMPLE" \
    --mode genome \
    --lineage_dataset "$LINEAGE" \
    --cpu "$SLURM_CPUS_PER_TASK" \
    --download_path "$BUSCO_DOWNLOADS" \
    --force

echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID done."

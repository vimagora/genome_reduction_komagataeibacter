#!/bin/bash
# RANGER-DTL reconciliation of one orthogroup per array task
# (task N = the Nth file of ranger_input/, sorted by name).
# Input : ${PROJECT_ROOT}/hgt_analysis/ranger_input/<OG>.input
# Output: ${PROJECT_ROOT}/hgt_analysis/ranger_output/<OG>.out
# Submit from the repository root, one task per input file:
#   mkdir -p logs/ranger
#   N=$(ls "$PROJECT_ROOT"/hgt_analysis/ranger_input/OG*.input | wc -l)
#   sbatch --account="$SLURM_ACCOUNT" --array=1-${N}%20 scripts/06_hgt/ranger_dtl_array.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=ranger_dtl
#SBATCH --output=logs/ranger/ranger_%A_%a.out
#SBATCH --error=logs/ranger/ranger_%A_%a.err
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

set -uo pipefail
set -a; source config/config.env; set +a

BASE=$PROJECT_ROOT/hgt_analysis
INPUT_DIR="$BASE/ranger_input"
OUTPUT_DIR="$BASE/ranger_output"
mkdir -p "$OUTPUT_DIR"

# Pick this task's input file (1-indexed sort of ranger_input files)
INPUT_FILE=$(ls "$INPUT_DIR"/OG*.input | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
if [[ -z "$INPUT_FILE" ]]; then
    echo "No input file at index $SLURM_ARRAY_TASK_ID" >&2
    exit 0
fi

OG=$(basename "$INPUT_FILE" .input)
OUT="$OUTPUT_DIR/${OG}.out"

echo "[$(date '+%F %T')] task $SLURM_ARRAY_TASK_ID -> $OG"
# $RANGER_ARGS is left unquoted on purpose: it holds zero or more extra options
"$RANGER_DTL" -i "$INPUT_FILE" -o "$OUT" $RANGER_ARGS
echo "[$(date '+%F %T')] done."

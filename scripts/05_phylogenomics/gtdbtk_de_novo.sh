#!/bin/bash
# GTDB-Tk de_novo_wf: species tree of the included genomes, rooted on $GTDBTK_OUTGROUP.
# Inputs: ${PROJECT_ROOT}/gtdbtk/input/, ${PROJECT_ROOT}/gtdbtk/custom_taxonomy.tsv
# Output: ${PROJECT_ROOT}/gtdbtk/results/de_novo/ (gtdbtk.bac120.decorated.tree is read next)
# Submit from the repository root:
#   mkdir -p logs/gtdbtk
#   sbatch --account="$SLURM_ACCOUNT" scripts/05_phylogenomics/gtdbtk_de_novo.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=gtdbtk_de_novo
#SBATCH --output=logs/gtdbtk/gtdbtk_de_novo_%j.out
#SBATCH --error=logs/gtdbtk/gtdbtk_de_novo_%j.err
#SBATCH --partition=small
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/gtdbtk-env/bin:$PATH"

INPUT=$PROJECT_ROOT/gtdbtk/input
OUTBASE=$PROJECT_ROOT/gtdbtk/results
TAXONOMY=$PROJECT_ROOT/gtdbtk/custom_taxonomy.tsv
TMPDIR_JOB=$TMPDIR/gtdbtk_${SLURM_JOB_ID}

mkdir -p "$OUTBASE" "$TMPDIR_JOB"
export TMPDIR="$TMPDIR_JOB"
trap 'rm -rf "$TMPDIR_JOB"' EXIT

# Preflight
[[ -d "$INPUT" ]] || { echo "ERROR: input dir '$INPUT' not found." >&2; exit 127; }
n_fna=$(ls "$INPUT"/*.fna 2>/dev/null | wc -l)
(( n_fna > 0 )) || { echo "ERROR: no .fna files in $INPUT" >&2; exit 127; }
[[ -f "$TAXONOMY" ]] || { echo "ERROR: $TAXONOMY not found (run make_custom_taxonomy.sh)" >&2; exit 127; }
command -v gtdbtk >/dev/null 2>&1 || { echo "ERROR: gtdbtk not on PATH." >&2; exit 127; }

echo "[$(date '+%F %T')] Running de_novo_wf"
echo "  input    : $INPUT ($n_fna genomes)"
echo "  db       : $GTDBTK_DATA_PATH"
echo "  outgroup : $GTDBTK_OUTGROUP"
echo "  cpus     : $SLURM_CPUS_PER_TASK"

gtdbtk de_novo_wf \
    --genome_dir "$INPUT" \
    --extension fna \
    --bacteria \
    --out_dir "$OUTBASE/de_novo" \
    --cpus "$SLURM_CPUS_PER_TASK" \
    --outgroup_taxon "$GTDBTK_OUTGROUP" \
    --skip_gtdb_refs \
    --custom_taxonomy_file "$TAXONOMY" \
    --tmpdir "$TMPDIR"

echo "[$(date '+%F %T')] de_novo_wf done."

#!/bin/bash
# Annotate the rotated target assembly with Bakta and store the genome, proteins
# and annotation where the downloaded genomes are:
#   ${PROJECT_ROOT}/genomes/fna/<target>.fna
#   ${PROJECT_ROOT}/genomes/faa/<target>.faa
#   ${PROJECT_ROOT}/genomes/gff/<target>.gff
# Input : ${PROJECT_ROOT}/dnaapler/<target>/<target>_reoriented.fasta (from rotate_target.sh)
#         genus, species and strain of the target from $GENOMES (target = the sample in $REGION)
# Output: ${PROJECT_ROOT}/bakta_annot/<target>/ (all Bakta output) and the three copies above
# Submit from the repository root:
#   mkdir -p logs/bakta
#   sbatch --account="$SLURM_ACCOUNT" scripts/01_genomes/annotate_target.sh
# SLURM does not expand variables in #SBATCH lines; --account is overridden at submission.
#SBATCH --account=${SLURM_ACCOUNT:-your_account}
#SBATCH --job-name=bakta
#SBATCH --output=logs/bakta/bakta_%j.out
#SBATCH --error=logs/bakta/bakta_%j.err
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/bakta-env/bin:$PATH"

TARGET=$(awk -F'\t' 'NR==2 {print $1}' "$REGION")
INFILE=$PROJECT_ROOT/dnaapler/$TARGET/${TARGET}_reoriented.fasta
OUTDIR=$PROJECT_ROOT/bakta_annot/$TARGET
GENOMES_DIR=$PROJECT_ROOT/genomes

# Genus, species and strain of the target, by column name
read -r GENUS SPECIES STRAIN < <(awk -F'\t' -v t="$TARGET" '
NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
$1 == t { print $col["genus"], $col["species"], $col["strain"]; found = 1 }
END { if (!found) exit 1 }' "$GENOMES") \
    || { echo "ERROR: $TARGET not found in $GENOMES" >&2; exit 1; }

[[ -f "$INFILE" ]] || { echo "ERROR: $INFILE not found (run rotate_target.sh first)" >&2; exit 127; }
command -v bakta >/dev/null 2>&1 || { echo "ERROR: bakta not on PATH." >&2; exit 127; }

echo "[$(date '+%F %T')] Annotating $TARGET ($GENUS $SPECIES $STRAIN)"
echo "  input  : $INFILE"
echo "  output : $OUTDIR"
echo "  db     : $BAKTA_DB"

mkdir -p "$OUTDIR"
bakta \
    --db "$BAKTA_DB" \
    --output "$OUTDIR" \
    --prefix "$TARGET" \
    --genus "$GENUS" \
    --species "$SPECIES" \
    --strain "$STRAIN" \
    --threads "$SLURM_CPUS_PER_TASK" \
    --force \
    --compliant \
    --keep-contig-headers \
    "$INFILE"

# Store the three files like a downloaded genome
mkdir -p "$GENOMES_DIR/fna" "$GENOMES_DIR/faa" "$GENOMES_DIR/gff"
cp "$OUTDIR/$TARGET.fna"  "$GENOMES_DIR/fna/$TARGET.fna"
cp "$OUTDIR/$TARGET.faa"  "$GENOMES_DIR/faa/$TARGET.faa"
cp "$OUTDIR/$TARGET.gff3" "$GENOMES_DIR/gff/$TARGET.gff"

echo "[$(date '+%F %T')] done."
printf "  %s\n" "$GENOMES_DIR/fna/$TARGET.fna" "$GENOMES_DIR/faa/$TARGET.faa" "$GENOMES_DIR/gff/$TARGET.gff"
echo "  proteins: $(grep -c '^>' "$GENOMES_DIR/faa/$TARGET.faa")"

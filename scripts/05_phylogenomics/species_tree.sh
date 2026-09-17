#!/bin/bash
# Convert the GTDB-Tk species tree: an iTOL version for viewing, and a plain
# newick tree (labels removed) used by the HGT analysis and the heatmap.
# Input : ${PROJECT_ROOT}/gtdbtk/results/de_novo/gtdbtk.bac120.decorated.tree
# Output: ${PROJECT_ROOT}/hgt_analysis/gtdbtk_itol.tree
#         ${PROJECT_ROOT}/hgt_analysis/species_tree.newick
# Run from the repository root:
#   bash scripts/05_phylogenomics/species_tree.sh

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/gtdbtk-env/bin:$PATH"

TREE=$PROJECT_ROOT/gtdbtk/results/de_novo/gtdbtk.bac120.decorated.tree
OUTDIR=$PROJECT_ROOT/hgt_analysis
mkdir -p "$OUTDIR"

gtdbtk convert_to_itol --input_tree "$TREE" --output_tree "$OUTDIR/gtdbtk_itol.tree"
gtdbtk remove_labels   --input_tree "$TREE" --output_tree "$OUTDIR/species_tree.newick"

echo "wrote $OUTDIR/gtdbtk_itol.tree and $OUTDIR/species_tree.newick"

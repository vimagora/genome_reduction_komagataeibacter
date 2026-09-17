#!/bin/bash
# Build the custom taxonomy file for GTDB-Tk de_novo_wf from the classify_wf
# results. Genomes that GTDB-Tk could not assign to a species (empty "s__")
# get "s__<genus> sp", using the genus of their own classification.
# Input : ${PROJECT_ROOT}/gtdbtk/results/classify/gtdbtk.bac120.summary.tsv
# Output: ${PROJECT_ROOT}/gtdbtk/custom_taxonomy.tsv
# Run from the repository root:
#   bash scripts/05_phylogenomics/make_custom_taxonomy.sh

set -euo pipefail
set -a; source config/config.env; set +a

CLASSIFY=$PROJECT_ROOT/gtdbtk/results/classify/gtdbtk.bac120.summary.tsv
OUT=$PROJECT_ROOT/gtdbtk/custom_taxonomy.tsv

awk -F'\t' 'BEGIN {OFS="\t"}
NR > 1 {
    taxonomy = $2
    if (taxonomy ~ /s__$/ && match(taxonomy, /g__[^;]+/)) {
        taxonomy = taxonomy substr(taxonomy, RSTART + 3, RLENGTH - 3) " sp"
        filled++
    }
    print $1, taxonomy
}
END { print "species filled from genus: " filled+0 > "/dev/stderr" }' "$CLASSIFY" > "$OUT"

echo "wrote $OUT ($(wc -l < "$OUT") genomes)"

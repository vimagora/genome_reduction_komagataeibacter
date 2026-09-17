#!/bin/bash
# Find the orthogroups of the genes of interest and collect their gene trees
# and gene counts.
# Inputs : ${PROJECT_ROOT}/region/genes_of_interest.txt
#          ${PROJECT_ROOT}/orthofinder/results/ (Orthogroups.tsv, Orthogroups.GeneCount.tsv, Gene_Trees/)
# Outputs (${PROJECT_ROOT}/hgt_analysis/):
#          genes.txt, gene_to_orthogroup.tsv, orthogroups_of_interest.txt,
#          gene_trees/<OG>_tree.txt, gene_count_matrix.tsv
# Run from the repository root:
#   bash scripts/06_hgt/extract_hgt_orthogroups.sh

set -uo pipefail
set -a; source config/config.env; set +a

# --- inputs -------------------------------------------------------------------
TARGET=$(awk -F'\t' 'NR==2 {print $1}' "$REGION")
GENES=$PROJECT_ROOT/region/genes_of_interest.txt
RESULTS=$PROJECT_ROOT/orthofinder/results
OUTDIR=$PROJECT_ROOT/hgt_analysis

OGTSV="$RESULTS/Orthogroups/Orthogroups.tsv"
GCTSV="$RESULTS/Orthogroups/Orthogroups.GeneCount.tsv"
GT="$RESULTS/Gene_Trees"

# --- preflight ----------------------------------------------------------------
for f in "$GENES" "$OGTSV" "$GCTSV"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: $f not found" >&2
        exit 1
    fi
done
if [[ ! -d "$GT" ]]; then
    echo "ERROR: $GT dir not found" >&2
    exit 1
fi

mkdir -p "$OUTDIR/gene_trees"

# --- prep: unique gene list ---------------------------------------------------
tr -d '\r' < "$GENES" | awk 'NF' | sort -u > "$OUTDIR/genes.txt"
n_genes=$(wc -l < "$OUTDIR/genes.txt")
echo "Target: $TARGET"
echo "Genes of interest: $n_genes"

# --- Step 1: map genes -> orthogroups -----------------------------------------
# Only the target's column is searched: NCBI protein IDs (WP_*) can be shared
# by several genomes.
echo "Mapping genes to orthogroups..."
awk -F'\t' -v genes="$OUTDIR/genes.txt" -v out="$OUTDIR/gene_to_orthogroup.tsv" -v target="$TARGET" '
BEGIN {
    while ((getline g < genes) > 0) if (g != "") goi[g] = 1
    close(genes)
    printf "gene\torthogroup\n" > out
}
NR == 1 {
    for (i = 2; i <= NF; i++) if ($i == target) col = i
    if (!col) { print "ERROR: " target " is not a column of Orthogroups.tsv" > "/dev/stderr"; exit 1 }
    next
}
{
    n = split($col, cell, ", ")
    for (j = 1; j <= n; j++) {
        gene = cell[j]
        if (gene in goi) mapping[gene] = $1
    }
}
END {
    if (!col) exit 1
    n = asorti(goi, sorted)
    for (k = 1; k <= n; k++) {
        g = sorted[k]
        if (g in mapping) printf "%s\t%s\n", g, mapping[g] >> out
        else              printf "%s\t<UNASSIGNED>\n", g >> out
    }
    n_assigned = 0
    for (g in mapping) n_assigned++
    print "  assigned to orthogroups: " n_assigned
    print "  unassigned (species-specific): " (length(goi) - n_assigned)
}' "$OGTSV" || exit 1

# --- Step 2: unique orthogroups of interest -----------------------------------
awk -F'\t' 'NR>1 && $2 != "<UNASSIGNED>" {print $2}' "$OUTDIR/gene_to_orthogroup.tsv" \
    | sort -u > "$OUTDIR/orthogroups_of_interest.txt"
n_ogs=$(wc -l < "$OUTDIR/orthogroups_of_interest.txt")
echo "Unique orthogroups involving the genes of interest: $n_ogs"

# --- Step 3: copy the gene tree for each orthogroup ---------------------------
missing=0
copied=0
while read -r og; do
    src="$GT/${og}_tree.txt"
    if [[ -f "$src" ]]; then
        cp "$src" "$OUTDIR/gene_trees/"
        copied=$((copied+1))
    else
        missing=$((missing+1))
    fi
done < "$OUTDIR/orthogroups_of_interest.txt"
echo "Gene trees copied: $copied  (missing $missing — OGs too small for a tree)"

# --- Step 4: gene count matrix ------------------------------------------------
head -1 "$GCTSV" > "$OUTDIR/gene_count_matrix.tsv"
grep -F -w -f "$OUTDIR/orthogroups_of_interest.txt" "$GCTSV" \
    >> "$OUTDIR/gene_count_matrix.tsv"
n_rows=$(( $(wc -l < "$OUTDIR/gene_count_matrix.tsv") - 1 ))
echo "Gene count matrix rows: $n_rows"

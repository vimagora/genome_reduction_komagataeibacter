#!/bin/bash
# Keep the transfers into the target whose donor is a genome outside the
# target's own taxon at $DONOR_EXCLUDE_RANK (e.g. Family); transfers from
# ancestral nodes (n1, n2, ...) are not kept. Then list the genes of interest in
# those orthogroups with their NCBI annotation.
# Inputs : ${PROJECT_ROOT}/hgt_analysis/hgt_candidates_detailed.tsv, gene_to_orthogroup.tsv
#          ${PROJECT_ROOT}/region/region_cds.tsv, $GENOMES, $REGION
# Outputs (${PROJECT_ROOT}/hgt_analysis/):
#          strong_hgt_ogs_with_donors.tsv  orthogroup, donor genome
#          hgt_candidate_ogs.txt           orthogroups with a distant donor
#          hgt_candidate_genes.tsv         orthogroup, gene
#          hgt_candidate_annotated.tsv     orthogroup, gene, locus tag, product
# Run from the repository root:
#   bash scripts/06_hgt/filter_distant_donors.sh

set -euo pipefail
set -a; source config/config.env; set +a

DIR=$PROJECT_ROOT/hgt_analysis
TARGET=$(awk -F'\t' 'NR==2 {print $1}' "$REGION")

# Genomes sharing the target's taxon at $DONOR_EXCLUDE_RANK
awk -F'\t' -v rank="$DONOR_EXCLUDE_RANK" -v target="$TARGET" '
NR == 1 { for (i = 1; i <= NF; i++) if ($i == rank) col = i; next }
{ taxon[$1] = $col }
END {
    if (!col) { print "ERROR: no column " rank " in genomes table" > "/dev/stderr"; exit 1 }
    for (s in taxon) if (taxon[s] == taxon[target]) print s
}' "$GENOMES" > "$DIR/excluded_donor_genomes.txt"
echo "Donors excluded ($DONOR_EXCLUDE_RANK of $TARGET): $(wc -l < "$DIR/excluded_donor_genomes.txt") genomes"

# Distant donors: named genomes (not ancestral nodes) outside the excluded set
awk -F'\t' '
NR == FNR { excluded[$1] = 1; next }
FNR > 1 && $5 != "" && !($5 in excluded) && $5 !~ /^n[0-9]+$/ { print $1 "\t" $5 }
' "$DIR/excluded_donor_genomes.txt" "$DIR/hgt_candidates_detailed.tsv" \
    | sort -u > "$DIR/strong_hgt_ogs_with_donors.tsv"

cut -f1 "$DIR/strong_hgt_ogs_with_donors.tsv" | sort -u > "$DIR/hgt_candidate_ogs.txt"

# Genes of interest in the candidate orthogroups
awk -F'\t' '
NR == FNR { ogs[$1] = 1; next }
$2 in ogs { print $2 "\t" $1 }
' "$DIR/hgt_candidate_ogs.txt" "$DIR/gene_to_orthogroup.tsv" \
    | sort > "$DIR/hgt_candidate_genes.tsv"

# Add locus tag and product from the region's NCBI annotation
{
    echo -e "orthogroup\tgene_id\tlocus_tag\tproduct"
    awk -F'\t' '
    NR == FNR { if (FNR > 1 && !($6 in locus)) { locus[$6] = $5; product[$6] = $7 }; next }
    { print $1 "\t" $2 "\t" locus[$2] "\t" product[$2] }
    ' "$PROJECT_ROOT/region/region_cds.tsv" "$DIR/hgt_candidate_genes.tsv" | sort
} > "$DIR/hgt_candidate_annotated.tsv"

echo "Transfers from distant donors: $(wc -l < "$DIR/strong_hgt_ogs_with_donors.tsv") (orthogroup, donor) pairs"
echo "Candidate orthogroups:         $(wc -l < "$DIR/hgt_candidate_ogs.txt")"
echo "Candidate genes:               $(wc -l < "$DIR/hgt_candidate_genes.tsv")"
echo "wrote $DIR/hgt_candidate_annotated.tsv"

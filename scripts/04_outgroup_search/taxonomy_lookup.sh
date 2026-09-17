#!/bin/bash
# Look up the taxonomic lineage of every taxon hit by blastp.
# Input : ${PROJECT_ROOT}/blast_results/blastp_nr.tsv
# Output: ${PROJECT_ROOT}/blast_results/unique_taxids.txt
#         ${PROJECT_ROOT}/blast_results/taxid_lineage_lookup.tsv
#           (taxid, full lineage, phylum, class, order, family, genus, species)
# Run from the repository root:
#   bash scripts/04_outgroup_search/taxonomy_lookup.sh

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/ncbi-tools-env/bin:$PATH"

DIR=$PROJECT_ROOT/blast_results

# Unique taxids of the blastp hits (column 11; several taxids are separated by ';')
cut -f11 "$DIR/blastp_nr.tsv" | tr ';' '\n' | awk 'NF' | sort -u > "$DIR/unique_taxids.txt"

taxonkit lineage --data-dir "$TAXONKIT_DB" "$DIR/unique_taxids.txt" \
    | taxonkit reformat --data-dir "$TAXONKIT_DB" -I 1 -f "{p}\t{c}\t{o}\t{f}\t{g}\t{s}" \
    > "$DIR/taxid_lineage_lookup.tsv"

echo "taxids: $(wc -l < "$DIR/unique_taxids.txt")"
echo "wrote $DIR/taxid_lineage_lookup.tsv"

#!/bin/bash
# Look up the taxonomic lineage of every taxon hit by blastp.
# Input : ${PROJECT_ROOT}/blast_results/blastp_nr.tsv
# Output: ${PROJECT_ROOT}/blast_results/unique_taxids.txt
#         ${PROJECT_ROOT}/blast_results/taxid_lineage_lookup.tsv
#           (taxid, full lineage, phylum, class, order, family, genus, species)
#         ${PROJECT_ROOT}/blast_results/taxonkit_warnings.log
# Taxids merged since the BLAST database was built are resolved to their new
# taxid; deleted taxids get an empty lineage (and are dropped by add_taxonomy.py).
# Run from the repository root:
#   bash scripts/04_outgroup_search/taxonomy_lookup.sh

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/ncbi-tools-env/bin:$PATH"

DIR=$PROJECT_ROOT/blast_results

# Unique taxids of the blastp hits (column 11; several taxids are separated by ';')
cut -f11 "$DIR/blastp_nr.tsv" | tr ';' '\n' | awk 'NF' | sort -u > "$DIR/unique_taxids.txt"

# Both taxonkit commands warn about merged/deleted taxids; keep the warnings in a log
taxonkit lineage --data-dir "$TAXONKIT_DB" "$DIR/unique_taxids.txt" 2> "$DIR/taxonkit_warnings.log" \
    | taxonkit reformat --data-dir "$TAXONKIT_DB" -I 1 -f "{p}\t{c}\t{o}\t{f}\t{g}\t{s}" 2>> "$DIR/taxonkit_warnings.log" \
    > "$DIR/taxid_lineage_lookup.tsv" \
    || { echo "ERROR: taxonkit failed; last messages:" >&2; tail -5 "$DIR/taxonkit_warnings.log" >&2; exit 1; }

n_merged=$( { grep -o 'taxid [0-9]* was merged' "$DIR/taxonkit_warnings.log" || true; } | sort -u | wc -l)
n_deleted=$( { grep -o 'taxid [0-9]* was deleted' "$DIR/taxonkit_warnings.log" || true; } | sort -u | wc -l)
echo "taxids: $(wc -l < "$DIR/unique_taxids.txt") ($n_merged merged into a newer taxid, $n_deleted deleted from the taxonomy)"
echo "warnings: $DIR/taxonkit_warnings.log"
echo "wrote $DIR/taxid_lineage_lookup.tsv"

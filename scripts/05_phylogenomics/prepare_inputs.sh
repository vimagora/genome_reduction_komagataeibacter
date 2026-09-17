#!/bin/bash
# Link the proteomes and genomes of the genomes marked include=yes in $GENOMES
# into the OrthoFinder and GTDB-Tk input folders.
# Inputs : ${PROJECT_ROOT}/genomes/faa/<sample>.faa, ${PROJECT_ROOT}/genomes/fna/<sample>.fna
# Outputs: ${PROJECT_ROOT}/orthofinder/proteins/<sample>.faa (links)
#          ${PROJECT_ROOT}/gtdbtk/input/<sample>.fna (links)
# Run from the repository root:
#   bash scripts/05_phylogenomics/prepare_inputs.sh

set -euo pipefail
set -a; source config/config.env; set +a

FAA_IN=$PROJECT_ROOT/genomes/faa
FNA_IN=$PROJECT_ROOT/genomes/fna
FAA_OUT=$PROJECT_ROOT/orthofinder/proteins
FNA_OUT=$PROJECT_ROOT/gtdbtk/input
mkdir -p "$FAA_OUT" "$FNA_OUT"

# Clear links from a previous selection so it doesn't leak into this run
find "$FAA_OUT" "$FNA_OUT" -maxdepth 1 -type l -delete

n_ok=0
n_miss=0
while read -r sample; do
    faa=$FAA_IN/$sample.faa
    fna=$FNA_IN/$sample.fna
    if [[ -f "$faa" && -f "$fna" ]]; then
        ln -sf "$faa" "$FAA_OUT/$sample.faa"
        ln -sf "$fna" "$FNA_OUT/$sample.fna"
        n_ok=$((n_ok+1))
    else
        printf "  MISS  %-45s (proteins or genome not downloaded)\n" "$sample"
        n_miss=$((n_miss+1))
    fi
done < <(awk -F'\t' 'NR>1 && NF && $11=="yes" {print $1}' "$GENOMES")

echo "Summary: $n_ok linked, $n_miss missing"
echo "  proteomes: $(find "$FAA_OUT" -maxdepth 1 -type l | wc -l)  ($FAA_OUT)"
echo "  genomes:   $(find "$FNA_OUT" -maxdepth 1 -type l | wc -l)  ($FNA_OUT)"

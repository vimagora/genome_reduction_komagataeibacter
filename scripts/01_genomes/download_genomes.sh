#!/bin/bash
# Download every genome listed in $GENOMES from NCBI (genome sequence, proteins
# and GFF3 annotation) and store the files by sample name:
#   ${PROJECT_ROOT}/genomes/fna/<sample>.fna
#   ${PROJECT_ROOT}/genomes/faa/<sample>.faa
#   ${PROJECT_ROOT}/genomes/gff/<sample>.gff
#
# Needs internet access (run it on a login node), from the repository root:
#   bash scripts/01_genomes/download_genomes.sh

set -euo pipefail
set -a; source config/config.env; set +a
export PATH="$ENV_ROOT/ncbi-tools-env/bin:$PATH"

OUT=$PROJECT_ROOT/genomes
mkdir -p "$OUT/fna" "$OUT/faa" "$OUT/gff"

# accession <tab> sample, skipping rows whose accession is still a <TODO> placeholder
awk -F'\t' 'NR>1 && NF {
    if ($2 ~ /^</) print "WARN: no accession for " $1 " - skipped" > "/dev/stderr"
    else           print $2 "\t" $1
}' "$GENOMES" > "$OUT/accession_to_sample.tsv"
cut -f1 "$OUT/accession_to_sample.tsv" > "$OUT/accessions.txt"
echo "Accessions to download: $(wc -l < "$OUT/accessions.txt")"

datasets download genome accession \
    --inputfile "$OUT/accessions.txt" \
    --include genome,protein,gff3 \
    --filename "$OUT/ncbi_download.zip"
rm -rf "$OUT/ncbi_download"
unzip -q "$OUT/ncbi_download.zip" -d "$OUT/ncbi_download"
DL=$OUT/ncbi_download/ncbi_dataset/data

n_ok=0
n_incomplete=0
while IFS=$'\t' read -r accession sample; do
    accdir=$DL/$accession
    genome=$(find "$accdir" -maxdepth 1 -name 'GC[AF]_*_genomic.fna' 2>/dev/null | head -1 || true)
    if [[ -n "$genome" ]];             then cp "$genome"             "$OUT/fna/$sample.fna"; fi
    if [[ -f "$accdir/protein.faa" ]]; then cp "$accdir/protein.faa" "$OUT/faa/$sample.faa"; fi
    if [[ -f "$accdir/genomic.gff" ]]; then cp "$accdir/genomic.gff" "$OUT/gff/$sample.gff"; fi

    missing=()
    [[ -f "$OUT/fna/$sample.fna" ]] || missing+=(genome)
    [[ -f "$OUT/faa/$sample.faa" ]] || missing+=(proteins)
    [[ -f "$OUT/gff/$sample.gff" ]] || missing+=(annotation)
    if (( ${#missing[@]} )); then
        printf "  MISS  %-45s %s: no %s\n" "$sample" "$accession" "${missing[*]}"
        n_incomplete=$((n_incomplete+1))
    else
        printf "  OK    %-45s %s (%d proteins)\n" "$sample" "$accession" "$(grep -c '^>' "$OUT/faa/$sample.faa")"
        n_ok=$((n_ok+1))
    fi
done < "$OUT/accession_to_sample.tsv"

echo
echo "Summary: complete=$n_ok incomplete=$n_incomplete"

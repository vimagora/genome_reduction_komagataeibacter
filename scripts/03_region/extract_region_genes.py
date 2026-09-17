#!/usr/bin/env python3
"""Define the genes of the region of interest from the target genome's NCBI annotation.

A CDS belongs to the region when its midpoint lies between the start and end
given in $REGION (the same rule as the permutation test). Transposases and
insertion sequences (products matching $IS_PATTERN), pseudogenes, and CDS
without a sequence in the protein FASTA are excluded from the genes of interest.

Inputs:
  $REGION, $GENOMES
  ${PROJECT_ROOT}/genomes/gff/<target>.gff
  ${PROJECT_ROOT}/genomes/faa/<target>.faa
Outputs (${PROJECT_ROOT}/region/):
  region_cds.tsv            every CDS in the region, with exclusion flags
  genes_of_interest.txt     protein IDs of the genes of interest, in genome order
  genes_of_interest.faa     their protein sequences

Run from the repository root, after loading config/config.env:
  python3 scripts/03_region/extract_region_genes.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import in_window, mobile_element_pattern, read_cds, read_region, work_path  # noqa: E402

sample, seqid, start, end = read_region()
gff = work_path("genomes", "gff", f"{sample}.gff")
faa = work_path("genomes", "faa", f"{sample}.faa")
outdir = work_path("region")
outdir.mkdir(parents=True, exist_ok=True)

mobile = mobile_element_pattern()
region = [c for c in read_cds(gff) if in_window(c, seqid, start, end)]
if not region:
    sys.exit(f"ERROR: no CDS found in {seqid}:{start}-{end} of {gff}")
wanted = {c["protein_id"] for c in region}

# Protein sequences of the region's CDS
sequences, current = {}, None
with open(faa) as fh:
    for line in fh:
        if line.startswith(">"):
            pid = line[1:].split()[0]
            current = pid if pid in wanted else None
            if current:
                sequences[current] = [line]
        elif current:
            sequences[current].append(line)

genes, seen, missing = [], set(), []
with open(outdir / "region_cds.tsv", "w") as fh:
    fh.write("seqid\tstart\tend\tstrand\tlocus_tag\tprotein_id\tproduct\tmobile_element\tpseudogene\n")
    for c in region:
        is_mobile = bool(mobile.search(c["product"]))
        fh.write(f"{c['seqid']}\t{c['start']}\t{c['end']}\t{c['strand']}\t{c['locus_tag']}"
                 f"\t{c['protein_id']}\t{c['product']}"
                 f"\t{'yes' if is_mobile else 'no'}\t{'yes' if c['pseudo'] else 'no'}\n")
        if is_mobile or c["pseudo"] or c["protein_id"] in seen:
            continue
        if c["protein_id"] not in sequences:
            missing.append(c["protein_id"] or c["locus_tag"])
            continue
        seen.add(c["protein_id"])
        genes.append(c["protein_id"])

(outdir / "genes_of_interest.txt").write_text("\n".join(genes) + "\n")
with open(outdir / "genes_of_interest.faa", "w") as fh:
    for g in genes:
        fh.writelines(sequences[g])

n_mobile = sum(1 for c in region if mobile.search(c["product"]))
n_pseudo = sum(1 for c in region if c["pseudo"])
print(f"Region: {sample} {seqid}:{start:,}-{end:,} ({end - start + 1:,} bp)")
print(f"  CDS in region:                 {len(region)}")
print(f"  transposases/IS:               {n_mobile}")
print(f"  pseudogenes:                   {n_pseudo}")
print(f"  genes of interest (proteins):  {len(genes)}")
if missing:
    print(f"WARN: {len(missing)} CDS without a sequence in {faa}, excluded: {', '.join(missing[:5])}", file=sys.stderr)
print(f"Wrote {outdir}/region_cds.tsv, genes_of_interest.txt, genes_of_interest.faa")

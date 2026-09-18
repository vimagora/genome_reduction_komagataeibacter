# Deprecated analysis: HGT inference with RANGER-DTL

This analysis is not part of the article. It infers horizontal gene transfer
by reconciling the gene tree of each orthogroup of the region with the species
tree using RANGER-DTL, then keeps transfers from donors outside the target's
family and plots them as a heatmap. The HGT results in the article were
obtained with other methods, described there.

It is kept for reference. It continues from `../data_analysis.md`: run steps
0–6.1 there first, and load the configuration as described in its step 0. The
section numbers below follow that recipe.

## 0. Setup

Install RANGER-DTL

| Tool | Version | Environment |
|---|---|---|
| ete3, pandas, matplotlib | not pinned (Python 3.10) | `envs/ete3.yaml` |
| RANGER-DTL | 2.0 | standalone binary |

RANGER-DTL is a standalone binary; set its path in `RANGER_DTL`.

From 6.1, the orthogroup list is read by 6.2;
`gene_to_orthogroup.tsv` by 6.5 and 6.6.

## 6. Horizontal gene transfer

```bash
mkdir -p logs/ranger
```
### 6.2 RANGER-DTL input

Run `prepare_ranger_input.py` to relabel the species tree and each gene tree with
short species codes, resolve polytomies, and write one RANGER-DTL input file
(species tree + gene tree) per orthogroup:

```bash
"$PY" scripts/06_hgt/prepare_ranger_input.py
```

Output: `hgt_analysis/species_code_map.tsv`, `species_tree_coded.newick` and
`ranger_input/<OG>.input`. The inputs are read by 6.3; the code map by 6.4.

### 6.3 Reconciliation

Run `ranger_dtl_array.sh` to reconcile each gene tree with the species tree
(one array task per input file):

```bash
N=$(ls "$PROJECT_ROOT"/hgt_analysis/ranger_input/OG*.input | wc -l)
sbatch --account="$SLURM_ACCOUNT" --array=1-${N}%20 scripts/06_hgt/ranger_dtl_array.sh
```

Output: `hgt_analysis/ranger_output/<OG>.out`. It is read by 6.4.

### 6.4 Transfer events

Run `parse_ranger_hgt.py` to extract the transfers whose recipient is the target
(direct) or an ancestor of the target:

```bash
"$PY" scripts/06_hgt/parse_ranger_hgt.py
```

Output: `hgt_analysis/hgt_summary_per_og.tsv` (counts per orthogroup) and
`hgt_candidates_detailed.tsv` (one row per transfer). The detailed table is read by 6.5.

### 6.5 Distant donors

Run `filter_distant_donors.sh` to keep the transfers from genomes outside the
target's taxon at `DONOR_EXCLUDE_RANK`, and annotate the genes in those
orthogroups with `region/region_cds.tsv`:

```bash
bash scripts/06_hgt/filter_distant_donors.sh
```

Output: `hgt_analysis/strong_hgt_ogs_with_donors.tsv`, `hgt_candidate_ogs.txt`,
`hgt_candidate_genes.tsv` and `hgt_candidate_annotated.tsv`. The candidate
orthogroups are read by step 6.6.

### 6.6 Heatmap

Run `heatmap_order.py` to order the genomes by the ladderized species tree and
the genes by genome position:

```bash
"$PY" scripts/06_hgt/heatmap_order.py
```

Output: `figures/strain_order.txt` and `figures/gene_order.txt`. They are read by
`plot_hgt_heatmap.py`, together with the gene counts from `orthofinder/results`,
`hgt_analysis/gene_to_orthogroup.tsv`, `hgt_analysis/hgt_candidate_ogs.txt`, the
species tree and the taxonomy in `genomes.tsv`:

```bash
"$PY" scripts/06_hgt/plot_hgt_heatmap.py
```

Output: `figures/hgt_heatmap.png` and `figures/hgt_heatmap.svg`.
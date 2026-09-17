#!/usr/bin/env python3
"""Row and column order for the heatmap.

Rows: genomes in the order of the ladderized species tree.
Columns: genes of interest in genome order (as written by extract_region_genes.py).

Inputs : ${PROJECT_ROOT}/hgt_analysis/species_tree.newick
         ${PROJECT_ROOT}/region/genes_of_interest.txt
Outputs: ${PROJECT_ROOT}/figures/strain_order.txt
         ${PROJECT_ROOT}/figures/gene_order.txt

Run from the repository root, after loading config/config.env:
  python3 scripts/07_figures/heatmap_order.py
"""
import sys
from pathlib import Path

from ete3 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import work_path  # noqa: E402

outdir = work_path("figures")
outdir.mkdir(parents=True, exist_ok=True)

t = Tree(str(work_path("hgt_analysis", "species_tree.newick")), format=1)
t.ladderize()  # ladderize gives visually nice ordering
order = [leaf.name for leaf in t.get_leaves()]
(outdir / "strain_order.txt").write_text("\n".join(order) + "\n")
print(f"wrote {len(order)} strains in tree order")

genes = work_path("region", "genes_of_interest.txt").read_text().split()
(outdir / "gene_order.txt").write_text("\n".join(genes) + "\n")
print(f"wrote {len(genes)} genes in genome order")

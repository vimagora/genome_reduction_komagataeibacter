#!/usr/bin/env python3
"""Count how many genes of interest have at least one blastp hit outside the
target's own family, order, class and phylum (taxonomy taken from $GENOMES).

Inputs : ${PROJECT_ROOT}/blast_results/results_with_taxonomy.tsv
         ${PROJECT_ROOT}/region/genes_of_interest.txt
Output : printed table

Run from the repository root, after loading config/config.env:
  python3 scripts/04_outgroup_search/taxonomy_coverage_check.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import genome_row, read_tsv, cfg, work_path  # noqa: E402

target = read_tsv(cfg("REGION"))[0]["sample"]
lineage = genome_row(target)
n_queries = len(work_path("region", "genes_of_interest.txt").read_text().split())

df = pd.read_csv(work_path("blast_results", "results_with_taxonomy.tsv"), sep="\t")

for level, column in [("family", "Family"), ("order", "Order"),
                      ("class", "Class"), ("phylum", "Phylum")]:
    name = lineage[column]
    outside = df[df[level] != name]
    n_queries_with_hits = outside["qseqid"].nunique()
    print(f"Excluding {level} ({name}): {n_queries_with_hits} / {n_queries} queries have at least one hit outside")

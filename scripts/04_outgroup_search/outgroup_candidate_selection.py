#!/usr/bin/env python3
"""Rank candidate outgroup species: for each taxonomic level of the target
(genus, family, order, class, phylum; taken from $GENOMES), list the species
outside that taxon that match the most genes of interest.

Input  : ${PROJECT_ROOT}/blast_results/results_with_taxonomy.tsv
Output : printed tables (top $OUTGROUP_TOP_N species per level)

Run from the repository root, after loading config/config.env:
  python3 scripts/04_outgroup_search/outgroup_candidate_selection.py
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cfg, genome_row, read_tsv, work_path  # noqa: E402

target = read_tsv(cfg("REGION"))[0]["sample"]
lineage = genome_row(target)
top_n = int(cfg("OUTGROUP_TOP_N"))

df = pd.read_csv(work_path("blast_results", "results_with_taxonomy.tsv"), sep="\t")

# Pattern catches MAG/rank-placeholder names: "X bacterium", "X archaeon", "Candidatus X", "uncultured X", bare "bacterium"
placeholder_pattern = re.compile(
    r"(?i)\bbacterium$|\barchaeon$|^candidatus|^uncultured|^unclassified|^metagenome",
)
df_clean = df[~df["species"].fillna("").str.contains(placeholder_pattern)]

levels = [
    ("genus",  lineage["genus"]),
    ("family", lineage["Family"]),
    ("order",  lineage["Order"]),
    ("class",  lineage["Class"]),
    ("phylum", lineage["Phylum"]),
]

for level, name in levels:
    outside = df_clean[df_clean[level] != name]
    top = outside.groupby("species")["qseqid"].nunique().sort_values(ascending=False).head(top_n)
    print(f"\n=== Top {top_n} outgroup species excluding {level} ({name}) ===")
    print(top.to_string())

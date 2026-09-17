#!/usr/bin/env python3
"""Attach the taxonomic lineage to every blastp hit.

Inputs : ${PROJECT_ROOT}/blast_results/blastp_nr.tsv
         ${PROJECT_ROOT}/blast_results/taxid_lineage_lookup.tsv
Output : ${PROJECT_ROOT}/blast_results/results_with_taxonomy.tsv

Run from the repository root, after loading config/config.env:
  python3 scripts/04_outgroup_search/add_taxonomy.py
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import work_path  # noqa: E402

blast_file = work_path("blast_results", "blastp_nr.tsv")
taxonomy_file = work_path("blast_results", "taxid_lineage_lookup.tsv")

cols = ["qseqid", "sseqid", "pident", "length", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "staxids", "stitle"]
blast = pd.read_csv(blast_file, sep="\t", names=cols, dtype={"staxids": str})

# split staxids on ';' and explode into one row per taxid
blast["staxids"] = blast["staxids"].str.split(";")
blast = blast.explode("staxids").reset_index(drop=True)

lineage = pd.read_csv(
    taxonomy_file, sep="\t",
    names=["taxid", "full_lineage", "phylum", "class", "order", "family", "genus", "species"],
    dtype={"taxid": str}
)
lineage = lineage.drop(columns=["full_lineage"])  # don't need the raw string version

merged = blast.merge(lineage, left_on="staxids", right_on="taxid", how="left")
merged = merged.drop(columns=["taxid"])

merged.to_csv(work_path("blast_results", "results_with_taxonomy.tsv"), sep="\t", index=False)
print(f"Total rows: {len(merged)}")
print(f"Rows with missing phylum: {merged['phylum'].isna().sum()}")

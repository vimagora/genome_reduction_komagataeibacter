#!/usr/bin/env python3
"""Attach the taxonomic lineage to every blastp hit.

Hits whose taxid has no lineage in the taxonomy dump (deleted since the BLAST
database was built) are dropped, so they are not counted as hits outside the
target's taxa.

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
merged = blast.merge(lineage, left_on="staxids", right_on="taxid", how="left")

unresolved = merged["full_lineage"].fillna("").eq("")
print(f"Rows without a lineage (taxid unknown to the taxonomy dump), dropped: {unresolved.sum()}"
      f" ({merged.loc[unresolved, 'staxids'].nunique()} taxids)")
merged = merged[~unresolved]
merged = merged.drop(columns=["taxid", "full_lineage"])  # don't need the raw string version

merged.to_csv(work_path("blast_results", "results_with_taxonomy.tsv"), sep="\t", index=False)
print(f"Total rows: {len(merged)}")
print(f"Rows with missing phylum: {merged['phylum'].isna().sum()}")

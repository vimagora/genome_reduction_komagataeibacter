#!/usr/bin/env python3
"""Prepare RANGER-DTL input files for HGT analysis.

Uses ete3 for tree manipulation. Handles:
  - Species tree preparation: relabel with short codes, resolve polytomies,
    strip everything except topology and leaf names.
  - Gene tree preparation: relabel each leaf as <SpeciesCode>_<index> so
    RANGER-DTL's underscore-based species detection works, resolve polytomies.
  - Persist a species code map for interpreting RANGER outputs later.

Inputs (under ${PROJECT_ROOT}):
  - hgt_analysis/species_tree.newick (from species_tree.sh)
  - orthofinder/results/Orthogroups/Orthogroups.tsv (species column order + gene-species mapping)
  - orthofinder/results/Gene_Trees/ (one .txt per orthogroup)
  - orthogroups_of_interest.txt (from extract_hgt_orthogroups.sh)

Outputs:
  - hgt_analysis/species_code_map.tsv
  - hgt_analysis/species_tree_coded.newick
  - hgt_analysis/ranger_input/<OG>.input (one per OG of interest)

Run from the repository root, after loading config/config.env:
  python3 scripts/06_hgt/prepare_ranger_input.py
"""
from collections import defaultdict
from pathlib import Path
from ete3 import Tree
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import work_path  # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = work_path("hgt_analysis")
RESULTS = work_path("orthofinder", "results")

RAW_SPECIES_TREE = BASE / "species_tree.newick"
OGTSV = RESULTS / "Orthogroups/Orthogroups.tsv"
GT_DIR = RESULTS / "Gene_Trees"
OGS_LIST = BASE / "orthogroups_of_interest.txt"

INPUT_DIR = BASE / "ranger_input"
CODED_SPECIES_TREE = BASE / "species_tree_coded.newick"
CODE_MAP = BASE / "species_code_map.tsv"

INPUT_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Step 1: build species -> short-code map, sorted alphabetically for stability
# ---------------------------------------------------------------------------
with open(OGTSV) as fh:
    species = sorted(fh.readline().rstrip("\n").split("\t")[1:])

name_to_code = {sp: f"S{i+1:02d}" for i, sp in enumerate(species)}

with open(CODE_MAP, "w") as fh:
    fh.write("code\toriginal_species_name\n")
    for name, code in name_to_code.items():
        fh.write(f"{code}\t{name}\n")
print(f"[1] Wrote code map ({len(species)} species) -> {CODE_MAP}")

# ---------------------------------------------------------------------------
# Step 2: prepare species tree — relabel, resolve polytomies, topology-only
# ---------------------------------------------------------------------------
t = Tree(str(RAW_SPECIES_TREE), format=1)  # format=1 handles internal node names
for leaf in t.get_leaves():
    if leaf.name not in name_to_code:
        raise ValueError(f"Species tree leaf '{leaf.name}' not in OrthoFinder headers")
    leaf.name = name_to_code[leaf.name]

t.resolve_polytomy(recursive=True)

# format=9 = leaf names only, no branch lengths, no support, no internal labels
t.write(outfile=str(CODED_SPECIES_TREE), format=9)
species_tree_text = CODED_SPECIES_TREE.read_text().strip()
print(f"[2] Wrote coded species tree ({len(t.get_leaves())} leaves) -> {CODED_SPECIES_TREE}")

# ---------------------------------------------------------------------------
# Step 3: prepare each gene tree
# ---------------------------------------------------------------------------
# OrthoFinder gene tree leaves are named "<SpeciesName>_<geneID>", so we match
# each leaf name against species names (longest-first) to find its species.
species_sorted_len = sorted(name_to_code, key=len, reverse=True)

def leaf_to_species(leaf_name):
    for sp in species_sorted_len:
        if leaf_name == sp or leaf_name.startswith(sp + "_"):
            return sp
    return None

ogs = OGS_LIST.read_text().split()
print(f"[3] Processing {len(ogs)} orthogroups of interest...")

n_ok = 0
n_no_tree = 0
n_no_species_match = 0
failed = []

for og_id in ogs:
    gt_path = GT_DIR / f"{og_id}_tree.txt"
    if not gt_path.exists():
        n_no_tree += 1
        continue

    try:
        # format=0 = flexible parsing, accepts support values and other quirks
        gt = Tree(str(gt_path), format=0)
    except Exception as e:
        failed.append((og_id, f"parse: {e}"))
        continue

    # Rename each leaf to <SpeciesCode>_<index>
    counter = defaultdict(int)
    unmatched = 0
    for leaf in gt.get_leaves():
        sp = leaf_to_species(leaf.name)
        if sp is None:
            unmatched += 1
            continue
        code = name_to_code[sp]
        counter[code] += 1
        leaf.name = f"{code}_{counter[code]}"

    if not counter:
        n_no_species_match += 1
        continue

    gt.resolve_polytomy(recursive=True)

    # Write RANGER-DTL input file: species tree line 1, gene tree line 2
    gene_tree_text = gt.write(format=9)
    (INPUT_DIR / f"{og_id}.input").write_text(
        f"{species_tree_text}\n{gene_tree_text}\n"
    )
    n_ok += 1

print(f"    written: {n_ok}")
print(f"    skipped (no gene tree file): {n_no_tree}")
print(f"    skipped (no leaves matched a species): {n_no_species_match}")
if failed:
    print(f"    failed to parse: {len(failed)}")
    for og, err in failed[:5]:
        print(f"      {og}: {err}")
    if len(failed) > 5:
        print(f"      ... and {len(failed)-5} more")

print(f"\nDone. Input files in: {INPUT_DIR}")
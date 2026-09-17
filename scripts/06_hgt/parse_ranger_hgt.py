#!/usr/bin/env python3
"""Parse all RANGER-DTL outputs, extract transfer events, and identify HGT
candidates involving the target genome (the sample in $REGION).

Inputs : ${PROJECT_ROOT}/hgt_analysis/ranger_output/<OG>.out, species_code_map.tsv
Outputs: ${PROJECT_ROOT}/hgt_analysis/hgt_summary_per_og.tsv, hgt_candidates_detailed.tsv

Run from the repository root, after loading config/config.env:
  python3 scripts/06_hgt/parse_ranger_hgt.py
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cfg, read_tsv, work_path  # noqa: E402

BASE = work_path("hgt_analysis")
OUTPUT_DIR = BASE / "ranger_output"
CODE_MAP = BASE / "species_code_map.tsv"

TARGET_SPECIES = read_tsv(cfg("REGION"))[0]["sample"]

# ---------------------------------------------------------------------------
# Load code map: S## <-> species name
# ---------------------------------------------------------------------------
code_to_name = {}
name_to_code = {}
with open(CODE_MAP) as fh:
    fh.readline()
    for line in fh:
        code, name = line.rstrip("\n").split("\t")
        code_to_name[code] = name
        name_to_code[name] = code

target_code = name_to_code.get(TARGET_SPECIES)
if not target_code:
    sys.exit(f"ERROR: '{TARGET_SPECIES}' not in code map")
print(f"Target: {TARGET_SPECIES} = {target_code}")

# ---------------------------------------------------------------------------
# Parse RANGER species tree into node -> descendant species codes
# ---------------------------------------------------------------------------
def parse_species_tree_topology(newick):
    """Return {node_label: set_of_descendant_species_codes} for internal + leaf nodes."""
    tree = newick.rstrip(";").strip()
    node_desc = {}
    while "(" in tree:
        match = re.search(r'\(([^()]+)\)(\w+)?', tree)
        if not match:
            break
        contents = match.group(1)
        label = match.group(2) or "unnamed"
        parts = [p.strip() for p in contents.split(",")]
        desc = set()
        for p in parts:
            if p in node_desc:
                desc.update(node_desc[p])
            else:
                desc.add(p)
        node_desc[label] = desc
        tree = tree[:match.start()] + label + tree[match.end():]
    return node_desc

def parse_ranger_output(path):
    text = path.read_text()

    # Species tree from RANGER's output (has internal node labels n1, n2, ...)
    m = re.search(r'Species Tree:\s*\n\s*(\S+)', text)
    if not m:
        return None, [], {}
    topology = parse_species_tree_topology(m.group(1))

    transfers = []
    for line in text.splitlines():
        m = re.match(
            r'\s*(m\d+)\s*=\s*LCA\[[^\]]+\]:\s*Transfer,\s*Mapping\s*-->\s*(\S+?),\s*Recipient\s*-->\s*(\S+)',
            line,
        )
        if m:
            transfers.append({
                "gene_node": m.group(1),
                "donor": m.group(2),
                "recipient": m.group(3),
            })

    cost = {}
    m = re.search(r'Duplications:\s*(\d+),\s*Transfers:\s*(\d+),\s*Losses:\s*(\d+)', text)
    if m:
        cost = {"D": int(m.group(1)), "T": int(m.group(2)), "L": int(m.group(3))}

    return topology, transfers, cost

# ---------------------------------------------------------------------------
# Process all OGs
# ---------------------------------------------------------------------------
og_files = sorted(OUTPUT_DIR.glob("OG*.out"))
print(f"\nParsing {len(og_files)} RANGER outputs...")

results_all = []
transfers_to_target = []

for of in og_files:
    og_id = of.stem
    topology, transfers, cost = parse_ranger_output(of)
    if topology is None:
        continue

    tier1 = []
    tier2 = []
    for t in transfers:
        recipient = t["recipient"]
        if recipient == target_code:
            tier1.append(t)
        elif recipient in topology and target_code in topology[recipient]:
            tier2.append(t)

    results_all.append({
        "og": og_id,
        "D": cost.get("D", 0), "T": cost.get("T", 0), "L": cost.get("L", 0),
        "tier1_direct": len(tier1), "tier2_ancestor": len(tier2),
    })
    for t in tier1 + tier2:
        transfers_to_target.append({
            "og": og_id,
            "gene_node": t["gene_node"],
            "donor": t["donor"],
            "recipient": t["recipient"],
            "tier": "direct" if t["recipient"] == target_code else "ancestor",
        })

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------
def resolve(code_or_node):
    return code_to_name.get(code_or_node, code_or_node)

summary_path = BASE / "hgt_summary_per_og.tsv"
with open(summary_path, "w") as fh:
    fh.write("orthogroup\tduplications\ttransfers\tlosses\ttransfers_to_target_direct\ttransfers_to_target_ancestor\n")
    for r in sorted(results_all, key=lambda x: (-x["tier1_direct"], -x["tier2_ancestor"])):
        fh.write(f"{r['og']}\t{r['D']}\t{r['T']}\t{r['L']}\t{r['tier1_direct']}\t{r['tier2_ancestor']}\n")

details_path = BASE / "hgt_candidates_detailed.tsv"
with open(details_path, "w") as fh:
    fh.write("orthogroup\ttier\tgene_node\tdonor_code\tdonor_species\trecipient_code\trecipient_species\n")
    for t in transfers_to_target:
        fh.write(f"{t['og']}\t{t['tier']}\t{t['gene_node']}"
                 f"\t{t['donor']}\t{resolve(t['donor'])}"
                 f"\t{t['recipient']}\t{resolve(t['recipient'])}\n")

print(f"\nWrote: {summary_path}")
print(f"Wrote: {details_path}")

ogs_with_direct = [r for r in results_all if r["tier1_direct"] > 0]
ogs_with_ancestor = [r for r in results_all if r["tier2_ancestor"] > 0]
print(f"\nOrthogroups with direct HGT into {TARGET_SPECIES}: {len(ogs_with_direct)}")
print(f"Orthogroups with HGT into ancestor of {TARGET_SPECIES}: {len(ogs_with_ancestor)}")

if ogs_with_direct:
    print(f"\nTop 10 by direct-transfer count:")
    print(f"{'OG':<12}{'Direct':<10}{'Ancestor':<10}{'Transfers':<12}")
    for r in ogs_with_direct[:10]:
        print(f"{r['og']:<12}{r['tier1_direct']:<10}{r['tier2_ancestor']:<10}{r['T']:<12}")
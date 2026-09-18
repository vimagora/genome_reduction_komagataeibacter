#!/usr/bin/env python3
"""Heatmap: presence/absence of the orthogroups of the region's genes across the
included genomes, ordered by species tree (rows) and genomic position (columns).

Layout (left -> right):
    [ species tree cladogram | taxonomy bars | heatmap ]
Legend goes on the right side of the figure.

Inputs : $GENOMES (taxonomy), ${PROJECT_ROOT}/orthofinder/results/Orthogroups/Orthogroups.GeneCount.tsv,
         ${PROJECT_ROOT}/hgt_analysis/ (species_tree.newick, gene_to_orthogroup.tsv, hgt_candidate_ogs.txt),
         ${PROJECT_ROOT}/figures/ (strain_order.txt, gene_order.txt)
Outputs: ${PROJECT_ROOT}/figures/hgt_heatmap.png, hgt_heatmap.svg

Run from the repository root, after loading config/config.env:
  python3 scripts/06_hgt/plot_hgt_heatmap.py
"""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import numpy as np
from ete3 import Tree
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import cfg, read_tsv, work_path  # noqa: E402

BASE = work_path("hgt_analysis")
FIGURES = work_path("figures")
RESULTS = work_path("orthofinder", "results")
METADATA = Path(cfg("GENOMES"))

TARGET_STRAIN = read_tsv(cfg("REGION"))[0]["sample"]
HGT_OGS_FILE = BASE / "hgt_candidate_ogs.txt"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
strain_order = (FIGURES / "strain_order.txt").read_text().split()
gene_order   = (FIGURES / "gene_order.txt").read_text().split()

gene_to_og = {}
with open(BASE / "gene_to_orthogroup.tsv") as fh:
    fh.readline()
    for line in fh:
        gene, og = line.rstrip("\n").split("\t")
        gene_to_og[gene] = og

hgt_ogs = set(HGT_OGS_FILE.read_text().split())

gc = pd.read_csv(RESULTS / "Orthogroups/Orthogroups.GeneCount.tsv", sep="\t", index_col=0)
gc = gc.drop(columns=["Total"], errors="ignore")

meta = pd.read_csv(METADATA, sep="\t")
meta = meta[["sample", "Phylum", "Class", "Order", "Family"]].set_index("sample")

# ---------------------------------------------------------------------------
# Presence/absence matrix
# ---------------------------------------------------------------------------
gene_ogs = [gene_to_og.get(g) for g in gene_order]
valid_cols = [(g, og) for g, og in zip(gene_order, gene_ogs) if og and og in gc.index]
genes_kept = [x[0] for x in valid_cols]
ogs_kept   = [x[1] for x in valid_cols]

pa = (gc.loc[ogs_kept, strain_order] > 0).astype(int)
pa.columns = strain_order
pa.index = [f"{g} / {og}" for g, og in valid_cols]

M = pa.T
n_rows, n_cols = M.shape
print(f"Matrix: {n_rows} strains x {n_cols} genes")

# ---------------------------------------------------------------------------
# Species tree — precompute coords for the cladogram
# ---------------------------------------------------------------------------
tree = Tree(str(BASE / "species_tree.newick"), format=1)
tree.ladderize()

# Assign y positions to leaves in the order they appear in strain_order
leaf_y = {name: i for i, name in enumerate(strain_order)}

# For internal nodes, y = midpoint of descendant leaves
# Depth from root = x position (rooted left, tips at right)
def compute_coords(node, depth=0):
    if node.is_leaf():
        node._x = depth
        node._y = leaf_y[node.name]
    else:
        for ch in node.children:
            compute_coords(ch, depth + 1)
        node._x = depth
        # y as midpoint of children's y's
        ys = [ch._y for ch in node.children]
        node._y = (min(ys) + max(ys)) / 2
compute_coords(tree)

# Normalize depths to [0, 1] for tree pane
max_depth = max(n._x for n in tree.traverse())
for n in tree.traverse():
    n._x = n._x / max_depth if max_depth else 0

# Collect line segments (parent->child) for drawing
segments = []  # each is (x0, y0, x1, y1)
for node in tree.traverse():
    if node.is_leaf():
        continue
    for ch in node.children:
        # Horizontal segment from parent x to child x at child y
        segments.append((node._x, ch._y, ch._x, ch._y))
    # Vertical segment at parent x from min child y to max child y
    ys = [ch._y for ch in node.children]
    segments.append((node._x, min(ys), node._x, max(ys)))

# ---------------------------------------------------------------------------
# Taxonomy annotations
# ---------------------------------------------------------------------------
tax_levels = ["Phylum", "Class", "Order", "Family"]
tax_df = meta.reindex(strain_order)[tax_levels]

def make_cmap(values):
    unique = list(dict.fromkeys(v for v in values if pd.notna(v) and v))
    palette = plt.cm.tab20(np.linspace(0, 1, max(len(unique), 20)))
    return {v: palette[i % 20] for i, v in enumerate(unique)}, unique

cmaps = {level: make_cmap(tax_df[level].tolist()) for level in tax_levels}

# ---------------------------------------------------------------------------
# Layout dimensions (in inches)
# ---------------------------------------------------------------------------
CELL_W = 0.10                 # cell width (was 0.20)
CELL_H = 0.30                 # cell height (unchanged)

TREE_W = 2.2                  # phylogeny pane
TAX_BAR_W = 0.28              # each taxonomy bar
N_TAX = len(tax_levels)
HM_W = CELL_W * n_cols
HM_H = CELL_H * n_rows

# margins
LEFT = 0.2
RIGHT_LEGEND_W = 3.2          # right-side legend area
TOP = 3.0                     # room for gene labels rotated
BOTTOM = 0.4

FIG_W = LEFT + TREE_W + N_TAX * TAX_BAR_W + 0.1 + HM_W + 0.1 + RIGHT_LEGEND_W
FIG_H = TOP + HM_H + BOTTOM

# Axes x-positions in figure fraction
tree_x0    = LEFT / FIG_W
tree_w     = TREE_W / FIG_W
tax_x0     = (LEFT + TREE_W) / FIG_W
tax_w      = (N_TAX * TAX_BAR_W) / FIG_W
hm_x0      = (LEFT + TREE_W + N_TAX * TAX_BAR_W + 0.1) / FIG_W
hm_w       = HM_W / FIG_W
legend_x0  = (LEFT + TREE_W + N_TAX * TAX_BAR_W + 0.1 + HM_W + 0.1) / FIG_W
legend_w   = RIGHT_LEGEND_W / FIG_W

y0 = BOTTOM / FIG_H
axes_h = HM_H / FIG_H

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(FIG_W, FIG_H))

# --- 1. Species tree cladogram ---
ax_tree = fig.add_axes([tree_x0, y0, tree_w, axes_h])
for (x0v, y0v, x1v, y1v) in segments:
    ax_tree.plot([x0v, x1v], [y0v, y1v], color="black", linewidth=0.8)
ax_tree.set_xlim(-0.02, 1.05)
ax_tree.set_ylim(n_rows - 0.5, -0.5)  # invert so top-of-tree matches heatmap top
ax_tree.set_xticks([])
ax_tree.set_yticks([])
for spine in ax_tree.spines.values():
    spine.set_visible(False)

# --- 2. Taxonomy bars ---
for j, level in enumerate(tax_levels):
    bar_x0 = (LEFT + TREE_W + j * TAX_BAR_W) / FIG_W
    bar_w  = (TAX_BAR_W * 0.9) / FIG_W
    ax = fig.add_axes([bar_x0, y0, bar_w, axes_h])
    color_map, _ = cmaps[level]
    for i, strain in enumerate(strain_order):
        val = tax_df.loc[strain, level]
        color = color_map.get(val, "#DDDDDD")
        ax.add_patch(mpatches.Rectangle((0, i - 0.5), 1, 1, color=color, lw=0))
    ax.set_xlim(0, 1)
    ax.set_ylim(n_rows - 0.5, -0.5)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    # Level name below bar (rotated so it doesn't overlap gene labels above)
    ax.text(0.5, n_rows - 0.5 + 0.7, level, ha="center", va="top",
            fontsize=7, fontweight="bold")

# --- 3. Heatmap ---
ax_hm = fig.add_axes([hm_x0, y0, hm_w, axes_h])
# Colors: white for absent, classic heatmap blue for present
cmap = ListedColormap(["#F7F7F7", "#2166AC"])
ax_hm.imshow(M.values, aspect="auto", cmap=cmap, vmin=0, vmax=1,
             interpolation="nearest")

ax_hm.set_yticks([])  # strains are already labeled by the tree on the left
ax_hm.set_xticks(range(n_cols))
ax_hm.set_xticklabels(M.columns, rotation=90, fontsize=5, ha="center", va="bottom")
ax_hm.xaxis.tick_top()

# Highlight HGT-candidate gene labels
for label in ax_hm.get_xticklabels():
    og = label.get_text().split(" / ")[-1]
    if og in hgt_ogs:
        label.set_color("#C22")
        label.set_fontweight("bold")

# Thin white grid
ax_hm.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
ax_hm.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
ax_hm.grid(which="minor", color="white", linewidth=0.15)
ax_hm.tick_params(which="minor", length=0)

# Highlight target strain: draw a red rectangle around its row
try:
    target_row = strain_order.index(TARGET_STRAIN)
    ax_hm.add_patch(mpatches.Rectangle(
        (-0.5, target_row - 0.5), n_cols, 1,
        fill=False, edgecolor="#C22", linewidth=1.2, zorder=10,
    ))
except ValueError:
    pass

# --- 4. Legend (right side, stacked vertically) ---
ax_leg = fig.add_axes([legend_x0, 0.05, legend_w - 0.02, 0.95])
ax_leg.axis("off")

# Present/absent legend
pa_handles = [
    mpatches.Patch(color="#2166AC", label="Orthogroup present"),
    mpatches.Patch(color="#F7F7F7", label="Orthogroup absent", ec="#999"),
]
leg_pa = ax_leg.legend(handles=pa_handles, title="Heatmap",
                       loc="upper left", bbox_to_anchor=(0, 1.0),
                       fontsize=8, title_fontsize=9, frameon=False)
ax_leg.add_artist(leg_pa)

# Taxonomy legends
y_anchor = 0.90
for level in tax_levels:
    color_map, unique = cmaps[level]
    if not unique:
        continue
    handles = [mpatches.Patch(color=color_map[v], label=v) for v in unique]
    leg = ax_leg.legend(handles=handles, title=level,
                        loc="upper left", bbox_to_anchor=(0, y_anchor),
                        fontsize=7, title_fontsize=8, frameon=False)
    ax_leg.add_artist(leg)
    # Rough spacing: leave ~0.04 per row + title
    y_anchor -= 0.03 * (len(unique) + 2)

# Target strain marker
target_handle = mpatches.Patch(edgecolor="#C22", facecolor="none", label=TARGET_STRAIN)
ax_leg.legend(handles=[target_handle], title="Target",
              loc="upper left", bbox_to_anchor=(0, y_anchor),
              fontsize=7, title_fontsize=8, frameon=False)

# --- Save ---
out_png = FIGURES / "hgt_heatmap.png"
out_svg = FIGURES / "hgt_heatmap.svg"
plt.savefig(out_png, dpi=200, bbox_inches="tight", facecolor="white")
plt.savefig(out_svg, bbox_inches="tight", facecolor="white")
print(f"Wrote:\n  {out_png}\n  {out_svg}")

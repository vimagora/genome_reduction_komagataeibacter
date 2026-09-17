#!/usr/bin/env python3
"""Permutation test: is the region enriched for transposases/IS relative to
random windows of the same size on the same sequence?

A CDS is counted in a window when its midpoint lies inside it; it is a
transposase/IS when its product matches $IS_PATTERN. Random windows never
overlap the region.

Inputs:
  $REGION, $IS_PATTERN, $N_PERMUTATIONS, $PERMUTATION_SEED
  ${PROJECT_ROOT}/genomes/gff/<target>.gff
  ${PROJECT_ROOT}/genomes/fna/<target>.fna
Output:
  ${PROJECT_ROOT}/statistics/permutation_null_distribution.txt

Run from the repository root, after loading config/config.env:
  python3 scripts/08_statistics/permutation_test_transposases.py
"""
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (cfg, in_window, mobile_element_pattern, read_cds,  # noqa: E402
                    read_region, sequence_lengths, work_path)

sample, region_seqid, region_start, region_end = read_region()
region_size = region_end - region_start + 1
N_PERMUTATIONS = int(cfg("N_PERMUTATIONS"))
SEED = int(cfg("PERMUTATION_SEED"))
random.seed(SEED)

print(f"Region: {region_seqid}:{region_start:,} - {region_end:,}  ({region_size:,} bp)")

lengths = sequence_lengths(work_path("genomes", "fna", f"{sample}.fna"))
if region_seqid not in lengths:
    sys.exit(f"ERROR: {region_seqid} not found in the {sample} genome FASTA")
L = lengths[region_seqid]
print(f"Sequence length: {L:,} bp")
if L < 2 * region_size:
    sys.exit(f"ERROR: {region_seqid} ({L:,} bp) too small for non-overlapping windows of {region_size:,} bp")

mobile = mobile_element_pattern()
cds = [(c, bool(mobile.search(c["product"])))
       for c in read_cds(work_path("genomes", "gff", f"{sample}.gff"))
       if c["seqid"] == region_seqid]


def count_te_in_window(start, end):
    inside = [is_te for c, is_te in cds if in_window(c, region_seqid, start, end)]
    return len(inside), sum(inside)


observed_total, observed_te = count_te_in_window(region_start, region_end)

random_counts = []
for _ in range(N_PERMUTATIONS):
    for _ in range(100):
        start = random.randint(1, L - region_size + 1)
        end = start + region_size - 1
        if end < region_start or start > region_end:
            break
    else:
        sys.exit("ERROR: could not find non-overlapping window after 100 attempts")
    _, n_te = count_te_in_window(start, end)
    random_counts.append(n_te)

n_ge = sum(1 for x in random_counts if x >= observed_te)
pval = n_ge / N_PERMUTATIONS
mean_random = sum(random_counts) / N_PERMUTATIONS
fold = observed_te / mean_random if mean_random > 0 else float("inf")

print()
print("Observed region")
print(f"  total CDS:        {observed_total}")
print(f"  transposases/IS:  {observed_te}")
print()
print(f"Permutation test ({N_PERMUTATIONS} non-overlapping random windows, seed {SEED})")
print(f"  mean transposases (random):  {mean_random:.2f}")
print(f"  fold enrichment:             {fold:.2f}x")
print(f"  windows with >= observed:    {n_ge} / {N_PERMUTATIONS}")
print(f"  p-value (one-sided upper):   {pval:.4f}")

out_dist = work_path("statistics", "permutation_null_distribution.txt")
out_dist.parent.mkdir(parents=True, exist_ok=True)
out_dist.write_text("\n".join(str(x) for x in random_counts) + "\n")
print(f"\nWrote null distribution to: {out_dist}")

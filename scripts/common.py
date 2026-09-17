"""Shared helpers for the Python scripts: configuration, input tables and NCBI GFF3 parsing."""

import csv
import os
import re
import sys
from pathlib import Path
from urllib.parse import unquote


def cfg(name):
    """Return a setting from config/config.env (loaded into the environment)."""
    value = os.environ.get(name)
    if not value:
        sys.exit(f"ERROR: {name} is not set. From the repository root run:\n"
                 f"  set -a; source config/config.env; set +a")
    return value


def work_path(*parts):
    """Path under PROJECT_ROOT."""
    return Path(cfg("PROJECT_ROOT"), *parts)


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_region():
    """Return (sample, seqid, start, end) from the REGION file."""
    rows = read_tsv(cfg("REGION"))
    if len(rows) != 1:
        sys.exit(f"ERROR: {cfg('REGION')} must contain exactly one region")
    row = rows[0]
    if row["seqid"].startswith("<"):
        sys.exit(f"ERROR: fill in the chromosome accession in {cfg('REGION')}")
    return row["sample"], row["seqid"], int(row["start"]), int(row["end"])


def genome_row(sample):
    """Return the GENOMES row of a sample."""
    for row in read_tsv(cfg("GENOMES")):
        if row["sample"] == sample:
            return row
    sys.exit(f"ERROR: sample '{sample}' not found in {cfg('GENOMES')}")


def mobile_element_pattern():
    return re.compile(cfg("IS_PATTERN"), re.IGNORECASE)


def read_cds(gff_path):
    """Return the CDS features of an NCBI GFF3 file, one entry per feature ID.

    CDS split over several lines (joined or frameshifted) are merged into one
    feature spanning all segments. Pseudogenes are kept and flagged.

    The protein ID is the protein_id attribute (NCBI annotations) or, when there
    is none, the locus tag (Bakta annotations, whose protein FASTA headers are
    locus tags). Pseudogenes are flagged by pseudo=true (NCBI) or a pseudogene
    attribute (Bakta).
    """
    features = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            start, end = int(f[3]), int(f[4])
            fid = attrs.get("ID") or f"{f[0]}:{start}-{end}"
            if fid in features:
                feat = features[fid]
                feat["start"] = min(feat["start"], start)
                feat["end"] = max(feat["end"], end)
                continue
            features[fid] = {
                "seqid": f[0],
                "start": start,
                "end": end,
                "strand": f[6],
                "locus_tag": unquote(attrs.get("locus_tag", "")),
                "protein_id": unquote(attrs.get("protein_id") or attrs.get("locus_tag", "")),
                "product": unquote(attrs.get("product", "")),
                "pseudo": attrs.get("pseudo", "") == "true" or "pseudogene" in attrs,
            }
    return sorted(features.values(), key=lambda c: (c["seqid"], c["start"], c["end"]))


def in_window(cds, seqid, start, end):
    """A CDS belongs to a window when its midpoint lies inside it."""
    return cds["seqid"] == seqid and start <= (cds["start"] + cds["end"]) // 2 <= end


def sequence_lengths(fasta_path):
    """Return {sequence id: length} from a FASTA file."""
    lengths, current = {}, None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif current:
                lengths[current] += len(line)
    return lengths

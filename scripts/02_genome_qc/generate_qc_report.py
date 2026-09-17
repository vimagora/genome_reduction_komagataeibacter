#!/usr/bin/env python3
"""Generate an interactive HTML report of genome completeness + assembly QC.

Combines BUSCO completeness scores with assembly stats (# contigs, size, N50,
GC%) and taxonomy metadata into a single self-contained HTML file that opens
in any browser.

Usage (from the repository root, after loading config/config.env):
    python3 scripts/02_genome_qc/generate_qc_report.py \
        --busco-dir     "$PROJECT_ROOT/busco" \
        --assembly-dir  "$PROJECT_ROOT/genomes/fna" \
        --metadata      "$GENOMES" \
        --output        "$PROJECT_ROOT/busco/genome_qc_report.html"
"""
import argparse
import json
import re
import sys
from pathlib import Path
from datetime import datetime


def parse_short_summary(txt_path):
    """Return dict with C/S/D/F/M/n and lineage from a BUSCO short_summary file."""
    text = txt_path.read_text()
    m = re.search(
        r'C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%,n:(\d+)',
        text,
    )
    if not m:
        return None
    lin = re.search(r'lineage dataset(?: is)?:\s*(\S+)', text)
    return {
        "complete":   float(m.group(1)),
        "single":     float(m.group(2)),
        "duplicated": float(m.group(3)),
        "fragmented": float(m.group(4)),
        "missing":    float(m.group(5)),
        "n_markers":  int(m.group(6)),
        "lineage":    lin.group(1) if lin else "unknown",
    }


def compute_assembly_stats(fasta_path):
    """Return dict with n_contigs, total_length, longest, N50, gc_percent."""
    lengths, gc_total, n_total, cur_len, cur_gc = [], 0, 0, 0, 0
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur_len:
                    lengths.append(cur_len)
                    gc_total += cur_gc
                    n_total += cur_len
                cur_len = cur_gc = 0
            else:
                s = line.strip().upper()
                cur_len += len(s)
                cur_gc  += s.count("G") + s.count("C")
        if cur_len:
            lengths.append(cur_len)
            gc_total += cur_gc
            n_total  += cur_len
    if not lengths:
        return None
    lengths.sort(reverse=True)
    cum, n50 = 0, 0
    for L in lengths:
        cum += L
        if cum >= n_total / 2:
            n50 = L
            break
    return {
        "n_contigs":       len(lengths),
        "total_length":    n_total,
        "longest_contig":  lengths[0],
        "n50":             n50,
        "gc_percent":      round(100 * gc_total / n_total, 2),
    }


def parse_metadata(tsv_path):
    """Return dict keyed by sample name -> dict of columns."""
    out = {}
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            out[fields[0]] = dict(zip(header, fields))
    return out


def collect(busco_dir, assembly_dir, metadata_file):
    busco_dir    = Path(busco_dir)
    assembly_dir = Path(assembly_dir)
    metadata     = parse_metadata(metadata_file) if metadata_file else {}

    samples = []
    for sample_dir in sorted(p for p in busco_dir.iterdir() if p.is_dir()):
        sample = sample_dir.name

        # BUSCO 6.x names summaries like:
        #   short_summary.specific.<lineage>.<sample>.txt
        # Fall back to older patterns just in case.
        summaries = list(sample_dir.glob("short_summary.specific.*.txt")) \
                 or list(sample_dir.glob("short_summary*.txt"))
        if not summaries:
            print(f"WARN: no short summary for {sample}", file=sys.stderr)
            continue

        busco = parse_short_summary(summaries[0])
        if not busco:
            print(f"WARN: could not parse {summaries[0]}", file=sys.stderr)
            continue

        meta = metadata.get(sample, {})
        fna_name = meta.get("file_name", f"{sample}.fna")
        fna_path = assembly_dir / fna_name
        asm = compute_assembly_stats(fna_path) if fna_path.exists() else {}

        samples.append({
            "sample":   sample,
            "genus":    meta.get("genus", ""),
            "species":  meta.get("species", ""),
            "strain":   meta.get("strain", ""),
            "family":   meta.get("Family", ""),
            "busco":    busco,
            "assembly": asm or dict.fromkeys(
                ["n_contigs", "total_length", "longest_contig", "n50", "gc_percent"], None
            ),
        })
    return samples


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Genome QC Report</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
         margin: 2rem; color: #222; max-width: 1600px; }
  h1 { border-bottom: 2px solid #333; padding-bottom: 0.3rem; }
  h2 { margin-top: 2.5rem; }
  .meta { color: #666; font-size: 0.9rem; }
  .summary { display: flex; gap: 1rem; margin: 1rem 0 2rem; flex-wrap: wrap; }
  .stat { background: #f5f5f5; padding: 1rem 1.5rem; border-radius: 6px; min-width: 120px; }
  .stat .n { font-size: 1.6rem; font-weight: bold; }
  .stat .label { color: #666; font-size: 0.8rem; text-transform: uppercase; }
  table { border-collapse: collapse; width: 100%; margin-top: 1rem; font-size: 0.88rem; }
  th, td { padding: 0.4rem 0.7rem; text-align: left; border-bottom: 1px solid #eee; }
  th { background: #333; color: white; cursor: pointer; user-select: none;
       position: sticky; top: 0; }
  th:hover { background: #555; }
  tr:hover { background: #fafafa; }
  .flag-good { color: #2a8a2a; font-weight: 600; }
  .flag-warn { color: #b57500; font-weight: 600; }
  .flag-fail { color: #c22;    font-weight: 600; }
  .num { text-align: right; font-variant-numeric: tabular-nums; }
</style>
</head>
<body>
<h1>Genome QC Report</h1>
<p class="meta">Generated __GENERATED_AT__ &middot; BUSCO completeness + assembly statistics</p>

<div class="summary" id="summary"></div>

<h2>Completeness overview</h2>
<div id="chart"></div>

<h2>Per-genome details</h2>
<p class="meta">Click any column header to sort. Completeness thresholds:
  <span class="flag-good">&ge;95%</span>,
  <span class="flag-warn">90&ndash;95%</span>,
  <span class="flag-fail">&lt;90%</span>.
</p>
<table id="tbl">
<thead><tr>
  <th>Sample</th><th>Family</th>
  <th class="num">Complete %</th><th class="num">Single</th><th class="num">Duplicated</th>
  <th class="num">Fragmented</th><th class="num">Missing</th><th class="num">Markers</th>
  <th class="num">Contigs</th><th class="num">Size (Mb)</th>
  <th class="num">N50 (kb)</th><th class="num">GC %</th>
  <th>Lineage</th>
</tr></thead>
<tbody></tbody>
</table>

<script>
const data = __DATA_JSON__;

const fmt = (n, d=1) => (n == null) ? "&mdash;" : n.toFixed(d);
const fmtI = n => (n == null) ? "&mdash;" : n.toLocaleString();
const flag = c => c >= 95 ? "flag-good" : c >= 90 ? "flag-warn" : "flag-fail";

// Summary tiles
const good = data.filter(d => d.busco.complete >= 95).length;
const warn = data.filter(d => d.busco.complete >= 90 && d.busco.complete < 95).length;
const fail = data.filter(d => d.busco.complete < 90).length;
const med  = [...data].map(d => d.busco.complete).sort((a,b)=>a-b)[Math.floor(data.length/2)];
document.getElementById("summary").innerHTML = `
  <div class="stat"><div class="n">${data.length}</div><div class="label">Genomes</div></div>
  <div class="stat"><div class="n flag-good">${good}</div><div class="label">&ge;95% complete</div></div>
  <div class="stat"><div class="n flag-warn">${warn}</div><div class="label">90&ndash;95%</div></div>
  <div class="stat"><div class="n flag-fail">${fail}</div><div class="label">&lt;90%</div></div>
  <div class="stat"><div class="n">${fmt(med)}%</div><div class="label">Median completeness</div></div>
`;

// Stacked bar chart, sorted by completeness
const sorted = [...data].sort((a,b) => a.busco.complete - b.busco.complete);
const makeTrace = (key, name, color) => ({
    x: sorted.map(d => d.busco[key]),
    y: sorted.map(d => d.sample),
    name, type: "bar", orientation: "h",
    marker: {color},
    hovertemplate: `%{y}<br>${name}: %{x:.1f}%<extra></extra>`,
});
Plotly.newPlot("chart", [
    makeTrace("single",     "Single-copy", "#2a8a2a"),
    makeTrace("duplicated", "Duplicated",  "#4aa4de"),
    makeTrace("fragmented", "Fragmented",  "#e6a41c"),
    makeTrace("missing",    "Missing",     "#cc2222"),
], {
    barmode: "stack",
    height: Math.max(400, data.length * 22),
    margin: {l: 320, r: 40, t: 20, b: 40},
    xaxis: {title: "% of markers", range: [0, 100]},
    legend: {orientation: "h", y: -0.05},
}, {responsive: true});

// Sortable table
const tbody = document.querySelector("#tbl tbody");
function renderTable() {
    tbody.innerHTML = data.map(d => {
        const b = d.busco, a = d.assembly;
        return `<tr>
            <td>${d.sample}</td>
            <td>${d.family || ""}</td>
            <td class="num ${flag(b.complete)}">${fmt(b.complete)}</td>
            <td class="num">${fmt(b.single)}</td>
            <td class="num">${fmt(b.duplicated)}</td>
            <td class="num">${fmt(b.fragmented)}</td>
            <td class="num">${fmt(b.missing)}</td>
            <td class="num">${b.n_markers}</td>
            <td class="num">${fmtI(a.n_contigs)}</td>
            <td class="num">${a.total_length ? (a.total_length/1e6).toFixed(2) : "&mdash;"}</td>
            <td class="num">${a.n50 ? (a.n50/1e3).toFixed(1) : "&mdash;"}</td>
            <td class="num">${fmt(a.gc_percent, 2)}</td>
            <td>${b.lineage}</td>
        </tr>`;
    }).join("");
}
renderTable();

const cols = [
    d => d.sample, d => d.family || "",
    d => d.busco.complete, d => d.busco.single, d => d.busco.duplicated,
    d => d.busco.fragmented, d => d.busco.missing, d => d.busco.n_markers,
    d => d.assembly.n_contigs ?? -1, d => d.assembly.total_length ?? -1,
    d => d.assembly.n50 ?? -1, d => d.assembly.gc_percent ?? -1,
    d => d.busco.lineage,
];
let sortCol = null, sortAsc = true;
document.querySelectorAll("#tbl th").forEach((th, i) => {
    th.addEventListener("click", () => {
        if (sortCol === i) sortAsc = !sortAsc;
        else { sortCol = i; sortAsc = true; }
        data.sort((a, b) => {
            const va = cols[i](a), vb = cols[i](b);
            if (typeof va === "string") return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
            return sortAsc ? va - vb : vb - va;
        });
        renderTable();
    });
});
</script>
</body>
</html>
"""


def render_html(samples, generated_at):
    return (HTML_TEMPLATE
        .replace("__DATA_JSON__", json.dumps(samples))
        .replace("__GENERATED_AT__", generated_at))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--busco-dir",    required=True, help="parent dir with one BUSCO run per sample")
    p.add_argument("--assembly-dir", required=True, help="dir containing the input .fna files")
    p.add_argument("--metadata",     default=None,  help="optional metadata TSV")
    p.add_argument("--output",       default="genome_qc_report.html")
    args = p.parse_args()

    samples = collect(args.busco_dir, args.assembly_dir, args.metadata)
    if not samples:
        print("ERROR: no samples parsed. Check --busco-dir path.", file=sys.stderr)
        sys.exit(1)

    Path(args.output).write_text(
        render_html(samples, datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    )
    print(f"wrote {args.output} ({len(samples)} genomes)")


if __name__ == "__main__":
    main()